/*author David Edward Condon, University of Pennsylvania*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE//strtok_r, getline, strdup, mkstemp; this has to come before every #include
#endif
#include "safe_str_num.h"
#include <math.h>//lgammal, pow, sqrt
#include <limits.h>//UINT_MAX definition
#include <strings.h>//strcasecmp
#include <sys/stat.h> // mkdir
#include <unistd.h>//close
#ifdef _OPENMP
#include <omp.h>//parallelization, only when compiled with -fopenmp
#endif
#include "annotation.h"//for reading annotations and a few string functions
#include "fdr.h"
#include "dmr.h"//for calculating and identifying DMRs
#include "read_bisulfite_data.h"//reading bisulfite data
#include "my_string.h"//count_occurences
#include "option_range.h"//so I can get a range of options if I want to vary the options
#include "file_exists.h"
#include "is_program_installed.h"
#include "layout_R.h"
#include "uint.h"//adding and clearing 
#include "double.h"
#include "help.h"

struct options_struct {//will not contain all options
	unsigned int allowed_nucleotide_gap_between_CpN, MINIMUM_CPN_RANGE;
	unsigned int MAX_SIMILAR_CPN;
	bool make_bed_files;//, overwrite_output;
	unsigned short int CPU;
	unsigned short int ALLOWED_NON_DEFAULTS;
	bool make_figures, printing_CpN_list;//, isolation;
	bool print_all_nucleotides;
	bool print_DMR_pvalue;
	bool fdr_on_all_cpn;
	bool paired;
	bool make_heatmap;
};

void require_argument (const int argc, char *argv[], const int loop) {
//stop if the option at argv[loop] is the last thing on the command line, as reading argv[loop+1] would be reading NULL
	if (loop >= (argc-1)) {
		fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
		exit(EXIT_FAILURE);
	}
}

unsigned int unsigned_argument (const int argc, char *argv[], const int loop) {
//the unsigned integer that follows the option at argv[loop]
	require_argument(argc, argv, loop);
	if (is_string_integer(argv[loop+1]) == false) {
		fprintf(stderr, "%s '%s' isn't a non-negative integer.\n", argv[loop], argv[loop+1]);
		exit(EXIT_FAILURE);
	}
	const unsigned long int VALUE = safe_strtoul(argv[loop+1]);
	if (VALUE > UINT_MAX) {
		fprintf(stderr, "%s '%s' is too big, the maximum is %u.\n", argv[loop], argv[loop+1], UINT_MAX);
		exit(EXIT_FAILURE);
	}
	return (unsigned int)VALUE;
}

char * get_FDR_type (const char *restrict ARGUMENT, bool *restrict argument_used) {
//used by '-fdr' and '-v', which take an optional p-value adjustment method; ARGUMENT can be NULL
//returns a string to be freed, and sets argument_used if ARGUMENT named a method
	const char *restrict METHOD = "BH";//the default
	*argument_used = true;
	if (ARGUMENT == NULL) {
		*argument_used = false;
	} else if ((strcasecmp(ARGUMENT, "BH") == 0) || (strcasecmp(ARGUMENT, "fdr") == 0)) {
		METHOD = "BH";
	} else if (strcasecmp(ARGUMENT, "by") == 0) {
		METHOD = "BY";
	} else if (strcasecmp(ARGUMENT, "Bonferroni") == 0) {
		METHOD = "Bonferroni";
	} else if (strcasecmp(ARGUMENT, "hochberg") == 0) {
		METHOD = "Hochberg";
	} else if (strcasecmp(ARGUMENT, "holm") == 0) {
		METHOD = "Holm";
	} else if (strcasecmp(ARGUMENT, "hommel") == 0) {
		METHOD = "Hommel";
	} else {
		*argument_used = false;//the next argument is another option
	}
	char *restrict type = strdup(METHOD);
	if (type == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return type;
}

char * shell_quote (const char *restrict STRING) {
//wrap a filename in single quotes so that spaces etc. survive system(); returns a string to be freed
	char *restrict quoted = malloc(4 * strlen(STRING) + 3);//worst case: every character is a single quote, which becomes '\''
	if (quoted == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	size_t q = 0;
	quoted[q++] = '\'';
	for (size_t i = 0; STRING[i] != '\0'; i++) {
		if (STRING[i] == '\'') {
			memcpy(quoted + q, "'\\''", 4);
			q += 4;
		} else {
			quoted[q++] = STRING[i];
		}
	}
	quoted[q++] = '\'';
	quoted[q] = '\0';
	return quoted;
}

bool chromosome_in_every_replicate (const char *restrict CHROMOSOME, CHROMOSOME_CHANGES_STRUCT ***restrict chromosome_changes, unsigned int **restrict number_of_chromosomes, const unsigned short int *restrict COMPARED_SETS, const unsigned short int *restrict NO_OF_REPLICATES, size_t **restrict index_min, size_t **restrict index_max) {
//find where CHROMOSOME starts (index_min) and ends (index_max, exclusive) in every replicate's data, as files aren't necessarily sorted by chromosome in the same order
//each replicate has its own number of chromosomes, e.g. one file may lack chrY
//returns false if any replicate doesn't have this chromosome
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set]; r++) {
			const CHROMOSOME_CHANGES_STRUCT *restrict changes = chromosome_changes[COMPARED_SETS[set]][r];
			bool found = false;
			for (unsigned int chr_loop = 0; chr_loop < number_of_chromosomes[COMPARED_SETS[set]][r]; chr_loop++) {
				if (strcmp(CHROMOSOME, changes[chr_loop].chromosome) == 0) {
					index_min[set][r] = changes[chr_loop].index;
					index_max[set][r] = changes[chr_loop+1].index;//there is always one more element, "END"
					found = true;
					break;
				}
			}
			if (found == false) {
				return false;
			}
		}
	}
	return true;
}

void fprint_for_filename (FILE *restrict fh, const char *restrict STRING) {
//print gene names as part of a figure's filename inside an R script: nothing that ends the R string or makes a directory, and not too long for a filename
	for (size_t i = 0; (STRING[i] != '\0') && (i < 100); i++) {
		if ((STRING[i] == '/') || (STRING[i] == '\\') || (STRING[i] == '\'') || (STRING[i] == '"')) {
			fputc('_', fh);
		} else {
			fputc(STRING[i], fh);
		}
	}
}

bool strands_can_pair (const char PREVIOUS_STRAND, const char STRAND) {
//the two lines of one CpG are the C on the + strand, then the G (the - strand's C) at the next nucleotide.
//two adjacent lines on the same strand, or - followed by +, are two different cytosines.
//0 means the file type doesn't give a strand, then being adjacent is all that can be checked
	return (PREVIOUS_STRAND != '-') && (STRAND != '+');
}

#ifdef __GNUC__
__attribute__((returns_nonnull, malloc))//never returns NULL, which also keeps the compiler from warning about it
#endif
void * checked_malloc (const size_t SIZE, const unsigned int LINE) {
	void *restrict pointer = malloc(SIZE > 0 ? SIZE : 1);
	if (pointer == NULL) {
		fprintf(stderr, "malloc of %zu bytes failed @ %s line %u\n", SIZE, __FILE__, LINE);
		perror("");
		exit(EXIT_FAILURE);
	}
	return pointer;
}

int main (int argc, char *argv[]) {//"char *restrict argv[]" isn't a valid signature for main, and clang refuses to compile it
	puts("\n    ____  ___________________    _   ________");
   puts("   / __ \\/ ____/ ____/  _/   |  / | / /_  __/");
   puts("  / / / / __/ / /_   / // /| | /  |/ / / /   ");
   puts(" / /_/ / /___/ __/ _/ // ___ |/ /|  / / /    ");
   puts("/_____/_____/_/   /___/_/  |_/_/ |_/ /_/\n\nDifferential methylation: Easy, Fast, Identification and ANnoTation\n\n by David E. Condon, University of Pennsylvania, 2015-2020\n\n");
	if (argc == 1) {//in case options or files aren't entered
		fprintf(stderr, "The correct syntax is ./defiant [OPTIONS]... -i [FILES]\n\twhere [FILES] is a comma-separated list of input files within each replicate, and space-separated between different sets, e.g. control vs. case.\n");
		help();
		exit(EXIT_FAILURE);
	}
//See all default option settings here
//set options to defaults
	struct options_struct options;
	options.allowed_nucleotide_gap_between_CpN = 20000;
	options.MINIMUM_CPN_RANGE                  = 0;
	options.make_bed_files                     = false;
	options.CPU                                = 1;
	options.ALLOWED_NON_DEFAULTS               = 6;//there are 6 cutoffs; with 5, setting all 6 to non-default values gave 0 runs
	options.MAX_SIMILAR_CPN                    = 6;
//	options.overwrite_output                   = true;
	options.make_figures                       = false;
	options.printing_CpN_list                  = false;
//	options.isolation                          = false;
	options.print_all_nucleotides              = false;
	options.print_DMR_pvalue                   = false;
	options.fdr_on_all_cpn                     = false;
	options.paired                             = false;
	options.make_heatmap                       = false;
//get structs of option ranges, these are the defaults
	option_range_struct coverage = get_option_range("10");
	option_range_struct CpN = get_option_range("5");
	option_range_struct d = get_option_range("1");
	option_range_struct p = get_option_range("0.05");
	option_range_struct percent = get_option_range("0.1");//10%
	option_range_struct skip = get_option_range("0");
	cutoff_struct default_option_range;
	default_option_range.coverage	= 10;
	default_option_range.CpN		= 5;
	default_option_range.d			= 1;
	default_option_range.p			= 0.05;
	default_option_range.percent	= 0.10;
	default_option_range.skip		= 0;
//
	bool DEBUG = false;
//	bool dry = false;
	int input_file_start = 0;
	bool random_chromosomes_read = false, un_chromosomes_read = false;
	size_t labels_size = 0;
	unsigned int number_of_labels = 0;
	char *restrict total_label = NULL;
	char **restrict labels = NULL;
	char *restrict xaxis_label = strdup("CpN");
	char *restrict run_info_filename = strdup("defiant_run_info.txt");
	char *restrict FDR_type = NULL;
	for (unsigned short int loop = 1; loop < argc; loop++) {
		if (strcmp(argv[loop],"-b") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			options.make_bed_files = true;
			puts("You specified the '-b' option, so there will be an output bed file.");
		} else if (strcmp(argv[loop],"-c") == 0) {
			if (DEBUG == true) {
				printf("@ %s %u\n", __FILE__, __LINE__);
			}
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			coverage = get_option_range(argv[loop+1]);
			printf("The minimum coverage used will be %.0f\n",coverage.min);
			if (coverage.min != coverage.max) {
				printf("Coverage cutoffs will be incremented from %u to %u in steps of %u\n", (unsigned)coverage.min, (unsigned)coverage.max, (unsigned)coverage.step);
			}
			loop++;
		} else if (strcasecmp(argv[loop],"-CpN") == 0) {//command line flag for minimum CpN count
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			CpN = get_option_range(argv[loop+1]);
			if (CpN.min != CpN.max) {
				printf("CpN cutoffs will be incremented from %u to %u in steps of %u\n", (unsigned)CpN.min, (unsigned)CpN.max, (unsigned)CpN.step);
			}
			loop++;
		} else if (strcasecmp(argv[loop],"-cpu") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (is_string_integer(argv[loop+1]) == true) {
				options.CPU = (unsigned short)safe_strtoul(argv[loop+1]);
				printf("Calculations will be split among %u CPU.\n", options.CPU);
				loop++;
			} else {
				fprintf(stderr, "-cpu '%s' isn't an integer.\n", argv[loop+1]);
				exit(EXIT_FAILURE);
			}
		} else if (strcmp(argv[loop],"-d") == 0) {//looking for match to command line flag, 'd' = 'differentially methylated nucleotides
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			d = get_option_range(argv[loop+1]);
			if (d.min != d.max) {
				printf("You specified the '-d' option, so the minimum # of differentially methylated CpN will be increased from %u to %u in steps of %u\n", (unsigned)d.min, (unsigned)d.max, (unsigned)d.step);
			} else {
				printf("there will be a minimum of %u differentially methylated nucleotides in each region.\n",(unsigned)d.min);
			}
			loop++;
		} else if (strcasecmp(argv[loop], "-debug") == 0) {
			DEBUG = true;
/*		} else if (strcmp(argv[loop], "-dry") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			dry = true;*/
		} else if (strcmp(argv[loop],"-D") == 0) {//D for Defaults
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			require_argument(argc, argv, loop);
			if (is_string_integer(argv[loop+1]) == true) {
				options.ALLOWED_NON_DEFAULTS = (unsigned short)safe_strtoul(argv[loop+1]);
				if (options.ALLOWED_NON_DEFAULTS > 6) {
					printf("You specified %u allowed non-default cutoffs in parallelization, which is the equivalent of 6 non-defaults, as only 6 parameters can be parallelized.\n", options.ALLOWED_NON_DEFAULTS);
				}
				printf("You specified the '-D' option, so Defiant will calculate up to %u non-default cutoffs.\nThis will take Defiant much longer to compute.\n", options.ALLOWED_NON_DEFAULTS);
			} else {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "-D '%s' isn't an integer.\n", argv[loop+1]);
				exit(EXIT_FAILURE);
			}
			loop++;
		} else if (strcmp(argv[loop],"-E") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			options.print_all_nucleotides = true;//defined in dmr.h
			puts("option \"-E\" is set: A p-value and percent will be printed for each nucleotide.");
		} else if (strcmp(argv[loop],"-f") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			if (is_program_installed("Rscript") == true) {
				options.make_figures = true;
			} else {
				fprintf(stderr, "You specified '-f' to make figures, but Rscript is not installed.\n");
				exit(EXIT_FAILURE);
			}
			puts("You specified the '-f' option, so figures will be made for each and every DMR.");
		} else if (strcasecmp(argv[loop], "-fdr") == 0) {
			options.fdr_on_all_cpn = true;
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			bool method_given = false;
			free(FDR_type);//in case of a previous -fdr or -v
			FDR_type = get_FDR_type(loop < (argc-1) ? argv[loop+1] : NULL, &method_given);
			if (method_given == false) {
				loop--;
			}
			printf("Adjusted p-value will be calculated using %s's method on each CpN.\n", FDR_type);
			printf("I %c[1mSTRONGLY DISCOURAGE%c[0m this option, as it will\n", 27, 27);
			puts("give many false negatives and poor precision.\n");
			loop++;
		} else if (strcmp(argv[loop],"-G") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			if (loop == (argc-1)) {//prevents segfault on argv
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (is_string_integer(argv[loop+1]) == true) {
				options.allowed_nucleotide_gap_between_CpN = (unsigned int)safe_strtoul(argv[loop+1]);
			} else {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "-G %s isn't an integer.\n", argv[loop+1]);
				exit(EXIT_FAILURE);
			}
			printf("The maximum allowed gap between CpN will be %u nucleotides.\n", options.allowed_nucleotide_gap_between_CpN);
			loop++;
		} else if (strcasecmp(argv[loop],"-h") == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			free(xaxis_label); xaxis_label = NULL;
			free(run_info_filename); run_info_filename = NULL;
			help();
			exit(EXIT_SUCCESS);
		} else if (strcasecmp(argv[loop],"-heatmap") == 0) {
			if (is_program_installed("Rscript") == true) {
				options.make_heatmap = true;
			} else {
				fprintf(stderr, "You specified '-heatmap' to make figures, but Rscript is not installed.\n");
				exit(EXIT_FAILURE);
			}
		} else if (strcmp(argv[loop],"-i") == 0) {//start input files
//			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			input_file_start = loop+1;
			if ((argc - input_file_start) == 0) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "No data files were specified after '-i'.\n");
				free(run_info_filename); run_info_filename = NULL;
				exit(EXIT_FAILURE);
			}
			printf("There are %u sets to be read\n", argc - input_file_start);
			break;//no point in looping through input file names and making additional needless comparisons
/*		} else if (strcmp("-I", argv[loop]) == 0) {
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			options.isolation = true;
			puts("You specified the '-I' option, so the distance between *detected* CpG will be printed.");
		} else if (strcmp(argv[loop], "-j") == 0) {*/			
		} else if (strcmp(argv[loop],"-l") == 0) {//looking for match to command line flag
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			require_argument(argc, argv, loop);
			free(total_label);//in case -l is given twice
			total_label = strdup(argv[loop+1]);
			if (total_label == NULL) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "mallof of total_label failed @ %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			run_info_filename = realloc(run_info_filename, sizeof(char*) + strlen("run_info_.txt") + strlen(argv[loop+1]));
			if (run_info_filename == NULL) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "realloc of run_info_filename failed @ %s line %u\n", __FILE__, __LINE__);
				perror("");
				exit(EXIT_FAILURE);
			}
			strcpy(run_info_filename, "run_info_");
			strcat(run_info_filename, argv[loop+1]);
			strcat(run_info_filename, ".txt");//add null terminating string
			printf("Answer files will have \"%s\" in their names, and run information will be written to %s\n", argv[loop+1], run_info_filename);
			loop++;
		} else if (strcmp(argv[loop],"-L") == 0) {//labels with -f option
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
//			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			char *restrict label = NULL;
			char * labels_pointer = NULL;
			label = strtok_r(argv[loop+1],",",&labels_pointer);
			while (label != NULL) {//an empty string, or one of only commas, has no labels
				labels = add_name(labels, &labels_size, &number_of_labels, label, false, true);
				label = strtok_r(NULL,",", &labels_pointer);
			}
			for (unsigned short int x = 0; x < number_of_labels; x++) {
				printf("Set %u will be labeled \"%s\"\n", x+1, labels[x]);
			}
			loop++;
		} else if (strcmp(argv[loop],"-N") == 0) {//print nucleotide list (Yu-Chin's request)
			if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}
			options.printing_CpN_list = true;
			puts("\nYou specified the '-N' option, so Defiant will print a list of individually differentially methylated nucleotides in the answers file.\n");
/*		} else if (strcmp(argv[loop],"-o") == 0) {
			options.overwrite_output = false;
			puts("You specified the '-o' option, so files will not be overwritten if they exist.");*/
		} else if (strcmp(argv[loop],"-p") == 0) {//looking for match to command line flag
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			p = get_option_range(argv[loop+1]);
			if ((p.min < 0.0) || (p.min > 1.0) || (p.max < 0.0) || (p.max > 1.0)) {
				fprintf(stderr, "Maximum p-value is defined to be 0 <= p <= 1, and you selected %s\n",argv[loop+1]);
				exit(EXIT_FAILURE);
			}
			if (p.min != p.max) {
				printf("p-value cutoffs will be incremented from %lf to %lf in steps of %lf\n", p.min, p.max, p.step);
			}
			loop++;//I already know what the next value will be
		} else if (strcmp(argv[loop],"-P") == 0) {//command line flag for minimum percent change
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			percent = get_option_range(argv[loop+1]);
			if ((percent.min < 0.0) || (percent.min > 100.0) || (percent.max < 0.0) || (percent.max > 100.0)|| (percent.step < 0.0) || (percent.step > 100.0)) {
				fprintf(stderr, "Minimum percent is defined to be 0 <= p <= 100, and you selected %s\n",argv[loop+1]);
				exit(EXIT_FAILURE);
			}
			if (percent.min != percent.max) {
				printf("Percentage cutoffs will be incremented from %lf to %lf in steps of %lf\n", percent.min, percent.max, percent.step);
			}
			percent.min /= 100;
			percent.max /= 100;
			percent.step /= 100;
			loop++;
		} else if (strcasecmp(argv[loop], "--paired") == 0) {
			options.paired = true;
			puts("You specified '--paired', so p-values will be paired.");
		} else if ((strcmp(argv[loop],"-q") == 0) || (strcmp(argv[loop],"-promoter_cutoff") == 0)) {
			promoter_cutoff = unsigned_argument(argc, argv, loop);
			printf("Intergenic genes will be assigned a promoter cutoff of %u nucleotides for identification.\n", promoter_cutoff);
			loop++;
		} else if (strcmp(argv[loop],"-r") == 0) {//command line flag for CpN range
			options.MINIMUM_CPN_RANGE = unsigned_argument(argc, argv, loop);//was cast to unsigned short, so 100000 became 34464
			printf("The minimum CpN range is set to %u nucleotides.\n",options.MINIMUM_CPN_RANGE);
			loop++;
		} else if (strcmp(argv[loop],"-R") == 0) {//command line flag for reading random chromosomes
			random_chromosomes_read = true;
			puts("Random nucleotides will be included in analysis.\n");
		} else if (strcmp(argv[loop],"-S") == 0) {
			require_argument(argc, argv, loop);
			if (strcasecmp(argv[loop+1],"inf") == 0) {//case insensitive comparison
				skip.min = UINT_MAX-1;
				skip.max = UINT_MAX-1;
				skip.step = 1.0;
				loop++;
				printf("You entered '-S %s', so the limit of skips will be infinite (set to %.0f, which is big enough)\n", argv[loop], skip.min);
			} else {
				skip = get_option_range(argv[loop+1]);
				if (skip.min != skip.max) {
					printf("Skip cutoffs will be incremented from %.0f to %.0f in steps of %lf\n", skip.min, skip.max, skip.step);
				} else {
					printf("There will be %.0f allowed consecutive skip(s) of low coverage nucleotides when reading data.\n",skip.min);
				}
				loop++;
			}
		} else if (strcmp(argv[loop],"-s") == 0) {
			if (loop == (argc-1)) {//prevents segfault on argv
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (is_string_integer(argv[loop+1]) == true) {
				options.MAX_SIMILAR_CPN = (unsigned)safe_strtoul(argv[loop+1]);
				printf("The maximum number of allowed similar CpN after differentially methylated CpN is %u.\n", options.MAX_SIMILAR_CPN);
				loop++;
			} else {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "-s '%s' isn't an integer.\n", argv[loop+1]);
				free(run_info_filename); run_info_filename = NULL;
				exit(EXIT_FAILURE);
			}
		} else if (strcmp(argv[loop],"-U") == 0) {
			un_chromosomes_read = true;
			puts("Un nucleotides will be read.");
		} else if (strcmp(argv[loop],"-v") == 0) {
			options.print_DMR_pvalue = true;
			bool method_given = false;
			free(FDR_type);//in case of a previous -fdr or -v
			FDR_type = get_FDR_type(loop < (argc-1) ? argv[loop+1] : NULL, &method_given);
			if (method_given == false) {
				loop--;
			}
			loop++;
			printf("You specified the '-v' option, so a p-value & an adjusted p-value using %s's method will be printed for each DMR.\n", FDR_type);
/*		} else if (strcmp(argv[loop],"-w") == 0) {
			options.pvalue_scaling = true;
			puts("you specified the '-w' option, so DMR p-values will be scaled according to the distribution.");*/
		} else if (strcmp(argv[loop], "-version") == 0) {
			free(xaxis_label); xaxis_label = NULL;
			free(run_info_filename); run_info_filename = NULL;
			version();
			exit(EXIT_SUCCESS);
		} else if (strcmp(argv[loop],"-x") == 0) {
			if (loop == (argc-1)) {//prevents segfault on argv
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (is_program_installed("Rscript") == false) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "You are attempting to make figures, but you don't have Rscript installed.\n");
				exit(EXIT_FAILURE);
			}
			free(xaxis_label); xaxis_label = NULL;
			xaxis_label = strdup(argv[loop+1]);
			options.make_figures = true;
			printf("You specified the \"-x\" option, so the xaxis in figures will be labeled \"%s\".  Figures will be made.\n", xaxis_label);
			loop++;
		} else if (strcmp(argv[loop],"-a") == 0) {//command line flag for annotation
			if (DEBUG == true) {
				printf("@ %s %u\n", __FILE__, __LINE__);
			}
			if (loop == (argc-1)) {//prevents segfault on argv
				fprintf(stderr, "the '%s' option requires something to come after it.\n", argv[loop]);
				fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			FILE *restrict file;//does this file exist? 
			if ((file = fopen(argv[loop+1], "r")) == NULL) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				fprintf(stderr, "gene annotation file %s (specified by -a flag) could not be read.\n",argv[loop+1]);
				perror("");
				exit(EXIT_FAILURE);
			}
			fclose(file);
//---------------------------------------------------------------------------------------------
//Read gene data to memory
//---------------------------------------------------------------------------------------------
//check file type
			unsigned short int annotation_file_type = 0;		
			FILE *restrict fh0;
			fh0 = fopen(argv[loop+1],"r");
			if (fh0 == NULL) {
				free(xaxis_label); xaxis_label = NULL;
				free(run_info_filename); run_info_filename = NULL;
				perror("");
				fprintf(stderr, "Failed to read %s\n", argv[loop+1]);
				exit(EXIT_FAILURE);
			}
	   	size_t len = 0;
	   	ssize_t read = 0;
			char *gene_line = NULL;
			bool file_read = false;
	   	while ((read = getline(&gene_line, &len, fh0)) != -1) {//http://man7.org/linux/man-pages/man3/getline.3.html
	   		file_read = true;
	   		annotation_file_type = get_annotation_type(gene_line);
	   		if (annotation_file_type > 0) {
//	   			printf("annotation_file_type = %u at %s line %u\n", annotation_file_type, __FILE__, __LINE__);
	   			break;
	   		}//keep going if the annotation file isn't defined (i.e. = 0)
			}
			free(gene_line); gene_line = NULL;
			fclose(fh0);
			if (file_read == false) {
				fprintf(stderr, "%s wasn't read.  You may have entered a file of 0 size.\n", argv[loop+1]);
				fprintf(stderr, "Failed @ %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			if (annotation_file_type == 0) {
				printf("Failed to find annotation type for %s @ %s line %u\n\n",argv[loop+1], __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
//
			char *restrict annotation_filename = NULL;
			bool sorted = true;
			if (annotation_file_type != 1) {//only type 1 (refFlat) has to be sorted first
				annotation_filename = strdup(argv[loop+1]);//if the file doesn't need sorting
				sorted = false;
			} else {//i.e., the annotation file needs to be sorted
//a temporary file with a unique name, so that two Defiant jobs running at once don't write over each other's sorted annotation
				const char *restrict TMPDIR = getenv("TMPDIR");
				if ((TMPDIR == NULL) || (TMPDIR[0] == '\0')) {
					TMPDIR = "/tmp";
				}
				annotation_filename = checked_malloc(strlen(TMPDIR) + strlen("/defiant_annotation_XXXXXX") + 1, __LINE__);
				sprintf(annotation_filename, "%s/defiant_annotation_XXXXXX", TMPDIR);
				const int TEMPORARY_FILE = mkstemp(annotation_filename);
				if (TEMPORARY_FILE == -1) {
					fprintf(stderr, "Couldn't make a temporary file in %s\n", TMPDIR);
					perror("");
					exit(EXIT_FAILURE);
				}
				close(TEMPORARY_FILE);
				char *restrict quoted_input = shell_quote(argv[loop+1]), *restrict quoted_output = shell_quote(annotation_filename);
				char *restrict command = checked_malloc(strlen(quoted_input) + strlen(quoted_output) + 200, __LINE__);//200 is about 2x the length of the grep command
				sprintf(command,"grep -v chrUn %s | grep -v random | sort -k 3,3 -k 5n,5n -k 6n,6n -k 7n,7n -k 8n,8n > %s", quoted_input, quoted_output);
				if (system(command) != 0) {
					fprintf(stderr, "sort command failed.\n");
					remove(annotation_filename);
					exit(EXIT_FAILURE);
				}
				free(command); command = NULL;
				free(quoted_input); quoted_input = NULL;
				free(quoted_output); quoted_output = NULL;
			}
			FILE *restrict fh = fopen(annotation_filename,"r");
			if (fh == NULL) {
				fprintf(stderr, "Could not open %s :( dying :,(\n",annotation_filename);
				perror("");
				exit(EXIT_FAILURE);
			}
			printf("Reading gene annotation %s\n", argv[loop+1]);
			gene_annotation_defined = true;
			char *restrict previous_chromosome = strdup("");//no chromosome has this name
//			printf("annotation_file_type = %u\n", annotation_file_type);
	   	GENES_chromosome_change_index = 0;//globally accessed variable declared before main
	   	size_t GENES_size = 0;
	   	while ((read = getline(&gene_line, &len, fh)) != -1) {//http://man7.org/linux/man-pages/man3/getline.3.html
				ANNOTATION_LINE_STRUCT line = read_gene_annotation_line(gene_line, annotation_file_type);
//				printf("line.sense = %d\n", line.sense);
//				printf("line.gene_start = %u\n", line.gene_start);
//				printf("gene_name = %s\n", line.gene_name);
				if (line.sense == 0) {//i.e. this is a transcript and not a gene
					free(line.chromosome); line.chromosome = NULL;
					free(line.gene_name);  line.gene_name  = NULL;
					continue;
				}//don't allocate if it's not necessary
//				printf("@%s line %u	%s\n", __FILE__, __LINE__, gene_line);//debugging
//				printf("gene name %s	start %u	end %u	sense %d	chr %s\n", line.gene_name, line.gene_start, line.gene_end, line.sense, line.chromosome);
				if ((number_of_annotated_genes > 0) && (strcmp(previous_chromosome, line.chromosome) == 0) && (strcmp(line.gene_name, GENES[number_of_annotated_genes-1].gene_name) == 0)) {//same gene name on the same chromosome
//				printf("%s %u-%u\n",GENES[number_of_annotated_genes-1].gene_name, GENES[number_of_annotated_genes-1].gene_start, GENES[number_of_annotated_genes-1].gene_end);
					GENES[number_of_annotated_genes-1].gene_start = line.gene_start < GENES[number_of_annotated_genes-1].gene_start ? line.gene_start : GENES[number_of_annotated_genes-1].gene_start;
					GENES[number_of_annotated_genes-1].gene_end = line.gene_end > GENES[number_of_annotated_genes-1].gene_end ? line.gene_end : GENES[number_of_annotated_genes-1].gene_end;
					free(line.chromosome); line.chromosome = NULL;
					free(line.gene_name);  line.gene_name  = NULL;
//					printf("%s %u-%u\n",GENES[number_of_annotated_genes-1].gene_name, GENES[number_of_annotated_genes-1].gene_start, GENES[number_of_annotated_genes-1].gene_end);
					continue;//some GTF files lack genes, have only transcripts, so exact gene boundaries aren't known
				}
				GENES_size += sizeof(GENES_STRUCT) + strlen(line.gene_name);
				GENES = realloc(GENES, GENES_size);//add memory for this gene
				if (GENES == NULL) {//error checking
					fprintf(stderr, "realloc of GENES failed at index %u\n",number_of_annotated_genes);
					fprintf(stderr, "GENES struct failed realloc for line %s; file %s line %u.\n",gene_line,__FILE__,__LINE__);//DEBUGGING
					perror("");
					exit(EXIT_FAILURE);
				}
				GENES[number_of_annotated_genes].gene_start = line.gene_start;
				GENES[number_of_annotated_genes].gene_end = line.gene_end;
				GENES[number_of_annotated_genes].sense = line.sense;
				GENES[number_of_annotated_genes].gene_name = strdup(line.gene_name);
				if (strcmp(previous_chromosome,line.chromosome) != 0) {
					number_of_chromosomes_in_GENES++;
					gene_chromosome_changes = realloc(gene_chromosome_changes, (size_t)number_of_chromosomes_in_GENES*sizeof(CHROMOSOME_CHANGES_STRUCT));
					if (gene_chromosome_changes == NULL) {
						printf("@line %u	realloc of gene_chromosome_changes failed.\n",__LINE__);
						exit(EXIT_FAILURE);
					}
					gene_chromosome_changes[GENES_chromosome_change_index].index = number_of_annotated_genes;
					gene_chromosome_changes[GENES_chromosome_change_index].chromosome = strdup(line.chromosome);
					if (DEBUG == true) {
						printf("chromosome %s, gene_chromosome_changes[%u] = %u\n", gene_chromosome_changes[GENES_chromosome_change_index].chromosome, GENES_chromosome_change_index, gene_chromosome_changes[GENES_chromosome_change_index].index);
					}
					GENES_chromosome_change_index++;
					if (DEBUG == true) {
						printf("@ %s %u\n", __FILE__, __LINE__);
					}
					free(previous_chromosome); previous_chromosome = NULL;
					previous_chromosome = strdup(line.chromosome);
					if (DEBUG == true) {
						printf("@ %s %u\n", __FILE__, __LINE__);
					}
					if (strcmp(previous_chromosome, line.chromosome) != 0) {
						fprintf(stderr, "%s failed copy from %s @ line %u\n",previous_chromosome,line.chromosome,__LINE__);
						exit(EXIT_FAILURE);
					}
				}
				free(line.chromosome); line.chromosome = NULL;
				free(line.gene_name);  line.gene_name  = NULL;
				number_of_annotated_genes++;//also by line
			}//end reading file
			fclose(fh);
			free(previous_chromosome); previous_chromosome = NULL;
			free(gene_line); gene_line = NULL;
			if (sorted == true) {//sorted_gene_filename will be хрень in this case
				if (remove(annotation_filename) != 0) {//deleting temporary sorted file
					fprintf(stderr, "Couldn't delete %s\n",annotation_filename);
					exit(EXIT_FAILURE);
				}
			}
			free(annotation_filename); annotation_filename = NULL;
//			printf("GENES[%u]:\tChromosome %s end -> [%u,%u]\n",index,gene_chromosome_changes[GENES_chromosome_change_index].chromosome,GENES[index].gene_start,GENES[index].gene_end);//debugging
			printf("There are %u genes and %u chromosomes read from %s.\n",number_of_annotated_genes, number_of_chromosomes_in_GENES, argv[loop+1]);
			if ((number_of_annotated_genes == 0) && (number_of_chromosomes_in_GENES == 0)) {
				gene_annotation_defined = false;//this prevents segfaults if you search for genes that don't exist
			}
//			printf("finished reading gene data from %s into RAM.\n",argv[loop+1]);//Tell user program has finished reading gene data into RAM
			if (DEBUG == true) {
				if (number_of_annotated_genes >= 1) {
					printf("last gene goes from %u-%u\n", GENES[number_of_annotated_genes-1].gene_start, GENES[number_of_annotated_genes-1].gene_end);
				} else {
					puts("There are 0 genes in GENES struct.");
				}
			}
			loop++;//the next string will be the file, no sense in strcmp the file
		} else {
			free(xaxis_label); xaxis_label = NULL;
			free(run_info_filename); run_info_filename = NULL;
			printf("\"%s\" isn't a recognized option.\n", argv[loop]);
//			help();
			exit(EXIT_FAILURE);
		}
	}
	check_uint_option_range(coverage, "-c");
	check_uint_option_range(CpN, "-CpN");
	check_uint_option_range(d, "-d");
	check_uint_option_range(skip, "-S");
	if (input_file_start == 0) {
		free(xaxis_label); xaxis_label = NULL;
		free(run_info_filename); run_info_filename = NULL;
		fprintf(stderr, "Input files were not defined with the '-i' option.\nThe correct syntax is `defiant [OPTIONS]... -i [FILES]'\n\twhere [FILES] is a comma-separated list of input files within each replicate, and space-separated between different sets, e.g. control vs. case.\n");
//		help();
		exit(EXIT_FAILURE);
	}
	const unsigned short int SETS = (unsigned short)(argc-input_file_start);
/*	if ((SETS > 2) && (FDR_type != NULL)) {
		fprintf(stderr, ">= 2 input sets cannot be compared with the -v option.\nPlease run the comparisons two at a time with the -v option.\n");
		exit(EXIT_FAILURE);
	}*/
	if (input_file_start > (argc-2)) {
		free(xaxis_label); xaxis_label = NULL;
		free(run_info_filename); run_info_filename = NULL;
		help();
		fprintf(stderr, "\nThere must be at least two sets of data to be compared.\n");
		exit(EXIT_FAILURE);
	}
/*	if (options.paired == true) {
		if (input_file_start != (argc - 2)) {
			free(xaxis_label); xaxis_label = NULL;
			free(run_info_filename); run_info_filename = NULL;
			fprintf(stderr, "\n\nIn a paired analysis, there can only be 2 sets compared.\n\n");
			exit(EXIT_FAILURE);
		}
	}
	if (DEBUG == true) {printf("@ %s %u\n", __FILE__, __LINE__);}*/
//------------------------------------------
//make struct of cutoffs
//-----------------------------------------
	const unsigned int RUNS = how_many_runs(coverage, CpN, d, p, percent, skip, default_option_range, options.ALLOWED_NON_DEFAULTS, options.CPU);
	printf("There will be %u runs.\n", RUNS);
	if (RUNS == 0) {
		fprintf(stderr, "No combination of cutoffs is within the %u non-default cutoffs allowed by '-D', so there is nothing to run.\n", options.ALLOWED_NON_DEFAULTS);
		exit(EXIT_FAILURE);
	}
/*	if (dry == true) {
		free(xaxis_label); xaxis_label = NULL;
		free(run_info_filename); run_info_filename = NULL;
		printf("c\tCpN\td\tp\tP\tS\n");
		for (unsigned int run = 0; run < RUNS; run++) {
			printf( "%u\t%u\t%u\t%g\t%.1f\t%u\n", cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p, 100.0*cutoff[run].percent, cutoff[run].skip);
		}
		free(cutoff); cutoff = NULL;
		exit(EXIT_SUCCESS);
	}*/
	if (options.make_figures == true) {
		if (file_exists("DMR_figures") == false) {
			if (mkdir("DMR_figures",0777) == -1) {
				perror("Cannot create directory DMR_figures");
				printf("Failing at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		if (RUNS > 1) {
			printf("You are making figures for %u runs, this will take a lot longer.\n", RUNS);
		}
	} else {
		free(xaxis_label); xaxis_label = NULL;
	}
//
	FILE *restrict run_info;
	printf("Run information will be in %s\n", run_info_filename);
	run_info = fopen(run_info_filename,"w");
	if (run_info == NULL) {
		fprintf(stderr, "Failed to write %s\n", run_info_filename);
		perror("");
		exit(EXIT_FAILURE);
	}
	fprintf(run_info,"The command line: ");
	if (DEBUG == true) {
		printf("at %s line %u\n", __FILE__, __LINE__);
	}
	for (unsigned short int x = 0; x < argc; x++) {//print the command line arguments to the run_info file
		fprintf(run_info," %s",argv[x]);
	}
	fprintf(run_info,"\n");
	if (SETS != number_of_labels) {
		puts("\n***WARNING***");
		printf("There is/are %u label(s), yet %u sets, just a warning that you may not have correctly labeled data.\n", number_of_labels, SETS);
		puts("***WARNING***\n");
	}
	unsigned short int *restrict number_of_replicates_per_set = checked_malloc(SETS*sizeof(unsigned short int), __LINE__);
	CHROMOSOME_CHANGES_STRUCT ***restrict bs_data_chromosome_changes = checked_malloc(SETS*sizeof(CHROMOSOME_CHANGES_STRUCT**), __LINE__);//trace chromosome changes to accelerate searching later on, chromosome_changes[set][replicate][chromosome_loop]
	unsigned int **restrict number_of_data_chromosomes = checked_malloc(SETS*sizeof(unsigned int*), __LINE__);//number_of_data_chromosomes[set][replicate], as each file can have a different number of chromosomes
	BISULFITE_DATA_STRUCT ***restrict BISULFITE_DATA = checked_malloc(SETS*sizeof(BISULFITE_DATA_STRUCT**), __LINE__);//where bisulfite sequencing data goes, BISULFITE_DATA[set][replicate][index]
//---------------------------------------------------------------------------------------------
//Begin reading data files
//---------------------------------------------------------------------------------------------
	size_t **restrict global_methylation = checked_malloc(SETS*sizeof(size_t*), __LINE__);//global coverage can have >= 1 set, so I cannot store values in DMR_struct, which can only have 2 sets
	size_t **restrict global_coverage = checked_malloc(SETS*sizeof(size_t*), __LINE__);
	for (int set_loop = input_file_start; set_loop < argc; set_loop++) {//for each set
		const unsigned short int set = (short unsigned)(set_loop-input_file_start);
//nucleotides below minimum coverage can only be thrown away while reading if no run will need them
		unsigned int MINIMUM_READ_COVERAGE = (unsigned)skip.max > 0 ? 0 : (unsigned)coverage.min;
		if (options.print_all_nucleotides == true) {
			MINIMUM_READ_COVERAGE = 0;
		}
		number_of_replicates_per_set[set] = count_tokens(argv[set_loop], ',');//the number of files that strtok_r will find below; "a,b," has 2 files, not 3
		if (number_of_replicates_per_set[set] == 0) {
			fprintf(stderr, "Set %u, \"%s\", doesn't have any files.\n", set+1, argv[set_loop]);
			exit(EXIT_FAILURE);
		}
		if (DEBUG == true) {
			printf("There are %u replicates in set %u @ line %u.\n",number_of_replicates_per_set[set],set, __LINE__);
		}
		global_methylation[set]         = calloc(number_of_replicates_per_set[set], sizeof(size_t));
		global_coverage[set]            = calloc(number_of_replicates_per_set[set], sizeof(size_t));
		BISULFITE_DATA[set]             = calloc(number_of_replicates_per_set[set], sizeof(BISULFITE_DATA_STRUCT*));
		bs_data_chromosome_changes[set] = calloc(number_of_replicates_per_set[set], sizeof(CHROMOSOME_CHANGES_STRUCT*));
		number_of_data_chromosomes[set] = calloc(number_of_replicates_per_set[set], sizeof(unsigned int));
		if ((global_methylation[set] == NULL) || (global_coverage[set] == NULL) || (BISULFITE_DATA[set] == NULL) || (bs_data_chromosome_changes[set] == NULL) || (number_of_data_chromosomes[set] == NULL)) {
			fprintf(stderr, "calloc failed for set %u @ %s line %u\n", set+1, __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		fprintf(run_info,"These files are considered set %u: %s\n\n",set+1,argv[set_loop]);
		char *restrict files;
		char *files_pointer;
		files = strtok_r(argv[set_loop],",",&files_pointer);//start string token extraction
		unsigned short int replicate = 0;
		while (files != NULL) {//parsing arguments from command line, for each file
			char *restrict previous_chromosome = strdup("");//the chromosome of the last nucleotide that was stored
			FILE *restrict BISULFITE_DATA_file;
			BISULFITE_DATA_file = fopen(files,"r");
			if (BISULFITE_DATA_file == NULL) {
				fprintf(stderr, "Could not open %s\n",files);
				perror("");
				exit(EXIT_FAILURE);
			}
			printf("Now reading file %s into memory...\n",files);//debugging, telling end user what's happening
			char *bs_line = NULL;
			size_t len = 0;
			unsigned int index = 0;//the number of nucleotides stored from this file
			ssize_t read;//for getline function
//the arrays double in size when they're full, and are trimmed to their final size after the file is read
			size_t data_capacity = 4096;
			BISULFITE_DATA_STRUCT *restrict data = checked_malloc(data_capacity * sizeof(BISULFITE_DATA_STRUCT), __LINE__);
			size_t changes_capacity = 64;
			CHROMOSOME_CHANGES_STRUCT *restrict changes = checked_malloc(changes_capacity * sizeof(CHROMOSOME_CHANGES_STRUCT), __LINE__);
			unsigned int chromosome_count = 0;
//about getting file types
			unsigned int failed_coverage = 0, passed_coverage = 0;
			bool two_strands = false;
			bool last_nucleotide_is_open = false;//two_strands only: the last stored nucleotide may still get its 2nd strand, so its coverage isn't final
			char previous_strand = 0;
			FILE_TYPE_STRUCT file_type;
			file_type.type = 0;
			file_type.header = false;
			while ((read = getline ( &bs_line, &len, BISULFITE_DATA_file)) != -1) {//splitting files based on whitespace
				if (file_type.type == 0) {//i.e. not defined
	 				file_type = get_file_type(bs_line);
	 				if (DEBUG == true) {
		 				printf("file_type = %u\n", file_type.type);
		 			}
	 				if ((file_type.type == 2) || (file_type.type == 5) || (file_type.type == 8) || (file_type.type == 9)) {
	 					two_strands = true;
	 				}
	 			}
	 			if (file_type.header == true) {
	 				file_type.header = false;
	 				continue;
	 			}
	 			LINE_STRUCT line = read_data_line(bs_line, file_type.type);
				if (line.nucleotide == 0) {
					free(line.chromosome); line.chromosome = NULL;
					continue;
				}
//"-R" keeps chromosomes like chr1_gl000191_random, "-U" keeps chromosomes like chrUn_gl000211; only the chromosome's name is searched, not the whole line
				if ((random_chromosomes_read == false) && (strstr(line.chromosome,"_random") != NULL)) {
					free(line.chromosome); line.chromosome = NULL;
					continue;
				}
				if ((un_chromosomes_read == false) && (strstr(line.chromosome,"Un") != NULL)) {
					free(line.chromosome); line.chromosome = NULL;
					continue;
				}
				if (strstr(line.chromosome,"_alt") != NULL) {//skip chromosomes with "_alt"
					free(line.chromosome); line.chromosome = NULL;
					continue;
				}
		 		const bool SAME_CHROMOSOME = (chromosome_count > 0) && (strcmp(line.chromosome,previous_chromosome) == 0);
				if (two_strands == false) {
					if (line.coverage < MINIMUM_READ_COVERAGE) {
						failed_coverage++;
						free(line.chromosome); line.chromosome = NULL;
						continue;
					}
				} else if ((last_nucleotide_is_open == true) && (SAME_CHROMOSOME == true) && (line.nucleotide == (data[index-1].nucleotide+1)) && (strands_can_pair(previous_strand, line.strand) == true)) {
					data[index-1].methylated_C += line.methylated_C;//the 2nd strand of the same CpG is added to the 1st
					data[index-1].coverage += line.coverage;
					previous_strand = line.strand;
					free(line.chromosome); line.chromosome = NULL;
					continue;
				}
				if (last_nucleotide_is_open == true) {//the last nucleotide can't get any more reads, so now its coverage is known
					if (data[index-1].coverage < MINIMUM_READ_COVERAGE) {
						failed_coverage++;
						index--;//this nucleotide will be written over
					} else {
						passed_coverage++;
						global_methylation[set][replicate] += data[index-1].methylated_C;
						global_coverage[set][replicate] += data[index-1].coverage;
					}
					last_nucleotide_is_open = false;
				}
				if (SAME_CHROMOSOME == false) {
					if (DEBUG == true) {
						printf("found a switch from chromosome %s to %s for element %u\n",previous_chromosome, line.chromosome, chromosome_count);
					}
					if ((chromosome_count + 2) > changes_capacity) {//+2 because of the "END" element
						changes_capacity *= 2;
						changes = realloc(changes, changes_capacity * sizeof(CHROMOSOME_CHANGES_STRUCT));
						if (changes == NULL) {
							fprintf(stderr, "@line %u	realloc of bs_data_chromosome_changes failed @ chromosome %s.\n", __LINE__, line.chromosome);
							exit(EXIT_FAILURE);
						}
					}
					changes[chromosome_count].index = index;
					changes[chromosome_count].chromosome = strdup(line.chromosome);
					if (DEBUG == true) {
						printf("%s -> bs_data_chromosome_changes[%u][%u][%u].index = %u (at %s line %u)\n",line.chromosome,set,replicate,chromosome_count, index, __FILE__, __LINE__);//debugging
					}
					free(previous_chromosome); previous_chromosome = NULL;
					previous_chromosome = strdup(line.chromosome);
					if ((changes[chromosome_count].chromosome == NULL) || (previous_chromosome == NULL)) {
						fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
						exit(EXIT_FAILURE);
					}
					chromosome_count++;
				}
				free(line.chromosome); line.chromosome = NULL;
				if (index == UINT_MAX) {
					fprintf(stderr, "%s has more than %u nucleotides, which is more than Defiant can index.\n", files, UINT_MAX);
					exit(EXIT_FAILURE);
				}
				if (index == data_capacity) {
					data_capacity *= 2;
					data = realloc(data, data_capacity * sizeof(BISULFITE_DATA_STRUCT));
					if (data == NULL) {
						fprintf(stderr, "realloc of BISULFITE_DATA[%u][%u] failed @ line %u.",set ,replicate, __LINE__);
						perror("");
						exit(EXIT_FAILURE);
					}
				}
				data[index].nucleotide = line.nucleotide;
				data[index].methylated_C = line.methylated_C;
				data[index].coverage = line.coverage;
				index++;
				previous_strand = line.strand;
				if (two_strands == true) {
					last_nucleotide_is_open = true;
				} else {
					passed_coverage++;
					global_methylation[set][replicate] += line.methylated_C;
					global_coverage[set][replicate] += line.coverage;
				}
			}
			if (last_nucleotide_is_open == true) {//the last nucleotide of the file
				if (data[index-1].coverage < MINIMUM_READ_COVERAGE) {
					failed_coverage++;
					index--;
				} else {
					passed_coverage++;
					global_methylation[set][replicate] += data[index-1].methylated_C;
					global_coverage[set][replicate] += data[index-1].coverage;
				}
			}
			fclose(BISULFITE_DATA_file);
			free(bs_line); bs_line = NULL;
			free(previous_chromosome); previous_chromosome = NULL;
			if (index == 0) {
				printf("\nWARNING\n\nno data points were added to BISULFITE_DATA[%u][%u]\n\n",set, replicate);
			} else if (DEBUG == true) {
				printf("last written data point was BISULFITE_DATA[%u][%u][%u]\n",set, replicate, index-1);
			}
			const unsigned int TOTAL_READ = failed_coverage + passed_coverage;
			const double PERCENT_PASSING = TOTAL_READ > 0 ? (double)100*passed_coverage/TOTAL_READ : 0.0;
			if (two_strands == true) {
				fprintf(run_info,"%s: %u nucleotide pairs failed coverage, and %u nucleotide pairs passed coverage = %.1f%% passing minimum coverage.\n",files,failed_coverage,passed_coverage,PERCENT_PASSING);
			} else {
				fprintf(run_info,"%s: %u nucleotides failed coverage, and %u nucleotides passed coverage = %.1f%% passing minimum coverage.\n",files,failed_coverage,passed_coverage,PERCENT_PASSING);
			}
			fprintf(run_info,"Global methylation percent = %.1f; %u CpN, sum(5mC) = %zu, sum(coverage) = %zu\n",global_coverage[set][replicate] > 0 ? 100.0*(double)global_methylation[set][replicate] / (double)global_coverage[set][replicate] : 0.0, index, global_methylation[set][replicate], global_coverage[set][replicate]);
//trim the arrays to the size that was used; there is always an "END" element after the last chromosome
			changes[chromosome_count].index = index;
			changes[chromosome_count].chromosome = strdup("END");
			changes = realloc(changes, ((size_t)chromosome_count + 1) * sizeof(CHROMOSOME_CHANGES_STRUCT));
			data = realloc(data, (index > 0 ? index : 1) * sizeof(BISULFITE_DATA_STRUCT));
			if ((changes == NULL) || (data == NULL)) {
				fprintf(stderr, "@line %u	realloc failed after file %s was read.\n", __LINE__, files);
				perror("");
				exit(EXIT_FAILURE);
			}
			BISULFITE_DATA[set][replicate] = data;
			bs_data_chromosome_changes[set][replicate] = changes;
			number_of_data_chromosomes[set][replicate] = chromosome_count;
			replicate++;
//			printf("there were %u chromosomes detected in set %u replicate %u\n", chromosome_count, set + 1, replicate);
//ABSOLUTELY **NOTHING** MORE IN THIS LOOP AFTER THE FOLLOWING COMMAND
			files = strtok_r(NULL,",",&files_pointer);
		}
	}
//print global methylation comparisons
	for (unsigned short int set1 = 0; set1 < SETS; set1++) {//global_methylation may have > 2 sets, so I must do a 2x2 comparison
		for (unsigned short int set2 = 1; set2 < SETS; set2++) {
			if (set2 <= set1) {
				continue;
	 		}
	 		DMR_struct **restrict global = malloc((size_t)(number_of_replicates_per_set[set2] + number_of_replicates_per_set[set1]) * sizeof(DMR_struct) + sizeof(char));
	 		global[0] = malloc(number_of_replicates_per_set[set1] * sizeof(DMR_struct));
	 		global[1] = malloc(number_of_replicates_per_set[set2] * sizeof(DMR_struct));
	 		const unsigned short int COMPARED_SETS[] = {set1,set2};//for references to global variables
			for (unsigned short int set = 0; set < 2; set++) {
				for (unsigned int replicate = 0; replicate < number_of_replicates_per_set[COMPARED_SETS[set]]; replicate++) {
					global[set][replicate].methylated_C = global_methylation[COMPARED_SETS[set]][replicate];
					global[set][replicate].coverage     = global_coverage[COMPARED_SETS[set]][replicate];
				}
			}
			const unsigned short int NUMBER_OF_REPLICATES_COMPARED[] = {number_of_replicates_per_set[COMPARED_SETS[0]], number_of_replicates_per_set[COMPARED_SETS[1]]};
	 		fprintf(run_info,"\nSet %u vs set %u global methylation:\tp-value = %g\t",set1+1,set2+1,Pvalue(global, NUMBER_OF_REPLICATES_COMPARED, options.paired));
	 		fprintf(run_info,"Percent Change = %g%%\n\n", 100.0*mean_percent_difference(global, NUMBER_OF_REPLICATES_COMPARED));
			free(global[0]); global[0] = NULL;
			free(global[1]); global[1] = NULL;
			free(global); global = NULL;
	 	}
		free(global_methylation[set1]); global_methylation[set1] = NULL;
		free(global_coverage[set1]);	global_coverage[set1] = NULL;
	}
	free(global_coverage); global_coverage = NULL;
	free(global_methylation); global_methylation = NULL;
	fprintf(run_info,"There will be %u runs split among %u workers.\n", RUNS, options.CPU);
	fclose(run_info);
	printf("Details of the data have been written to %s.\n",run_info_filename);
	free(run_info_filename); run_info_filename = NULL;
	printf("There will be %u runs split among %u workers.\n", RUNS, options.CPU);
//---------------------------------------------------------------------------------------------
//Analysis (heart of program)
//---------------------------------------------------------------------------------------------
#ifdef _OPENMP
	omp_set_num_threads(options.CPU);
#endif
	puts("\n------\nBeginning analysis.\n");
	for (unsigned short int set1 = 0; set1 < SETS; set1++) {
		for (unsigned short int set2 = 1; set2 < SETS; set2++) {
			if (set2 <= set1) {
				continue;
	 		}
			cutoff_struct *restrict cutoff = make_cutoff_struct(coverage, CpN, d, p, percent, skip, default_option_range, options.ALLOWED_NON_DEFAULTS, options.CPU);
	 		char *restrict label_set1 = malloc(16), *restrict label_set2 = malloc(16);
	 		if (label_set1 == NULL) {
	 			fprintf(stderr, "Failed to malloc at %s line %u\n", __FILE__, __LINE__);
	 			perror("");
	 			exit(EXIT_FAILURE);
	 		}
	 		if (label_set2 == NULL) {
	 			fprintf(stderr, "Failed to malloc at %s line %u\n", __FILE__, __LINE__);
	 			perror("");
	 			exit(EXIT_FAILURE);
	 		}
	 		if (number_of_labels >= (unsigned)(set1+1)) {//does this set have a label?
	 			label_set1 = realloc(label_set1, sizeof(char*) + strlen(labels[set1]));
	 			if (label_set1 == NULL) {
	 				printf("realloc failed at %s line %u\n", __FILE__, __LINE__);
	 				perror("");
	 				exit(EXIT_FAILURE);
	 			}
	 			safe_string_copy(label_set1, labels[set1]);
	 		} else {
	 			sprintf(label_set1, "set%u", set1+1);
	 		}
	 		if (number_of_labels >= (unsigned)(set2+1)) {//does this set have a label?
	 			label_set2 = realloc(label_set2, sizeof(char*) + strlen(labels[set2]));
	 			if (label_set2 == NULL) {
	 				printf("realloc failed at %s line %u\n", __FILE__, __LINE__);
	 				perror("");
	 				exit(EXIT_FAILURE);
	 			}
	 			safe_string_copy(label_set2, labels[set2]);
	 		} else {
	 			sprintf(label_set2, "set%u", set2+1);
	 		}
	 		if (DEBUG == true) {
	 			printf("Now comparing set %u (%s) against set %u (%s) @ %s line %u\n", set1, label_set1, set2, label_set2, __FILE__, __LINE__);
	 		}
			if ((options.paired == true) && (number_of_replicates_per_set[set1] != number_of_replicates_per_set[set2])) {
				puts("In a paired analysis, the number of replicates must be equal.");
				fprintf(stderr, "%u replicates in set %u were attempted to compare with %u replicates in set %u\n", number_of_replicates_per_set[set1], set1+1, number_of_replicates_per_set[set2], set2+1);
				exit(EXIT_FAILURE);
			}
	 		const unsigned short int COMPARED_SETS[] = {set1,set2};//for references to global variables
			const unsigned short int NUMBER_OF_REPLICATES_COMPARED[] = {number_of_replicates_per_set[set1], number_of_replicates_per_set[set2]};
			const unsigned short int TOTAL_REPLICATES = (unsigned short)(NUMBER_OF_REPLICATES_COMPARED[0] + NUMBER_OF_REPLICATES_COMPARED[1]);
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (unsigned int run = 0; run < RUNS; run++) {//for each cutoff group
//-------------------------------------------------
//open answers file and write header
//-------------------------------------------------
//alter label
/*				if (DEBUG == true) {
					printf("now at run %u/%u at %s line %u.\n", run, RUNS, __FILE__, __LINE__);
				}*/
				char options_string[128];
//each piece is appended at the end of what is already written: sprintf(options_string, "%s...", options_string) is undefined behavior, and loses the start of the string when compiled with _FORTIFY_SOURCE
				int written = snprintf(options_string, sizeof(options_string), "c%u_CpN%u_d%u", cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d);
				if (options.allowed_nucleotide_gap_between_CpN != 20000) {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_G%u", options.allowed_nucleotide_gap_between_CpN);
				}
//p = 0.05 is "p0.05" & 10% is "P10" as always, but p = 0.001 is "p0.001" instead of "p0.00", so that runs with different cutoffs get different filenames
				if (doubles_equal(cutoff[run].p * 100.0, round(cutoff[run].p * 100.0)) == true) {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_p%.2f", cutoff[run].p);
				} else {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_p%g", cutoff[run].p);
				}
				if (doubles_equal(cutoff[run].percent * 100.0, round(cutoff[run].percent * 100.0)) == true) {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_P%.0f", cutoff[run].percent * 100.0);
				} else {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_P%g", cutoff[run].percent * 100.0);
				}
				if (cutoff[run].skip >= (UINT_MAX-1)) {//"-S inf"
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_Sinf");
				} else if (cutoff[run].skip > 0) {
					written += snprintf(options_string + written, sizeof(options_string) - (size_t)written, "_S%u", cutoff[run].skip);
				}
//file_suffix is "<-l label>_<cutoffs>", and is in the name of every file that this run writes
				const size_t FILE_SUFFIX_SIZE = strlen(options_string) + (total_label != NULL ? strlen(total_label) + 1 : 0) + 1;
				char *restrict file_suffix = checked_malloc(FILE_SUFFIX_SIZE, __LINE__);
				if (total_label != NULL) {
					snprintf(file_suffix, FILE_SUFFIX_SIZE, "%s_%s", total_label, options_string);
				} else {
					snprintf(file_suffix, FILE_SUFFIX_SIZE, "%s", options_string);
				}
				const size_t FILENAME_BASE_SIZE = strlen(label_set1) + strlen("_vs_") + strlen(label_set2) + strlen("_") + strlen(file_suffix) + 1;
				const size_t ANSWERS_FILENAME_SIZE = FILENAME_BASE_SIZE + strlen(".tsv");
				char *restrict answers_filename = checked_malloc(ANSWERS_FILENAME_SIZE, __LINE__);
				snprintf(answers_filename, ANSWERS_FILENAME_SIZE, "%s_vs_%s_%s.tsv", label_set1, label_set2, file_suffix);
				char *restrict heatmap_tsv_filename = NULL;
				if (options.make_heatmap == true) {
					const size_t HEATMAP_FILENAME_SIZE = FILENAME_BASE_SIZE + strlen(".heatmap.tsv");
					heatmap_tsv_filename = checked_malloc(HEATMAP_FILENAME_SIZE, __LINE__);
					snprintf(heatmap_tsv_filename, HEATMAP_FILENAME_SIZE, "%s_vs_%s_%s.heatmap.tsv", label_set1, label_set2, file_suffix);
				}
/*				if ((options.overwrite_output == false) && (file_exists(answers_filename) == true)) {
					printf("Answers file %s already exists, so Defiant is skipping this iteration.\n", answers_filename);
					free(answers_filename); answers_filename = NULL;
					continue;
				}*/
//----make answers file-----------------
				FILE *restrict rscript = NULL;
				char rscript_filename[48];
				if (options.make_figures == true) {
					sprintf(rscript_filename, "DMR_figures/defiant_%u.R", run);
					rscript = fopen(rscript_filename,"w");
					if (rscript == NULL) {
						perror("failed to write DMR_figures/defiant.R.\n");
						exit(EXIT_FAILURE);
					}
					fprintf(rscript, "legend_space <- -0.26\n");
					fprintf(rscript, "right <- 10.5\n");
					fprintf(rscript, "bottom <- 0\n");
					fprintf(rscript, "left <- 5\n");
					fprintf(rscript, "top <- 0\ncex_main = 1\n");
				}
				FILE *restrict answers;
				answers = fopen(answers_filename, "w");
				if (answers == NULL) {
					printf("\nfailed to write %s\n", answers_filename);
					perror("Failed to open file.\n");
					free(answers_filename); answers_filename = NULL;
					exit(EXIT_FAILURE);
				}
				FILE *restrict fh_heatmap = NULL;
				if (options.make_heatmap == true) {
					fh_heatmap = fopen(heatmap_tsv_filename, "w");
					if (fh_heatmap == NULL) {
						printf("Failed to write %s at %s line %u\n", heatmap_tsv_filename, __FILE__, __LINE__);
						perror("");
						free(answers_filename); answers_filename = NULL;
						free(heatmap_tsv_filename); heatmap_tsv_filename = NULL;
						exit(EXIT_FAILURE);
					}
					//now write header for heatmap input TSV
					fprintf(fh_heatmap, "Gene.Name");
					if (label_set1 == NULL) {
						for (unsigned short r = 0; r < NUMBER_OF_REPLICATES_COMPARED[0]; r++) {
							fprintf(fh_heatmap, "\tSet1.R%u", r+1);
						}
					} else {
						for (unsigned short r = 0; r < NUMBER_OF_REPLICATES_COMPARED[0]; r++) {
							fprintf(fh_heatmap, "\t%s.R%u", label_set1, r+1);
						}
					}
					if (label_set2 == NULL) {
						for (unsigned short r = 0; r < NUMBER_OF_REPLICATES_COMPARED[1]; r++) {
							fprintf(fh_heatmap, "\tSet2.R%u", r+1);
						}
					} else {
						for (unsigned short r = 0; r < NUMBER_OF_REPLICATES_COMPARED[1]; r++) {
							fprintf(fh_heatmap, "\t%s.R%u", label_set2, r+1);
						}
					}
					fprintf(fh_heatmap, "\n");
				}
				fprintf(answers,"Chromosome\tStart\tEnd\t#mCpN\t#Diff.CpN");//despite all options, this will be in the header
				fprintf(answers,"\tMean_Difference\t%s\t%s",label_set1, label_set2);//answers file header
				if (cutoff[run].skip > 0) {
					fprintf(answers,"\tTotal_Skips");
				}
				if (gene_annotation_defined == true) {//part of the answer file header
					fprintf(answers,"\tInside_Genes\tBetween_Genes\tGene_Promoter_Cutoff_%u_Nucleotides", promoter_cutoff);
				}
				if (options.printing_CpN_list == true) {//part of the answer file header
					fprintf(answers,"\tDiff.CpN.list");
				}
/*				if (options.isolation == true) {//part of the answer file header
					fprintf(answers,"\tNearest.CpN.Distances");
				}*/
				fprintf(answers,"\n");//part of the answer file header
#ifdef _OPENMP
				printf("\n\nWorker %u/%u is analyzing run %u/%u:	sets %u vs. %u	with cutoffs	c = %u	CpN = %u	d = %u	p = %.2f	P = %lf	S = %u\n", omp_get_thread_num()+1, options.CPU, run+1, RUNS, set1+1, set2+1, cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p, cutoff[run].percent*100, cutoff[run].skip);
#endif
				unsigned int min_gene_search_index = 0, max_gene_search_index = 0, first_nucleotide = 0;
				size_t methylation_up_gene_list_size = 1;
				size_t methylation_down_gene_list_size = 1;
				size_t gene_list_size = 1;
				char **restrict methylation_up_gene_list = NULL;
				char **restrict methylation_down_gene_list = NULL;
				char **restrict gene_list = NULL;
				if (gene_annotation_defined == true) {
					methylation_up_gene_list = malloc(1);
					methylation_down_gene_list = malloc(1);
					gene_list = malloc(1);
				}
				unsigned int number_of_methylation_up_list_genes = 0, number_of_methylation_down_list_genes = 0;
				unsigned int number_of_list_genes = 0;
				unsigned int MIN_COVERAGE = cutoff[run].skip > 0 ? 0 : cutoff[run].coverage;
/*				if (options.isolation == true) {
					MIN_COVERAGE = 0;
				} else */
				if (options.print_all_nucleotides == true) {
					MIN_COVERAGE = 0;
				}
				FILE *restrict fh_print_all_nucleotides = NULL, *restrict bedfile_fh = NULL;
				if (options.print_all_nucleotides == true) {//printing to each individual nucleotide information
//with more than 1 run the cutoffs are in the name, so that the runs don't write over each other's file
					const size_t FILENAME_SIZE = strlen("nucleotide_") + strlen(label_set1) + strlen(label_set2) + strlen(file_suffix) + strlen(".tsv") + 3;
					char *restrict each_nucleotide_filename = checked_malloc(FILENAME_SIZE, __LINE__);
					if (RUNS > 1) {
						snprintf(each_nucleotide_filename, FILENAME_SIZE, "nucleotide_%s_%s_%s.tsv", label_set1, label_set2, file_suffix);
					} else if (total_label != NULL) {
						snprintf(each_nucleotide_filename, FILENAME_SIZE, "nucleotide_%s_%s_%s.tsv",label_set1, label_set2, total_label);
					} else {
						snprintf(each_nucleotide_filename, FILENAME_SIZE, "nucleotide_%s_%s.tsv",label_set1, label_set2);
					}
					fh_print_all_nucleotides = fopen(each_nucleotide_filename, "w");
					if (fh_print_all_nucleotides == NULL) {
						fprintf(stderr, "\nfailed to write %s\n\n", each_nucleotide_filename);
						perror("");
						free(answers_filename); answers_filename = NULL;
						exit(EXIT_FAILURE);
					}
					if (options.fdr_on_all_cpn == false) {
						fprintf(fh_print_all_nucleotides,"Chr\tNucleotide\t%%Difference\tp\t%s\t%s\tmC_count\tSimilar.Nucleotides\tSkipped.Nucleotides\tRunning.%%Diff\tRunning.p\n",label_set1, label_set2);
					} else {
						fprintf(fh_print_all_nucleotides,"Chr\tNucleotide\t%%Difference\tp\tq (%s)\t%s\t%s\tmC_count\tSimilar.Nucleotides\tSkipped.Nucleotides\tRunning.%%Diff\tRunning.p\n", FDR_type, label_set1, label_set2);
					}
					printf("\nNow writing %s\n\n",each_nucleotide_filename);
					free(each_nucleotide_filename); each_nucleotide_filename = NULL;
				}
				char *restrict bedfile_filename = NULL;
				if (options.make_bed_files == true) {
					const size_t BEDFILE_FILENAME_SIZE = strlen(label_set1) + strlen(label_set2) + strlen(file_suffix) + strlen(".bed") + 3;
					bedfile_filename = checked_malloc(BEDFILE_FILENAME_SIZE, __LINE__);
					if (RUNS > 1) {
						snprintf(bedfile_filename, BEDFILE_FILENAME_SIZE, "%s_%s_%s.bed", label_set1, label_set2, file_suffix);
					} else if (total_label != NULL) {
						snprintf(bedfile_filename, BEDFILE_FILENAME_SIZE, "%s_%s_%s.bed", label_set1, label_set2, total_label);
					} else {
						snprintf(bedfile_filename, BEDFILE_FILENAME_SIZE, "%s_%s.bed",label_set1, label_set2);
					}
					bedfile_fh = fopen(bedfile_filename, "w");
					if (bedfile_fh == NULL) {
						perror("\nFailed to open bedfile\n");
						fprintf(stderr, "Failed to open %s\n", bedfile_filename);
						fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
						free(answers_filename); answers_filename = NULL;
						free(bedfile_filename); bedfile_filename = NULL;
						exit(EXIT_FAILURE);
					}
					fprintf(bedfile_fh, "chrom chromStart chromEnd\n");
				}
//-------------------------------------------------
//finished opening answers file and writing header
//-------------------------------------------------
//		may need an array of p-values to determine the q-value ahead of time, if fdr == true
//-------------------------------------------------
//-------------------------------------------------
				double *restrict pvalues = NULL;//only used if fdr option is set to 'true'
				double *restrict qvalues = NULL;//only used if fdr option is set to 'true'
				unsigned int pvalue_index = 0;
				if (options.fdr_on_all_cpn == true) {//get an array of p-values for FDR computation
//pvalues & qvalues have one element for every nucleotide that the analysis below goes through, in the same order.
//only nucleotides with minimum coverage in every replicate can be tested, so only those are counted as tests when adjusting; the others get q = 1
					bool *restrict testable = NULL;
					size_t testable_capacity = 0;
					unsigned int number_of_testable = 0;
					for (unsigned int chromosome = 0; chromosome < number_of_data_chromosomes[set1][0]; chromosome++) {//this is done for each comparison
						printf("Gathering p-values for p-value adjustment/FDR for chromosome %s...\n", bs_data_chromosome_changes[set1][0][chromosome].chromosome);
						DMR_struct **restrict single_nucleotide = checked_malloc(2 * sizeof(DMR_struct*), __LINE__);
						size_t **restrict bisulfite_nucleotide_min = checked_malloc(2 * sizeof(size_t*), __LINE__);
						size_t **restrict bisulfite_nucleotide_max = checked_malloc(2 * sizeof(size_t*), __LINE__);
						for (unsigned short int set = 0; set < 2; set++) {
							single_nucleotide[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(DMR_struct), __LINE__);
							bisulfite_nucleotide_min[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(size_t), __LINE__);
							bisulfite_nucleotide_max[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(size_t), __LINE__);
						}
//now get the minimum index in BISULFITE_DATA where the chromosome matches, as Bismark does not sort by chromosome
						const bool CHROMOSOME_PRESENT = chromosome_in_every_replicate(bs_data_chromosome_changes[set1][0][chromosome].chromosome, bs_data_chromosome_changes, number_of_data_chromosomes, COMPARED_SETS, NUMBER_OF_REPLICATES_COMPARED, bisulfite_nucleotide_min, bisulfite_nucleotide_max);
						if (CHROMOSOME_PRESENT == false) {//I put this here because of Ayma's data, which had 0 coverage in an entire chromosome
							printf("Data is too sparse to consider chromosome %s, continuing to next chromosome.\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome);
						}
						while (CHROMOSOME_PRESENT == true) {//keeps looping until a replicate runs out of nucleotides
							bool breaker = false;
							zero_DMR_struct(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);//Initialize to 0
							unsigned int nucleotide = UINT_MAX;//as big as I can get it
//go through all sets to get the minimum nucleotide
							for (unsigned short int set = 0; set < 2; set++) {//get the new minimum nucleotide 
								for (unsigned short int r = 0; r < NUMBER_OF_REPLICATES_COMPARED[set]; r++) {
									if (bisulfite_nucleotide_min[set][r] >= bisulfite_nucleotide_max[set][r]) {
										breaker = true;
										break;//prevent segmentation faults
									}
									if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide < nucleotide) {
										nucleotide = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide;
									}
								}
							}
							if (breaker == true) {//end this chromosome
								break;
							}
							unsigned short int replicates_meeting_minimum_coverage = 0;
				 			for (unsigned short int set = 0; set < 2; set++) {//how many organisms have this nucleotide?
								for (unsigned short int r = 0; r < NUMBER_OF_REPLICATES_COMPARED[set]; r++) {
									if (nucleotide < BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide) {
										continue;//I haven't yet reached this nucleotide, continue to next individual
									}
									if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide == nucleotide) {
										if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage >= MIN_COVERAGE) {
											if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage == 0) {
												printf("got a bad point: set = %u, COMPARED_SETS[%u] = %u\n", set, set, COMPARED_SETS[set]);
												exit(EXIT_FAILURE);
											}
											single_nucleotide[set][r].methylated_C = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].methylated_C;
											single_nucleotide[set][r].coverage     = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage;
											if (single_nucleotide[set][r].coverage >= cutoff[run].coverage) {
												replicates_meeting_minimum_coverage++;
											}
										}
									}
									bisulfite_nucleotide_min[set][r]++;
								}
							}
							if (pvalue_index == UINT_MAX) {
								fprintf(stderr, "There are more than %u nucleotides, which is too many for '-fdr'.\n", UINT_MAX);
								exit(EXIT_FAILURE);
							}
							pvalues = add_double(pvalues, Pvalue(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED, options.paired), pvalue_index);
							if (pvalue_index >= testable_capacity) {
								testable_capacity = testable_capacity == 0 ? 65536 : 2 * testable_capacity;
								testable = realloc(testable, testable_capacity * sizeof(bool));
								if (testable == NULL) {
									fprintf(stderr, "realloc failed @ %s line %u\n", __FILE__, __LINE__);
									exit(EXIT_FAILURE);
								}
							}
							testable[pvalue_index] = (replicates_meeting_minimum_coverage == TOTAL_REPLICATES);
							if (testable[pvalue_index] == true) {
								number_of_testable++;
							}
							pvalue_index++;
						}
						for (unsigned short int set = 0; set < 2; set++) {//there are only 2 sets compared at a time
							free(bisulfite_nucleotide_min[set]); bisulfite_nucleotide_min[set] = NULL;
							free(bisulfite_nucleotide_max[set]); bisulfite_nucleotide_max[set] = NULL;
							free(single_nucleotide[set]); single_nucleotide[set] = NULL;
						}						
						free(single_nucleotide); single_nucleotide = NULL;
						free(bisulfite_nucleotide_min); bisulfite_nucleotide_min = NULL;
						free(bisulfite_nucleotide_max); bisulfite_nucleotide_max = NULL;
					}//end of each chromosome
					qvalues = checked_malloc((size_t)pvalue_index * sizeof(double), __LINE__);
					for (unsigned int i = 0; i < pvalue_index; i++) {
						qvalues[i] = 1.0;
					}
					if (number_of_testable > 0) {
						double *restrict testable_pvalues = checked_malloc((size_t)number_of_testable * sizeof(double), __LINE__);
						unsigned int t = 0;
						for (unsigned int i = 0; i < pvalue_index; i++) {
							if (testable[i] == true) {
								testable_pvalues[t] = pvalues[i];
								t++;
							}
						}
						double *restrict testable_qvalues = p_adjust(testable_pvalues, number_of_testable, FDR_type);
						t = 0;
						for (unsigned int i = 0; i < pvalue_index; i++) {
							if (testable[i] == true) {
								qvalues[i] = testable_qvalues[t];
								t++;
							}
						}
						free(testable_pvalues); testable_pvalues = NULL;
						free(testable_qvalues); testable_qvalues = NULL;
					}
					printf("%u of %u nucleotides have minimum coverage in every replicate, and were adjusted with %s's method.\n", number_of_testable, pvalue_index, FDR_type);
					free(testable); testable = NULL;
				}//end of FDR 
//-------------------------------------------------
// Finished getting q-values
//-------------------------------------------------
				unsigned int nucleotides_good_in_every_organism = 0, total_nucleotides = 0;
				if (number_of_data_chromosomes[set1][0] == 0) {
					fprintf(stderr, "there are 0 chromosomes.\n");
					fprintf(stderr, "\nFAILED at %s line %u\n", __FILE__, __LINE__);
					free(answers_filename); answers_filename = NULL;
					free(bedfile_filename); bedfile_filename = NULL;
					exit(EXIT_FAILURE);
				}
				DMR_struct ***restrict all_DMRs = NULL;
				unsigned int *restrict DMR_chromosome_indices = NULL;//only used if DMR FDR is 'true'
				unsigned int *restrict DMR_starts = NULL;
				unsigned int *restrict DMR_ends   = NULL;
				unsigned int *restrict nCpN       = NULL;
				unsigned int *restrict diffCpN    = NULL;
				unsigned int *restrict CpN_skips_in_this_DMR_all = NULL;
				if (options.print_DMR_pvalue == true) {
					if (DEBUG == true) {printf("at %s line %u\n", __FILE__, __LINE__);}
					all_DMRs = malloc_nucleotides(all_DMRs, NUMBER_OF_REPLICATES_COMPARED);
					if (DEBUG == true) {printf("at %s line %u\n", __FILE__, __LINE__);}
				}
				unsigned int number_of_DMRs = 0;
				for (unsigned int chromosome = 0; chromosome < number_of_data_chromosomes[set1][0]; chromosome++) {//this is done for each comparison
					if (options.CPU == 1) {
						printf("Now analyzing chromosome %s, at iteration loop %u/%u\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome, chromosome+1, number_of_data_chromosomes[set1][0]);
					}
					size_t **restrict bisulfite_nucleotide_min = checked_malloc(2 * sizeof(size_t*), __LINE__);
					size_t **restrict bisulfite_nucleotide_max = checked_malloc(2 * sizeof(size_t*), __LINE__);
					for (unsigned short int set = 0; set < 2; set++) {
						bisulfite_nucleotide_min[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(size_t), __LINE__);
						bisulfite_nucleotide_max[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(size_t), __LINE__);
					}
//now get the minimum index in BISULFITE_DATA where the chromosome matches, as Bismark does not sort by chromosome
					if (chromosome_in_every_replicate(bs_data_chromosome_changes[set1][0][chromosome].chromosome, bs_data_chromosome_changes, number_of_data_chromosomes, COMPARED_SETS, NUMBER_OF_REPLICATES_COMPARED, bisulfite_nucleotide_min, bisulfite_nucleotide_max) == false) {//I put this here because of Ayma's data, which had 0 coverage in an entire chromosome
						printf("Data is too sparse to consider chromosome %s, continuing to next chromosome.\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome);
						for (unsigned short int set = 0; set < 2; set++) {
							free(bisulfite_nucleotide_min[set]); bisulfite_nucleotide_min[set] = NULL;
							free(bisulfite_nucleotide_max[set]); bisulfite_nucleotide_max[set] = NULL;
						}
						free(bisulfite_nucleotide_min); bisulfite_nucleotide_min = NULL;
						free(bisulfite_nucleotide_max); bisulfite_nucleotide_max = NULL;
						continue;
					}
					DMR_struct **restrict DMR_sum           = checked_malloc(2 * sizeof(DMR_struct*), __LINE__);
					DMR_struct **restrict single_nucleotide = checked_malloc(2 * sizeof(DMR_struct*), __LINE__);
					DMR_struct **restrict DMR_similar_sum   = checked_malloc(2 * sizeof(DMR_struct*), __LINE__);
					DMR_struct ***restrict all_nucleotides  = NULL;
					if (options.make_figures == true) {
						all_nucleotides = malloc_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED);
					}
					for (unsigned short int set = 0; set < 2; set++) {
						DMR_sum[set]           = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(DMR_struct), __LINE__);
						single_nucleotide[set] = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(DMR_struct), __LINE__);
						DMR_similar_sum[set]   = checked_malloc(NUMBER_OF_REPLICATES_COMPARED[set] * sizeof(DMR_struct), __LINE__);
					}
//ensure chromosome loop in GENES and BISULFITE_DATA matches
					bool chromosome_annotation_found = false;//genes are only looked up on chromosomes that the annotation has
					if (gene_annotation_defined == true) {//find which indices of GENES belong to this chromosome
		 				for (unsigned int chr_loop = 0; chr_loop < number_of_chromosomes_in_GENES; chr_loop++) {
							if (strcmp(gene_chromosome_changes[chr_loop].chromosome,bs_data_chromosome_changes[set1][0][chromosome].chromosome) == 0) {
								chromosome_annotation_found = true;
								GENES_chromosome_change_index = chr_loop;
								min_gene_search_index = gene_chromosome_changes[chr_loop].index;
//this chromosome's genes end right before the next annotated chromosome's genes start, and the last annotated chromosome's genes end with the last gene.
//chr_loop is the chromosome's index in the annotation, which isn't its index in the data
								if ((chr_loop + 1) < number_of_chromosomes_in_GENES) {
									max_gene_search_index = gene_chromosome_changes[chr_loop + 1].index - 1;
								} else {
									max_gene_search_index = number_of_annotated_genes - 1;
								}
								first_nucleotide = GENES[min_gene_search_index].gene_start;
								if (DEBUG == true) {
									printf("@ line %u, GENES_chromosome_change_index = %u -> searching %u-%u indices\n", __LINE__, GENES_chromosome_change_index, min_gene_search_index, max_gene_search_index);
									printf("@ line %u -> 1st nucleotide = %u; max nucleotide in chromosome = %u\n",__LINE__,first_nucleotide, GENES[max_gene_search_index].gene_end);
								}
								break;
							}
						}
						if (chromosome_annotation_found == false) {
							printf("Warning: chromosome %s has no annotation data, its DMRs will have \"n/a\" for genes.\n", bs_data_chromosome_changes[set1][0][chromosome].chromosome);
						}
					}
					if (DEBUG == true) {printf("@line %u\n", __LINE__);}
//following arrays are only for making figures
					unsigned int *restrict CpN_list = NULL,	*restrict all_CpN = NULL;
					double *restrict all_p = NULL,				*restrict all_percent = NULL;
					double *restrict all_sum_p = NULL,			*restrict all_sum_percent = NULL;
					unsigned int mC_start = 0, mC_end = 0;//, last_CpN_before_DMR = 0;
//tail_count: CpN (similar or skipped) since the last differentially methylated CpN; these come after the DMR's end, so they aren't part of it
//tail_skips: how many of those were skipped
					unsigned int tail_count = 0, tail_skips = 0;
					unsigned int mC_count = 0,	similar_nucleotide_count = 0, good_nucleotide_count = 0, consecutive_skipped_CpN = 0, CpN_skips_in_this_DMR = 0;
							tail_count = 0;
							tail_skips = 0;
					zero_DMR_struct(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);//Initialize to 0
					zero_DMR_struct(DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);//Initialize to 0
					zero_DMR_struct(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);//Initialize to 0
//					if (DEBUG == true) {printf("@line %u\n", __LINE__);}
					while (true) {//always true, to keep it looping
						bool breaker = false;
						zero_DMR_struct(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);//Initialize to 0
						unsigned int nucleotide = UINT_MAX;//as big as I can get it
//go through all sets to get the minimum nucleotide
//					if (DEBUG == true) {printf("@line %u\n", __LINE__);}
						for (unsigned short int set = 0; set < 2; set++) {//get the new minimum nucleotide 
							for (unsigned short int r = 0; r < NUMBER_OF_REPLICATES_COMPARED[set]; r++) {
//							 	printf("line%u-> index for BISULFITE_DATA bisulfite_nucleotide_min[%u][%u] = %zu	bisulfite_nucleotide_max = %zu\n",__LINE__, COMPARED_SETS[set], r, bisulfite_nucleotide_min[set][r], bisulfite_nucleotide_max[set][r]);
								if (bisulfite_nucleotide_min[set][r] >= bisulfite_nucleotide_max[set][r]) {
									breaker = true;
									break;//prevent segmentation faults
								}
//							 	printf("line%u-> nucleotide from above: %u\n",__LINE__,BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide);
//							 	printf("@line %u bisulfite_nucleotide_max[%u][%u] index = \n",__LINE__,set,r);
//								printf("@line %u BISULFITE_DATA[%u][%u][%zu].nucleotide = %u\n",__LINE__,COMPARED_SETS[set],r,bisulfite_nucleotide_min[set][r],BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide);
								if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide < nucleotide) {
									nucleotide = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide;
								}
							}
						}
						if (breaker == true) {//end this chromosome
							if (((mC_count - tail_count) >= cutoff[run].CpN) && (good_nucleotide_count >= cutoff[run].d)) {
								const unsigned int NO_COAT_TAILS = mC_count - tail_count;
								const unsigned int SKIPS_IN_DMR = CpN_skips_in_this_DMR - tail_skips;
								if ((NO_COAT_TAILS >= cutoff[run].CpN) && ((mC_end-mC_start) >= options.MINIMUM_CPN_RANGE) && (good_nucleotide_count >= cutoff[run].d)) {
									subtract_2nd_dmr_from_1st_dmr(DMR_sum, DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
									const double SUM_PERCENT_DIFFERENCE = mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED);
									if (fabs(SUM_PERCENT_DIFFERENCE) >= cutoff[run].percent) {
										if (options.make_bed_files == true) {
											fprintf(bedfile_fh,"chr%s	%u	%u\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome,mC_start ,mC_end);
										}
										fprintf(answers,"%s	%u	%u\t%u	%u",bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start , mC_end, NO_COAT_TAILS, good_nucleotide_count);
										fprintf(answers,"\t%.1f", 100.0*SUM_PERCENT_DIFFERENCE);
										if (options.print_DMR_pvalue == true) {
											if (DEBUG == true) {printf("@line %u\n", __LINE__);}
											all_DMRs = add_nucleotides(all_DMRs, NUMBER_OF_REPLICATES_COMPARED, number_of_DMRs, DMR_sum);
											DMR_chromosome_indices = add_uint(DMR_chromosome_indices, chromosome, number_of_DMRs);
											DMR_starts             = add_uint(DMR_starts, mC_start, number_of_DMRs);
											DMR_ends               = add_uint(DMR_ends, mC_end, number_of_DMRs);
											nCpN                   = add_uint(nCpN, NO_COAT_TAILS, number_of_DMRs);
											diffCpN                = add_uint(diffCpN, good_nucleotide_count, number_of_DMRs);
											CpN_skips_in_this_DMR_all  = add_uint(CpN_skips_in_this_DMR_all, SKIPS_IN_DMR, number_of_DMRs);
											if (DEBUG == true) {printf("@line %u\n", __LINE__);}
										}
										print_methylation_percentages(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, answers);
										if (options.make_heatmap == true) {
											print_methylation_heatmap(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, fh_heatmap, number_of_DMRs);
										}
										if (cutoff[run].skip > 0) {
											fprintf(answers,"\t%u", SKIPS_IN_DMR);
										}
										char gene[GENE_STRING_SIZE] = "n/a", intergene[GENE_STRING_SIZE] = "n/a", list_gene[GENE_STRING_SIZE] = "n/a";
										if (gene_annotation_defined == true) {
	//									printf("@ %s line %u list_gene = %s\n", __FILE__, __LINE__, list_gene);
											if (chromosome_annotation_found == true) {identify_gene(mC_start, mC_end, gene, intergene, list_gene, GENE_STRING_SIZE, &min_gene_search_index, &max_gene_search_index, first_nucleotide);}
											fprintf(answers,"\t%s\t%s\t%s", gene, intergene, list_gene);
											gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, gene, false, true);
											gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, intergene, false, true);
											gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, list_gene, false, true);
								   	//void *restrict add_name(char **restrict array, size_t *restrict size, unsigned int *restrict f_number_of_list_genes, const char *restrict NAME, const bool UNIQUE) 
											if (SUM_PERCENT_DIFFERENCE > 0) {
												methylation_up_gene_list = add_name(methylation_up_gene_list, &methylation_up_gene_list_size, &number_of_methylation_up_list_genes, list_gene, true, false);
											} else {
												methylation_down_gene_list = add_name(methylation_down_gene_list, &methylation_down_gene_list_size, &number_of_methylation_down_list_genes, list_gene, true, false);
											}
										}
										if (options.printing_CpN_list == true) {
											print_CpN_list(answers, CpN_list, good_nucleotide_count);
										}
//										if (options.isolation == true) {
//											fprintf(answers, "\t%u;%u", mC_start - last_CpN_before_DMR, get_next_CpN(mC_end, (const size_t **restrict)bisulfite_nucleotide_min, (const size_t **restrict)bisulfite_nucleotide_max, BISULFITE_DATA, NUMBER_OF_REPLICATES_COMPARED, COMPARED_SETS) - mC_end);
	//									}
										fprintf(answers,"\n");
										if (options.make_figures == true) {
											fprintf(rscript, "setEPS()\n");
											fprintf(rscript, "postscript('DMR_figures/%s_%s_c%u_CpN%u_d%u_p0_%02.0f_P%.0f_S%u_chr%s_%u_%u", label_set1, label_set2, cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p*100, 100.0*cutoff[run].percent, cutoff[run].skip, bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start, mC_end);
											if (strcmp(list_gene,"n/a") != 0) {
												fprintf(rscript, "_");
											fprint_for_filename(rscript, list_gene);
											}
											fprintf(rscript, ".eps')\n");
										
											R_figure(rscript, all_CpN, NO_COAT_TAILS, label_set1, label_set2, xaxis_label, list_gene, all_p, all_percent, all_sum_p, all_sum_percent, all_nucleotides, NUMBER_OF_REPLICATES_COMPARED, mC_start, mC_end);
											all_CpN = uint_clear(all_CpN);
											all_p = double_clear(all_p);
											all_percent = double_clear(all_percent);
											all_sum_p = double_clear(all_sum_p);
											all_sum_percent = double_clear(all_sum_percent);
											all_nucleotides = clear_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED);
										}
										number_of_DMRs++;
									}
								}
							}
							break;
						}
//						if (DEBUG == true) {printf("@line %u\n", __LINE__);}
//						printf("@ %s line %u nucleotide = %u\n", __FILE__, __LINE__, nucleotide);
//						printf("%u/%u-chr %s:%u@line %u\n",replicates_meeting_minimum_coverage, TOTAL_REPLICATES, gene_chromosome_changes[chr].chromosome, nucleotide,__LINE__);
						unsigned short int replicates_meeting_minimum_coverage = 0;
			 			for (unsigned short int set = 0; set < 2; set++) {//how many organisms have this nucleotide?
							for (unsigned short int r = 0; r < NUMBER_OF_REPLICATES_COMPARED[set]; r++) {
								if (nucleotide < BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide) {
									continue;//I haven't yet reached this nucleotide, continue to next individual
								}
								if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].nucleotide == nucleotide) {
									if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage >= MIN_COVERAGE) {
/*										if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].methylated_C > BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage) {
											printf("@line %u nucleotide %u has methylated_C %u > coverage %u\n", __LINE__, nucleotide, BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].methylated_C, BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage);
											exit(EXIT_FAILURE);
										}*/
										if (BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage == 0) {
										printf("@ %s line %u\n", __FILE__, __LINE__);
											printf("got a bad point: set = %u, COMPARED_SETS[%u] = %u\n", set, set, COMPARED_SETS[set]);
											exit(EXIT_FAILURE);
										}
										single_nucleotide[set][r].methylated_C = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].methylated_C;
										single_nucleotide[set][r].coverage     = BISULFITE_DATA[COMPARED_SETS[set]][r][bisulfite_nucleotide_min[set][r]].coverage;
										if (single_nucleotide[set][r].coverage >= cutoff[run].coverage) {
											replicates_meeting_minimum_coverage++;
										}
									}
								}
								bisulfite_nucleotide_min[set][r]++;
							}
						}
						total_nucleotides++;
//--------------------------
//PRINTING ALL NUCLEOTIDES
//-------------------------
						if (options.print_all_nucleotides == true) {//printing to each individual nucleotide information, only the columns about this nucleotide, last columns printed later
	//						if (replicates_meeting_minimum_coverage == TOTAL_REPLICATES){//If I can print out p-value and percent change
							if (options.fdr_on_all_cpn == false) {
								fprintf(fh_print_all_nucleotides, "%s\t%u\t%.2f\t%.3g", bs_data_chromosome_changes[set1][0][chromosome].chromosome, nucleotide, 100.0*mean_percent_difference(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED), Pvalue(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							} else {
								fprintf(fh_print_all_nucleotides, "%s\t%u\t%.2f\t%.3g\t%.3g", bs_data_chromosome_changes[set1][0][chromosome].chromosome, nucleotide, 100.0*mean_percent_difference(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED), pvalues[total_nucleotides-1], qvalues[total_nucleotides-1]);
							}
//							} else {//If I don't have enough info to calculate p-value & percent change
	//							fprintf(fh_print_all_nucleotides, "%s\t%u\tn/a\tn/a", bs_data_chromosome_changes[set1][0][chromosome].chromosome, nucleotide);
//							}
//							fprintf(fh_print_all_nucleotides, "\t%u/%u", replicates_meeting_minimum_coverage, TOTAL_REPLICATES);
							for (unsigned short int nset = 0; nset < 2; nset++) {
								if (single_nucleotide[nset][0].coverage > 0) {//nucleotide is present
									fprintf(fh_print_all_nucleotides, "\t%zu/%zu",  single_nucleotide[nset][0].methylated_C, single_nucleotide[nset][0].coverage);
								} else {//nucleotide isn't present
									fprintf(fh_print_all_nucleotides, "\tX");
								}
								for (unsigned short int r = 1; r < NUMBER_OF_REPLICATES_COMPARED[nset]; r++) {
									if (single_nucleotide[nset][r].coverage > 0) {//nucleotide is present
	/*									if (single_nucleotide[0][r].methylated_C > single_nucleotide[0][r].coverage) {
											printf("%zu > %zu\n", single_nucleotide[0][r].methylated_C, single_nucleotide[0][r].coverage);
											exit(EXIT_FAILURE);
										}*/
										fprintf(fh_print_all_nucleotides, ",%zu/%zu",  single_nucleotide[nset][r].methylated_C, single_nucleotide[nset][r].coverage);
									} else {//nucleotide isn't present
										fprintf(fh_print_all_nucleotides, ",X");
									}
								}
							}
//the following conditions will be in further code, but I'm printing here for debuggingoptions.MINIMUM_CPN_RANGE) && (good_nucleotide_count >= cutoff[run].d)) {
						}
//						if (DEBUG == true) {printf("@line %u\n", __LINE__);}
//						printf("%u/%u-chr %s:%u@line %u\n",replicates_meeting_minimum_coverage, TOTAL_REPLICATES, bs_data_chromosome_changes[set1][0][chromosome].chromosome, nucleotide,__LINE__);
//printf("@line %u\n", __LINE__);	
						if ((cutoff[run].skip == 0) && (replicates_meeting_minimum_coverage == 0)) {
//not a single replicate has minimum coverage here. Usually such a nucleotide is thrown away when the files are read, and never gets here,
//but -E, a range of coverages, or a range of skips keep it in memory. It is ignored here too, so that those options don't change which DMRs are found
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t0\t%u\t%u\t%.2f\t%.3g\n", similar_nucleotide_count, consecutive_skipped_CpN, 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
							continue;
						}
						if (((mC_count > 0) && ((nucleotide - mC_end) > options.allowed_nucleotide_gap_between_CpN)) || (similar_nucleotide_count > options.MAX_SIMILAR_CPN) || ((consecutive_skipped_CpN >= cutoff[run].skip) && (replicates_meeting_minimum_coverage < TOTAL_REPLICATES))) {//if the gap is too big, or if this nucleotide doesn't have minimum coverage in every sample, don't use it
//	printf("@line %u\n", __LINE__);
							const unsigned int NO_COAT_TAILS = mC_count - tail_count;
							const unsigned int SKIPS_IN_DMR = CpN_skips_in_this_DMR - tail_skips;
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t0\t%u\t%u\t%.2f\t%.3g\n", similar_nucleotide_count, consecutive_skipped_CpN, 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
//							if (DEBUG == true) {printf("@line %u\n", __LINE__);}
							if ((NO_COAT_TAILS >= cutoff[run].CpN) && ((mC_end-mC_start) >= options.MINIMUM_CPN_RANGE) && (good_nucleotide_count >= cutoff[run].d)) {
								subtract_2nd_dmr_from_1st_dmr(DMR_sum, DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
								const double SUM_PERCENT_DIFFERENCE = mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED);
								if (fabs(SUM_PERCENT_DIFFERENCE) >= cutoff[run].percent) {
									if (options.make_bed_files == true) {
										fprintf(bedfile_fh,"chr%s	%u	%u\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome,mC_start ,mC_end);
									}
									fprintf(answers,"%s	%u	%u\t%u	%u",bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start , mC_end, NO_COAT_TAILS, good_nucleotide_count);
									fprintf(answers,"\t%.1f", 100.0*SUM_PERCENT_DIFFERENCE);
									if (options.print_DMR_pvalue == true) {
										all_DMRs = add_nucleotides(all_DMRs, NUMBER_OF_REPLICATES_COMPARED, number_of_DMRs, DMR_sum);
										DMR_chromosome_indices = add_uint(DMR_chromosome_indices, chromosome, number_of_DMRs);
										DMR_starts             = add_uint(DMR_starts, mC_start, number_of_DMRs);
										DMR_ends               = add_uint(DMR_ends, mC_end, number_of_DMRs);
										nCpN                   = add_uint(nCpN, NO_COAT_TAILS, number_of_DMRs);
										diffCpN               = add_uint(diffCpN, good_nucleotide_count, number_of_DMRs);
										CpN_skips_in_this_DMR_all  = add_uint(CpN_skips_in_this_DMR_all, SKIPS_IN_DMR, number_of_DMRs);
									}
									print_methylation_percentages(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, answers);
									if (options.make_heatmap == true) {
										print_methylation_heatmap(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, fh_heatmap, number_of_DMRs);
									}
//									if (DEBUG == true) {printf("@line %u\n", __LINE__);}
									if (cutoff[run].skip > 0) {
										fprintf(answers,"\t%u", SKIPS_IN_DMR);
									}
									char gene[GENE_STRING_SIZE] = "n/a", intergene[GENE_STRING_SIZE] = "n/a", list_gene[GENE_STRING_SIZE] = "n/a";
									if (gene_annotation_defined == true) {
										
//										if (DEBUG == true) {printf("@line %u\n", __LINE__);}
//									printf("@ %s line %u list_gene = %s\n", __FILE__, __LINE__, list_gene);
										if (chromosome_annotation_found == true) {identify_gene(mC_start, mC_end, gene, intergene, list_gene, GENE_STRING_SIZE, &min_gene_search_index, &max_gene_search_index, first_nucleotide);}
										fprintf(answers,"\t%s\t%s\t%s", gene, intergene, list_gene);
										gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, gene, false, true);
										gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, intergene, false, true);
										gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, list_gene, false, true);
										if (SUM_PERCENT_DIFFERENCE > 0) {
											methylation_up_gene_list = add_name(methylation_up_gene_list, &methylation_up_gene_list_size, &number_of_methylation_up_list_genes, list_gene, true, false);
										} else {
											methylation_down_gene_list = add_name(methylation_down_gene_list, &methylation_down_gene_list_size, &number_of_methylation_down_list_genes, list_gene, true, false);
										}
										
//										if (DEBUG == true) {printf("@line %u\n", __LINE__);}
									}
									if (options.printing_CpN_list == true) {
										print_CpN_list(answers, CpN_list, good_nucleotide_count);
									}
/*									if (options.isolation == true) {
										fprintf(answers, "\t%u;%u", mC_start - last_CpN_before_DMR, get_next_CpN(mC_end, (const size_t **restrict)bisulfite_nucleotide_min, (const size_t **restrict)bisulfite_nucleotide_max, BISULFITE_DATA, NUMBER_OF_REPLICATES_COMPARED, COMPARED_SETS) - mC_end);
									}*/
									fprintf(answers,"\n");
									if (options.make_figures == true) {
										fprintf(rscript, "setEPS()\n");
										fprintf(rscript, "postscript('DMR_figures/%s_%s_c%u_CpN%u_d%u_p0_%02.0f_P%.0f_S%u_chr%s_%u_%u", label_set1, label_set2, cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p*100, 100.0*cutoff[run].percent, cutoff[run].skip, bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start, mC_end);
										if (strcmp(list_gene,"n/a") != 0) {
											fprintf(rscript, "_");
											fprint_for_filename(rscript, list_gene);
										}
										fprintf(rscript, ".eps')\n");
										
										R_figure(rscript, all_CpN, NO_COAT_TAILS, label_set1, label_set2, xaxis_label, list_gene, all_p, all_percent, all_sum_p, all_sum_percent, all_nucleotides, NUMBER_OF_REPLICATES_COMPARED, mC_start, mC_end);
										all_CpN = uint_clear(all_CpN);
										all_p = double_clear(all_p);
										all_percent = double_clear(all_percent);
										all_sum_p = double_clear(all_sum_p);
										all_sum_percent = double_clear(all_sum_percent);
										all_nucleotides = clear_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED);
									}
									number_of_DMRs++;
								}
							}
//							if (DEBUG == true) {printf("@line %u\n", __LINE__);}
							zero_DMR_struct(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);
							zero_DMR_struct(DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
//							last_CpN_before_DMR = nucleotide;
							mC_count = 0;
							similar_nucleotide_count = 0;
							good_nucleotide_count = 0;
							consecutive_skipped_CpN = 0;
							CpN_skips_in_this_DMR = 0;
							tail_count = 0;
							tail_skips = 0;
							if (options.printing_CpN_list == true) {
								CpN_list = uint_clear(CpN_list);
							}
//							if (DEBUG == true) {printf("@line %u\n", __LINE__);}
							continue;
						} else if ((mC_count > 0) && (consecutive_skipped_CpN < cutoff[run].skip) && (replicates_meeting_minimum_coverage < TOTAL_REPLICATES)) {
//if mC_count > 0, then I don't start on similar CpN
							add_2nd_dmr_to_1st_dmr(DMR_similar_sum, single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
							add_2nd_dmr_to_1st_dmr(DMR_sum, single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
							consecutive_skipped_CpN++;
							CpN_skips_in_this_DMR++;
							tail_count++;
							tail_skips++;
							mC_count++;
							if (options.make_figures == true) {//every CpN counted by mC_count needs an element in these arrays, R_figure reads all of them
								all_CpN = add_uint(all_CpN, nucleotide, mC_count - 1);
								all_p = add_double(all_p, Pvalue(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED, options.paired), mC_count - 1);
								all_percent = add_double(all_percent, 100.0*mean_percent_difference(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED), mC_count - 1);
								all_sum_p = add_double(all_sum_p, Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired), mC_count - 1);
								all_sum_percent = add_double(all_sum_percent, 100.0*mean_percent_difference(DMR_sum, NUMBER_OF_REPLICATES_COMPARED), mC_count - 1);
								all_nucleotides = add_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED, mC_count - 1, single_nucleotide);
							}
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t%u\t%u\t%u\t%.2f\t%g\n", mC_count, similar_nucleotide_count, consecutive_skipped_CpN, 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
							continue;
						} else if (replicates_meeting_minimum_coverage < TOTAL_REPLICATES) {
//skips are allowed, but there isn't a DMR to continue (mC_count = 0), and a nucleotide without minimum coverage in every replicate can't start one
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t0\t0\t0\t%.2f\t%.3g\n", 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
							continue;
						}
//						if (DEBUG == true) {printf("@line %u\n", __LINE__);}
						nucleotides_good_in_every_organism++;
	//					printf("@line %u, single_nucleotide: \n", __LINE__);
	//					print_DMR_struct(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
	//					printf("@line %u, DMR_sum: \n", __LINE__);
//						print_DMR_struct(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);
						add_2nd_dmr_to_1st_dmr(DMR_sum, single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
						const double SINGLE_NUCLEOTIDE_MEAN_DIFFERENCE = mean_percent_difference(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
						mC_count++;
						if (options.make_figures == true) {
							all_CpN = add_uint(all_CpN, nucleotide, mC_count - 1);
							double pvalue = Pvalue(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED, options.paired);
//							if (pvalue < 0) {
	//							printf("@ %s line %u, p < 0\n", __FILE__, __LINE__);
		//					}
							all_p = add_double(all_p, pvalue, mC_count - 1);
							pvalue = Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired);
//							if (pvalue < 0) {
//								printf("@ %s line %u, p < 0\n", __FILE__, __LINE__);
//							}
							all_percent = add_double(all_percent, 100.0*mean_percent_difference(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED), mC_count - 1);
							all_sum_p = add_double(all_sum_p, pvalue , mC_count - 1);
							all_sum_percent = add_double(all_sum_percent, 100.0*mean_percent_difference(DMR_sum, NUMBER_OF_REPLICATES_COMPARED), mC_count - 1);
							all_nucleotides = add_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED, mC_count - 1, single_nucleotide);
						}
						
//						if (DEBUG == true) {printf("@line %u\n", __LINE__);}
	//					printf("Mean Difference = %lf\n",MEAN_DIFFERENCE);
						if ((fabs(SINGLE_NUCLEOTIDE_MEAN_DIFFERENCE) < cutoff[run].percent) || ((options.fdr_on_all_cpn == false) && (Pvalue(single_nucleotide, NUMBER_OF_REPLICATES_COMPARED, options.paired)) > cutoff[run].p) || ((options.fdr_on_all_cpn == true) && (qvalues[total_nucleotides-1] > cutoff[run].p))) {//i.e. this nucleotide is not differentially methylated
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t0\t%u\t%u\t%.2f\t%g\n", similar_nucleotide_count, consecutive_skipped_CpN, 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
							similar_nucleotide_count++;
							tail_count++;
							consecutive_skipped_CpN = 0;
							add_2nd_dmr_to_1st_dmr(DMR_similar_sum, single_nucleotide, NUMBER_OF_REPLICATES_COMPARED);
						} else {//this nucleotide *is* differentially methylated
							similar_nucleotide_count = 0;
							tail_count = 0;
							tail_skips = 0;
							consecutive_skipped_CpN = 0;
	//printf("@line %u\n",__LINE__);
							zero_DMR_struct(DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
							if (mC_count == 1) {
								mC_start = nucleotide;
							}
							mC_end = nucleotide;
							if (options.printing_CpN_list == true) {
								CpN_list = add_uint(CpN_list, nucleotide, good_nucleotide_count);
							}
							good_nucleotide_count++;
							if (options.print_all_nucleotides == true) {
								fprintf(fh_print_all_nucleotides,"\t%u\t0\t0\t%.2f\t%.3g\n", mC_count, 100.0*mean_percent_difference(DMR_sum,NUMBER_OF_REPLICATES_COMPARED), Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired));
							}
							continue;//this nucleotide is differentially methylated, can only improve score, why calculate p-value of set?
						}
	//					printf("chr.%s:%u @ line %u\n", gene_chromosome_changes[chr].chromosome, nucleotide, __LINE__);
//-----------------------------------------------------------
//CASE II: trial array shows differential methylation, so just add to the overall array
//-----------------------------------------------------------
//now I decide if I've got an end region or I'm continuing within a DMR
						const double SUM_MEAN_DIFFERENCE = mean_percent_difference(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);
						
//						if (DEBUG == true) {printf("@line %u\n", __LINE__);}
						const double P_VALUE = fabs(SUM_MEAN_DIFFERENCE) >= cutoff[run].percent ? Pvalue(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, options.paired) : 1.0;
		//				print_2d_double_array(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);//debugging
	//					printf("mean difference = %lf\tp = %lf\n\n",MEAN_DIFFERENCE, P_VALUE);//debugging
						if (P_VALUE < cutoff[run].p) {
							continue;//no point in doing the rest“
						}
//-----------------------------------------------------------
//CASE III: Differential methylation stops
//-----------------------------------------------------------
						const unsigned int NO_COAT_TAILS = mC_count - tail_count;
						const unsigned int SKIPS_IN_DMR = CpN_skips_in_this_DMR - tail_skips;
						if ((NO_COAT_TAILS >= cutoff[run].CpN) && ((mC_end-mC_start) >= options.MINIMUM_CPN_RANGE) && (good_nucleotide_count >= cutoff[run].d)) {
							subtract_2nd_dmr_from_1st_dmr(DMR_sum, DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
							const double SUM_PERCENT_DIFFERENCE = mean_percent_difference(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);
							if (fabs(SUM_PERCENT_DIFFERENCE) >= cutoff[run].percent) {
								if (options.make_bed_files == true) {
									fprintf(bedfile_fh,"chr%s	%u	%u\n",bs_data_chromosome_changes[set1][0][chromosome].chromosome,mC_start ,mC_end);
								}
								fprintf(answers,"%s	%u	%u\t%u	%u", bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start,mC_end, NO_COAT_TAILS, good_nucleotide_count);
								fprintf(answers,"\t%.1f", 100.0*SUM_PERCENT_DIFFERENCE);

								print_methylation_percentages(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, answers);
								if (options.make_heatmap == true) {
									print_methylation_heatmap(DMR_sum, NUMBER_OF_REPLICATES_COMPARED, fh_heatmap, number_of_DMRs);
								}
								if (cutoff[run].skip > 0) {
									fprintf(answers,"\t%u", SKIPS_IN_DMR);
								}
								char gene[GENE_STRING_SIZE] = "n/a", intergene[GENE_STRING_SIZE] = "n/a", list_gene[GENE_STRING_SIZE] = "n/a";
								if (gene_annotation_defined == true) {
									if (chromosome_annotation_found == true) {identify_gene(mC_start, mC_end, gene, intergene, list_gene, GENE_STRING_SIZE, &min_gene_search_index, &max_gene_search_index, first_nucleotide);}
									if (SUM_PERCENT_DIFFERENCE > 0) {
										methylation_up_gene_list = add_name(methylation_up_gene_list, &methylation_up_gene_list_size, &number_of_methylation_up_list_genes, list_gene, true, false);
									} else {
										methylation_down_gene_list = add_name(methylation_down_gene_list, &methylation_down_gene_list_size, &number_of_methylation_down_list_genes, list_gene, true, false);
									}
									fprintf(answers,"\t%s\t%s\t%s", gene, intergene, list_gene);
									gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, gene, false, true);
									gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, intergene, false, true);
									gene_list = add_name(gene_list, &gene_list_size, &number_of_list_genes, list_gene, false, true);
								}
								if (options.make_figures == true) {
									fprintf(rscript, "setEPS()\n");
									fprintf(rscript, "postscript('DMR_figures/%s_%s_c%u_CpN%u_d%u_p0_%02.0f_P%.0f_S%u_chr%s_%u_%u", label_set1, label_set2, cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p*100, 100.0*cutoff[run].percent, cutoff[run].skip, bs_data_chromosome_changes[set1][0][chromosome].chromosome, mC_start, mC_end);
									if (strcmp(list_gene,"n/a") != 0) {
										fprintf(rscript, "_");
											fprint_for_filename(rscript, list_gene);
									}
									fprintf(rscript, ".eps')\n");
									R_figure(rscript, all_CpN, NO_COAT_TAILS, label_set1, label_set2, xaxis_label, list_gene, all_p, all_percent, all_sum_p, all_sum_percent, all_nucleotides, NUMBER_OF_REPLICATES_COMPARED, mC_start, mC_end);
									all_CpN = uint_clear(all_CpN);
									all_p = double_clear(all_p);
									all_percent = double_clear(all_percent);
									all_sum_p = double_clear(all_sum_p);
									all_sum_percent = double_clear(all_sum_percent);
									all_nucleotides = clear_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED);
								}
								if (options.printing_CpN_list == true) {
									print_CpN_list(answers, CpN_list, good_nucleotide_count);
								}
								if (options.print_DMR_pvalue == true) {
									all_DMRs = add_nucleotides(all_DMRs, NUMBER_OF_REPLICATES_COMPARED, number_of_DMRs, DMR_sum);
									if (all_DMRs[0][0][0].coverage == 0) {
										printf("got coverage = 0 at %s line %u\n", __FILE__, __LINE__);
										exit(EXIT_FAILURE);
									}
									DMR_chromosome_indices = add_uint(DMR_chromosome_indices, chromosome, number_of_DMRs);
									DMR_starts             = add_uint(DMR_starts, mC_start, number_of_DMRs);
									DMR_ends               = add_uint(DMR_ends, mC_end, number_of_DMRs);
									nCpN                   = add_uint(nCpN, NO_COAT_TAILS, number_of_DMRs);
									diffCpN               = add_uint(diffCpN, good_nucleotide_count, number_of_DMRs);
									CpN_skips_in_this_DMR_all  = add_uint(CpN_skips_in_this_DMR_all, SKIPS_IN_DMR, number_of_DMRs);
									if (DEBUG == true) {printf("@line %u with %u DMRs\n", __LINE__, number_of_DMRs);}
								}
								number_of_DMRs++;
/*								if (options.isolation == true) {
									fprintf(answers, "\t%u;%u", mC_start - last_CpN_before_DMR, get_next_CpN(mC_end,(const size_t **restrict) bisulfite_nucleotide_min, (const size_t **restrict)bisulfite_nucleotide_max, BISULFITE_DATA, NUMBER_OF_REPLICATES_COMPARED, COMPARED_SETS) - mC_end);
								}*/
								fprintf(answers,"\n");
							}
						}
						mC_count = 0;
//						last_CpN_before_DMR = nucleotide;
						similar_nucleotide_count = 0;
						good_nucleotide_count = 0;
						consecutive_skipped_CpN = 0;
						CpN_skips_in_this_DMR = 0;
							tail_count = 0;
							tail_skips = 0;
						if (options.printing_CpN_list == true) {
							CpN_list = uint_clear(CpN_list);
						}
						zero_DMR_struct(DMR_sum, NUMBER_OF_REPLICATES_COMPARED);
						zero_DMR_struct(DMR_similar_sum, NUMBER_OF_REPLICATES_COMPARED);
					}//end of while loop
					CpN_list = uint_clear(CpN_list);
					all_CpN = uint_clear(all_CpN);
					free(all_CpN); all_CpN = NULL;
					free(CpN_list); CpN_list = NULL;
					free(all_p); all_p = NULL;
					free(all_percent); all_percent = NULL;
					free(all_sum_p); all_sum_p = NULL;
					free(all_sum_percent); all_sum_percent = NULL;
					for (unsigned short int set = 0; set < 2; set++) {//there are only 2 sets
						free(DMR_sum[set]); DMR_sum[set] = NULL;
						free(bisulfite_nucleotide_min[set]); bisulfite_nucleotide_min[set] = NULL;
						free(bisulfite_nucleotide_max[set]); bisulfite_nucleotide_max[set] = NULL;
						free(single_nucleotide[set]); single_nucleotide[set] = NULL;
						free(DMR_similar_sum[set]); DMR_similar_sum[set] = NULL;
					}
					free(DMR_sum);								DMR_sum = NULL;
					free(bisulfite_nucleotide_min);		bisulfite_nucleotide_min = NULL;
					free(bisulfite_nucleotide_max);		bisulfite_nucleotide_max = NULL;
					free(single_nucleotide);				single_nucleotide = NULL;
					free(DMR_similar_sum);					DMR_similar_sum = NULL;
					if (all_nucleotides != NULL) {
						all_nucleotides = final_free_nucleotides(all_nucleotides, NUMBER_OF_REPLICATES_COMPARED);
					}
					if (DEBUG == true) {
						printf("just finished chromosome %s @ line %u\n", bs_data_chromosome_changes[set1][0][chromosome].chromosome, __LINE__);
					}
				}//end of chromosome
				qvalues = double_clear(qvalues);//will eat up RAM, free ASAP
				pvalues = double_clear(pvalues);
				free(qvalues); qvalues = NULL;
				free(pvalues); pvalues = NULL;
				if (options.make_heatmap == true) {
					fclose(fh_heatmap);
				}
				if (number_of_DMRs == 0) {
					fprintf(answers, "\nNo_DMR_was_found_with_current_parameters\n\n");
					puts("\n\nNo_DMR_was_found_with_current_parameters\n");
				}
				fclose(answers);
				if (options.make_bed_files == true) {
					fclose(bedfile_fh);
				}
				if (options.print_all_nucleotides == true) {
					fclose(fh_print_all_nucleotides);
				}
				printf("Answers written to %s\n", answers_filename);//, answers_filename);helpful to run the command: sort -k4n %s\n
				if (options.print_DMR_pvalue == true) {
					puts(answers_filename);
					print_DMRs(answers_filename, all_DMRs, set1, DMR_chromosome_indices, bs_data_chromosome_changes, DMR_starts, DMR_ends, number_of_DMRs, NUMBER_OF_REPLICATES_COMPARED, FDR_type, label_set1, label_set2, gene_annotation_defined, promoter_cutoff, cutoff[run].skip, options.printing_CpN_list, nCpN, diffCpN, CpN_skips_in_this_DMR_all, gene_list, options.paired);
				all_DMRs = final_free_nucleotides(all_DMRs, NUMBER_OF_REPLICATES_COMPARED);
				}
				free(answers_filename); answers_filename = NULL;
				free(DMR_chromosome_indices); DMR_chromosome_indices = NULL;//only used if DMR FDR is 'true'
				free(DMR_starts); DMR_starts = NULL;
				free(DMR_ends);   DMR_ends = NULL;
				free(nCpN);       nCpN = NULL;
				free(diffCpN);    diffCpN = NULL;
				free(CpN_skips_in_this_DMR_all);   CpN_skips_in_this_DMR_all   = NULL;
				gene_list = clear_names(gene_list, &gene_list_size, &number_of_list_genes);
				free(gene_list); gene_list = NULL;
				if ((options.make_bed_files == true) && (number_of_DMRs > 0)) {
					printf("A bed file has been written to %s\n",bedfile_filename);
				}
				free(bedfile_filename); bedfile_filename = NULL;
				if (options.make_figures == true) {
					fclose(rscript);
					char command[64] = "Rscript ";
					strcat(command, rscript_filename);
					if (system(command) != 0) {
						fprintf(stderr, "%s FAILED.\n", command);
						perror("");
						exit(EXIT_FAILURE);
					}
				}
				if ((options.make_heatmap == true) && (number_of_DMRs > 2)) {
					make_heatmap_eps(heatmap_tsv_filename, label_set1, label_set2, file_suffix);
				}
				free(heatmap_tsv_filename); heatmap_tsv_filename = NULL;
//-------------------------------------------------
//finished opening gene list file
//-------------------------------------------------
				if (gene_annotation_defined == true) {
					if (number_of_methylation_up_list_genes > 0) {//only print the names if there are names to print
						FILE *restrict methylation_up_gene_list_fh = NULL;
						char *restrict gene_list_filename = checked_malloc(strlen(label_set1)+strlen(label_set2)+sizeof(char*)+strlen("_vs____methylation_up_gene_list.txt")+strlen(file_suffix), __LINE__);
						sprintf(gene_list_filename,"%s_vs_%s",label_set1, label_set2);
						strcat(gene_list_filename,"_");
						strcat(gene_list_filename,file_suffix);
						strcat(gene_list_filename,"_methylation_up_gene_list.txt");
						methylation_up_gene_list_fh = fopen(gene_list_filename,"w");
						if (methylation_up_gene_list_fh == NULL) {
							perror("");
							fprintf(stderr, "Failed to write %s\n", gene_list_filename);
							exit(EXIT_FAILURE);
						}
						for (unsigned int gene_index = 0; gene_index < number_of_methylation_up_list_genes; gene_index++) {
							fprintf(methylation_up_gene_list_fh,"%s\n",methylation_up_gene_list[gene_index]);
						}
						fclose(methylation_up_gene_list_fh);
						printf("A gene list, with intergenic genes identified with promoter cutoff = %u nucleotides has been written to %s\n", promoter_cutoff, gene_list_filename);
						free(gene_list_filename); gene_list_filename = NULL;
					} else {
						printf("For %s vs. %s, there were no detected genes with increases in methylation\n", label_set1, label_set2);
					}
					if (number_of_methylation_down_list_genes > 0) {//only print the names if there are names to print
						FILE *restrict gene_list_fh = NULL;
						char *restrict gene_list_filename = checked_malloc(strlen(label_set1)+strlen(label_set2)+sizeof(char*)+strlen("_vs____methylation_down_gene_list.txt")+strlen(file_suffix), __LINE__);
						sprintf(gene_list_filename,"%s_vs_%s",label_set1, label_set2);
						strcat(gene_list_filename,"_");
						strcat(gene_list_filename,file_suffix);
						strcat(gene_list_filename,"_methylation_down_gene_list.txt");
						gene_list_fh = fopen(gene_list_filename,"w");
						if (gene_list_fh == NULL) {
							perror("");
							fprintf(stderr, "Failed to write %s\n", gene_list_filename);
							exit(EXIT_FAILURE);
						}
						for (unsigned int gene_index = 0; gene_index < number_of_methylation_down_list_genes; gene_index++) {
							fprintf(gene_list_fh,"%s\n",methylation_down_gene_list[gene_index]);
						}
						fclose(gene_list_fh);
						printf("A gene list, with intergenic genes identified with promoter cutoff = %u nucleotides has been written to %s\n", promoter_cutoff, gene_list_filename);
						free(gene_list_filename); gene_list_filename = NULL;
					} else {
						printf("For %s vs. %s, there were no detected genes with decreases in methylation\n", label_set1, label_set2);
					}
				}
//begin clearing out memory
				methylation_up_gene_list = clear_names(methylation_up_gene_list, &methylation_up_gene_list_size, &number_of_methylation_up_list_genes);
				free(methylation_up_gene_list); methylation_up_gene_list = NULL;
				methylation_down_gene_list = clear_names(methylation_down_gene_list, &methylation_down_gene_list_size, &number_of_methylation_down_list_genes);
				free(methylation_down_gene_list); methylation_down_gene_list = NULL;

				printf("\t%u/%u = %.1f%% of nucleotides met minimum coverage and were present with minimum coverage = %u in every sample.\n", nucleotides_good_in_every_organism, total_nucleotides, total_nucleotides > 0 ? 100.0*(double)nucleotides_good_in_every_organism /(double)total_nucleotides : 0.0, cutoff[run].coverage);
				free(file_suffix); file_suffix = NULL;
				cutoff[run].number_of_DMRs = number_of_DMRs;
			} // end each run
			if (RUNS > 1) {	
//first need to allocate memory for the filename string
				size_t TOTAL_LABEL_STRING_LENGTH = 0;
				if (total_label != NULL) {//strlen on null string gives segfault
					TOTAL_LABEL_STRING_LENGTH = strlen(total_label);
				}
				char *restrict dmr_count_filename = checked_malloc(strlen(label_set1) + strlen(label_set2) + strlen("_vs___dmr_count.tsv") + TOTAL_LABEL_STRING_LENGTH + 1, __LINE__);//this used sizeof(label_set1), which is the size of a pointer, not the length of the label
				if (TOTAL_LABEL_STRING_LENGTH > 0) {
					sprintf(dmr_count_filename,"%s_vs_%s_%s_dmr_count.tsv", label_set1, label_set2, total_label);//make the filename
				} else {
					sprintf(dmr_count_filename,"%s_vs_%s_dmr_count.tsv", label_set1, label_set2);//make the filename
				}
				FILE *restrict dmr_count_fp = fopen(dmr_count_filename, "w");//open file for writing
				if (dmr_count_fp == NULL) {//in case I can't write the file, abort the execution
					free(label_set1); label_set1 = NULL;
					free(label_set2); label_set2 = NULL;
					fprintf(stderr, "Couldn't write %s\n", dmr_count_filename);
					perror("");
					exit(EXIT_FAILURE);
				}
				fprintf(dmr_count_fp, "Coverage\tMin.CpN\tMin.Diff.CpN\tMax.pValue\tMin.percent\tSkips\t#DMRs\n");//print the header
				for (unsigned int run = 0; run < RUNS; run++) {
					fprintf(dmr_count_fp, "%u\t%u\t%u\t%g\t%lf\t%u\t%u\n", cutoff[run].coverage, cutoff[run].CpN, cutoff[run].d, cutoff[run].p, 100.0*cutoff[run].percent, cutoff[run].skip, cutoff[run].number_of_DMRs);
				}
				fclose(dmr_count_fp);//close the file pointer
				printf("DMR counts with varying cutoffs have been written to %s\n", dmr_count_filename);
				free(dmr_count_filename); dmr_count_filename = NULL;//free the memory used for the dmr_count_filename
			}
			free(cutoff); cutoff = NULL;
			free(label_set1); label_set1 = NULL;
			free(label_set2); label_set2 = NULL;
		}
	}
//---------------------------------------------------------------------------------------------
//Begin Freeing memory
//---------------------------------------------------------------------------------------------
	free(FDR_type); FDR_type = NULL;
	for (unsigned short int set_loop = 0; set_loop < SETS; set_loop++) {
		for (unsigned short int replicate_loop = 0; replicate_loop < number_of_replicates_per_set[set_loop]; replicate_loop++) {
			if (DEBUG == true) {
		 		printf("Now freeing BISULFITE_DATA[%u][%u] on line %u\n",set_loop, replicate_loop, __LINE__);
		 	}
			free(BISULFITE_DATA[set_loop][replicate_loop]); BISULFITE_DATA[set_loop][replicate_loop] = NULL;
			if (DEBUG == true) {
		 		printf("Now freeing bs_data_chromosome_changes[%u][%u] on line %u\n",set_loop, replicate_loop, __LINE__);
		 	}
		 	for (unsigned int chr_loop = 0; chr_loop <= number_of_data_chromosomes[set_loop][replicate_loop]; chr_loop++) {
//		 		printf("freeing bs_data_chromosome_changes[%u][%u][%u].chromosome\n", set_loop, replicate_loop, chr_loop);
		 		free(bs_data_chromosome_changes[set_loop][replicate_loop][chr_loop].chromosome);
		 		bs_data_chromosome_changes[set_loop][replicate_loop][chr_loop].chromosome = NULL;
		 	}
			free(bs_data_chromosome_changes[set_loop][replicate_loop]); bs_data_chromosome_changes[set_loop][replicate_loop] = NULL;
		}
		if (DEBUG == true) {
			printf("about to free BISULFITE_LOOP[%u] on line %u\n",set_loop,__LINE__);
		}
		free(BISULFITE_DATA[set_loop]); BISULFITE_DATA[set_loop] = NULL;
		if (DEBUG == true) {
			printf("about to free bs_data_chromosome_changes[%u] on line %u\n",set_loop,__LINE__);
		}
		free(bs_data_chromosome_changes[set_loop]); bs_data_chromosome_changes[set_loop] = NULL;
		free(number_of_data_chromosomes[set_loop]); number_of_data_chromosomes[set_loop] = NULL;
	}
	free(BISULFITE_DATA); BISULFITE_DATA = NULL;
	free(number_of_data_chromosomes); number_of_data_chromosomes = NULL;
	for (unsigned int chr_loop = 0; chr_loop < number_of_chromosomes_in_GENES; chr_loop++) {
		free(gene_chromosome_changes[chr_loop].chromosome); gene_chromosome_changes[chr_loop].chromosome = NULL;
	}
	free(gene_chromosome_changes); gene_chromosome_changes = NULL;
	free(number_of_replicates_per_set); number_of_replicates_per_set = NULL;
	free(bs_data_chromosome_changes); bs_data_chromosome_changes = NULL;
	for (unsigned int g = 0; g < number_of_annotated_genes; g++) {
		free(GENES[g].gene_name); GENES[g].gene_name = NULL;
	}
	free(GENES); GENES = NULL;
	labels = clear_names(labels, &labels_size, &number_of_labels);
	free(labels); labels = NULL;
	free(xaxis_label); xaxis_label = NULL;
	free(total_label); total_label = NULL;
	puts("Finished analysis.");
	return 0;
}
