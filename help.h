#ifndef DEFIANT_HELP_H
#define DEFIANT_HELP_H
#include <stdio.h>
#include <limits.h>

void help (void) {
	puts("Help for DEFIANT\n");
	printf("%c[1mSynopsis%c[0m\n", 27, 27);
	puts("./defiant [OPTIONS]... [FILES]...\n\nwhere [FILES] is a list of files with spaces separating groups and commas separating replicates within groups. Options are like \"-a refFlat.gtf\" etc. (cf. Options).");
	puts("\nMinimal example: ./defiant -i control1,control2 case1,case2\n");
	puts("All files in should be specified after \"-i\", where commas separate replicates in the same group.\n");
	printf("%c[1mOptions%c[0m\n\n", 27, 27);//"%c[1m" starts a bold face, and "%c[0m" ends bold face, and must have 27 on the right end of the printf statement
	printf("%c[1m-a%c[0m\tSpecify annotation file, e.g. \"-a mm10.gtf\"\n\n", 27, 27);
	printf("%c[1m-b%c[0m\tOutput DMRs in bed file. This option does not take an argument.\n\n", 27, 27);
	printf("%c[1m-c%c[0m\tminimum coverage, e.g. \"-c 10\".  This option accepts positive integers and can be parallelized to test multiple options.\n\n", 27, 27);
	printf("%c[1m-CpN%c[0m\tminimum CpN/CpG/CH/CHH in a DMR, e.g. \"-CpN 10\".  This option accepts positive integers and can parallelized.  \"CpN\" is case insensitive.\n\n", 27, 27);
	printf("%c[1m-cpu%c[0m\tSet number of CPU when running multiple options, e.g. \"-cpu 4\". \"CPU\" is case insensitive and accepts integers > 0.  Runs are executed in parallel only if Defiant was compiled with `make OPENMP=1`, otherwise one after another; with more than 1 CPU the default value of every varied option is added to the runs.\n\n", 27, 27);
	printf("%c[1m-d%c[0m\tMinimum differential nucleotide count in a DMR, e.g. \"-d 3\".  This option can be parallelized.\n\n", 27, 27);
	printf("%c[1m-D%c[0m\tMaximum non-default options in a parallel run, e.g. \"-D 4\"\n\n", 27, 27);
	printf("%c[1m-debug%c[0m\tTurn on debugging mode.  This slows down the execution significantly, but can help diagnose problems if they arise.  This option does not accept any arguments.\n\n", 27, 27);
	printf("%c[1m-E%c[0m\tprint statistics for every CpN. This option does not take an argument. This slows Defiant down significantly.\n\n", 27,27);
	printf("%c[1m-f%c[0m\tmake EPS figures for each DMR. Warning: requires R installation. This option does not take an argument, and will slow defiant's execution.\n\n", 27, 27);
	printf("%c[1m-fdr%c[0m\tCalculate FDR-adjusted q-value for each CpN.  'FDR' is case insensitive.  This option can take case-insensitive arguments 'fdr' or 'bh' for Benjamini-Hothberg method, 'Bonferroni', 'Hochberg', 'Hommel', 'Holm', or 'BY' for Benjamini & Yekutieli.  If no argument is given, 'BH' is assumed.  This function is a translation of R's 'p.adjust'.  I recommend against using this option as for genome-scale CpG measurements, almost everything will be q = 1 and no DMRs will be obtained in any case.  This option will substantially increase RAM use and slow execution.  'Hommel' is so slow I strongly recommend against it.\n\n", 27, 27);
	printf("%c[1m-G%c[0m\tMaximum allowed gap between CpN, e.g. \"-G 1000\"\n\n", 27, 27);
	printf("%c[1m-h%c[0m\tPrint this help menu; \"-h\" is case insensitive. This option does not take an argument, and defiant exits after this option is read.\n\n", 27, 27);
	printf("%c[1m-heatmap%c[0m\tMake heatmap.  Requires R be available from command line, and that `library(gplots)` will work.\n\n", 27, 27);
	printf("%c[1m-i%c[0m\tStart reading input files.  This is the only required argument.  All further entries to the command line are assumed to be files.\n\n", 27, 27);
//	printf("%c[1m-I%c[0m\tPrint isolation of nucleotides in output files. This option does not take an argument.\n\n", 27, 27);
	printf("%c[1m-l%c[0m\tSet output file(s) label, e.g. \"-l new\"\n\n", 27, 27);
	printf("%c[1m-L%c[0m\tgive labels for each set in a comma-delimited string, e.g. \"-L case,control\"\n\n", 27, 27);
	printf("%c[1m-N%c[0m\tlist CpG Nucleotides in the DMR in output file. This option does not take an argument.\n\n", 27, 27);
//	printf("%c[1m-o%c[0m\toverwrite existing files if present.  This option does not take an argument.\n\n", 27, 27);
	printf("%c[1m-p%c[0m\tMaximum p-value, which is 0<=p<=1.  This option can be parallelized to test multiple options.  Default 0.05.\n\n", 27, 27);
	printf("%c[1m-paired%c[0m\t--paired. P-values will be paired.\n", 27, 27);
	printf("%c[1m-P%c[0m\tMinimum Percent methylation difference (0 <= P <= 100). This option can be parallelized to test multiple options (default 10%%).\n\n", 27, 27);
	printf("%c[1m-q%c[0m\tPromoter cutoff for gene assignment of intergenic DMRs (default 10,000 nucleotides).  This option accepts positive integers, e.g. \"-q 15000\".\n\n", 27, 27);
	printf("%c[1m-r%c[0m\tMinimum nucleotide range, which accepts a non-negative integer.  Default range is 0 nucleotides.\n\n", 27, 27);
	printf("%c[1m-R%c[0m\tinclude \"Random\" chromosomes.  This option does not accept an argument.\n\n", 27, 27);
	printf("%c[1m-s%c[0m\tMaximum allowed consecutive similar CpN, default is 6 CpN.  This accepts non-negative integers, e.g. \"-s 3\".\n\n", 27, 27);
	printf("%c[1m-S%c[0m\tAllow some number of consecutive skips of low coverage, default is 0.  This accepts positive integers, e.g. \"-S 1\", or \"-S inf\" for no limit.\n\n", 27, 27);
	printf("%c[1m-U%c[0m\tInclude \"Un\" chromosomes (default is to ignore them).  This option does not accept an argument.\n\n", 27, 27);
	printf("%c[1m-v%c[0m\tPrint a p-value & FDR-adjusted p-value for each DMR. This option accepts the same arguments that the '-FDR' option does.\n\n", 27, 27);
		printf("%c[1m-version%c[0m\tPrint Defiant's version.\n\n", 27, 27);
//	printf("%c[1m-w%c[0m\tApproximate proportion of CpN p-value products that are less than this DMR p-value product. This option doesn't accept an argument.\n\n", 27, 27);
	printf("%c[1m-x%c[0m\tThis option accepts a string which will become the x-axis in figures.  \"-x\" activates \"-f\" option and requires an R installation.\n\n", 27, 27);
//--------------------------------------------
	printf("%c[1mParallelization%c[0m\n\n", 27, 27);//"%c[1m" starts a bold face, and "%c[0m" ends bold face, and must have 27 on the right end of the printf statement
	puts("As each experiment is different, a different set of parameters may be appropriate for each experiment. You may not know these parameters ahead of time. Thus, defiant has been set to easily test multiple parameters in parallel via a shared-memory model. Parameters underlined and bold faced in Table 1, e.g. p, can be written like a C-style for loop and delimited with commas: ./defiant -p < min >, < max >, < step > which would increment the p-value from min to max in steps of step.");
	
	puts("For example, \n\n./defiant -c 5,15,5\n\nwhich would run minimum coverage from 5 to 15 in steps of 5.\nIf you only wish to run two parameters, you can simply write a comma in between the two parameters you wish to vary, e.g. \n");
	puts("./defiant -p 0.01,0.05\n");
	puts("will run p = 0.01 and then p = 0.05.\n");
	puts("The data is read off of the hard drive and into memory, which will then be shared among all the CPU. This is done to make 3D graphs, i.e. x vs. y with the 3rd dimension in color. However, the -D option can be used to vary all parameters as a nested for loop. I strongly recommend not to use the \"-f\" or \"-x\" options with multiple runs. One of defiant’s advantages is speed and low resource use, using both \"-f\" and \"-x\" options will make the runs take much much longer and potentially create a lot of files which will make I/O on your computer very slow. All DMR counts are then saved to a table, which will end in something like dmr_count.tsv\n");
//--------------------------------------------
//--------------- input --------------------------------
//--------------------------------------------
	printf("\n%c[1mInput Formats%c[0m\nDefiant is set up to automatically identify and read the following input formats:\n\n", 27, 27);
	printf("\n%c[1mInput Type 1%c[0m\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1	762	763	0.1764	37\n",27,27);//Jon Schug's format
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tnucleotide, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column3:\tignored.");
	puts("Column4:\tmethylation percent, a floating point in [0,1]");
	printf("Column5:\tcoverage, an unsigned integer, an unsigned integer in [0,%u]\n", UINT_MAX);
//--------------------------------------------
	printf("\n%c[1mInput Type 2%c[0m\tknown for MethylKit input\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1.762	chr1	762	R	100	17.64	82.36\n", 27, 27);
	puts("Column1:\tunique name, this is ignored.");
	puts("Column2:\tchromosome, which is a string.");
	printf("Column3:\tnucleotide, an unsigned integer [0,%u]\n", UINT_MAX);
	puts("Column4:\tsense, this is ignored.");
	printf("Column5:\tcoverage, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column6:\tmethylation percent, a floating point in [0,100].");
	puts("Column7:\tcytosine percent, a floating point in [0,100].");
//--------------------------------------------
	printf("\n%c[1mInput Type 3%c[0m\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1	762	763	0.1764\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tnucleotide, an unsigned integer [0,%u]\n", UINT_MAX);
	puts("Column3:\tignored.");
	puts("Column4:\tmethylation percent, a floating point in [0,1].");
//--------------------------------------------	
	printf("\n%c[1mInput Type 4%c[0m\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1	762	6	14\n",27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tnucleotide, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column3:\tmethylated C count, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column4:\tcount unmethylated, an unsigned integer [0,%u]\n", UINT_MAX);
//--------------------------------------------//return 5:
	printf("\n%c[1mInput Type 5%c[0m\tBismark bismark2bedGraph format:\n<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>(https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1\t762\t+\t17\t64\tCG\tCGA\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tnucleotide position, an unsigned integer [0,%u]\n", UINT_MAX);
	puts("Column3:\tstrand (+ or -)");
	printf("Column4:\tmethylated C count, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column5:\tcount unmethylated, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column6:\tC-context, e.g. CG, CH, CHH.");
	puts("Column7:\tC-context, e.g. CGA, CGT, etc.");
//return 6 <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>//Bismark bismark2bedGraph format
	printf("\n%c[1mInput Type 6%c[0m\tBismark coverage2cytosine format:\n<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1\t762\t763\t0.223684\t17\t76\n", 27, 27);
	puts("Column1:\tchromosome, which is a string");
	printf("Column2:\tnucleotide/start position, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column3:\tnucleotide/end position, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column4:\tmethylation percentage in [0,1]");
	printf("Column5:\tmethylated C count, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column6:\tcount unmethylated, an unsigned integer in [0,%u]\n", UINT_MAX);
//return 7::tid chr pos m-score conf
	printf("\n%c[1mInput Type 7%c[0m\tHELP-Tag data.  This can have a header.\n", 27, 27);
	printf("%c[1mExample%c[0m:\t1 chr1\t762\t763\t0.2656\t0.1776\n", 27, 27);
	puts("Column1:\tignored");
	puts("Column2:\tchromosome, a string");
	printf("Column3:\tposition, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column4:\tmethylation percent: a floating point number in [0,1]");
	puts("Column5:\tConf. ignored.");
//return 8: chr1	3010874	3010876	'46/47'	978	+ (EPP)Epigenome Processing Pipeline
	printf("\n%c[1mInput Type 8%c[0m\t(EPP)Epigenome Processing Pipeline\n",27,27);
	printf("%c[1mExample%c[0m:\tchr1\t762\t763\t'17/76'\t999\t+\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tstart nucleotide, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column3:\tend nucleotide, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column4:\tmethylation percent as a fraction, two unsigned integers.  Coverage is given as the denominator. Everything after column 4 is ignored.");
//return 9: bsmooth input: 10	60025	+	CG	25	26
	printf("\n%c[1mInput Type 9%c[0m\tBsmooth Input\n",27,27);
	printf("%c[1mExample%c[0m:\tX	762	+	CG	17	76\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tNucleotide, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column3:\tstrand sense, ignored");
	puts("Column4:\tcontext, ignored");
	printf("Column5:\tmethylated C count, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column6:\tcoverage, an unsigned integer in [0,%u]\n", UINT_MAX);
//return 10:	chr1	10496	10497	79.69	64	+	10496	10497	180,60,0 (BisSNP), found in RnBeads
	printf("\n%c[1mInput Type 10%c[0m\tBisSNP (found in RnBeads)\n",27,27);
	printf("%c[1mExample%c[0m:\tX\t762\t763\t10.96\t1918\t762\t763\t180,60,0\t0\t0\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	printf("Column2:\tNucleotide start, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column3:\tNucleotide end, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column4:\tmethylation value in [0:100].");
	printf("Column5:\tCoverage, an unsigned integer in [0,%u]. Everything after this column is ignored.\n", UINT_MAX);
//return 11: chr1	G	3000851	CHH	CC	0.1	1	10 (https://github.com/BSSeeker/BSseeker2)
	printf("\n%c[1mInput Type 11%c[0m\t CGmap file from BSSeeker https://github.com/BSSeeker/BSseeker2\n", 27, 27);
	printf("%c[1mExample%c[0m:\tchr1	C	3001631	CG	CG	1.0	5	5\n", 27, 27);
	puts("Column1:\tchromosome, which is a string.");
	puts("Column2:\tnucleotide on Watson (+) strand (ignored)");
	printf("Column3:\tNucleotide end, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column4:\tcontext (CG/CHG/CHH) (ignored)");
	puts("Column5:\tdinucleotide-context (CA/CC/CG/CT) (ignored)");
	puts("Column6:\tmethylation-level = #_of_C / (#_of_C + #_of_T)  (ignored)");
	printf("Column7:\t#_of_C (methylated C, the count of reads showing C here), an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column8:\t#_of_C + #_of_T (all Cytosines, the count of reads showing C or T here), an unsigned integer in [0,%u]\n", UINT_MAX);
// return 12: <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
/*	printf("\n%c[1mInput Type 12%c[0m\t<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context> Bismark coverage2cytosine format\n", 27, 27);
	puts("Column1:	chromosome, doesn't have \"chr\"\n");
	printf("Column2:	position an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column3:	strand (+ or -)");
	printf("Column4:	count methylated, an unsigned integer in [0,%u]\n", UINT_MAX);
	printf("Column5:	count unmethylated, an unsigned integer in [0,%u]\n", UINT_MAX);
	puts("Column6: C context");
	puts("Column7: trinucleotide context");*/
//--------------- output -------------------------
//--------------------------------------------
	printf("\n%c[1mOutput%c[0m\n\n", 27, 27);
	puts("The file names are formatted according to the options set at the command line. For example, consider the output file:");
	puts("\ncontrol_vs_case_c10_CpN5_d1_p0.01_P10.tsv\n");

	puts("• the start of the filename indicates that the \"control\" group was compared against the \"case\" group, as indicated by the \"-L\" option on the command line");
	puts("• minimum coverage of 10,");
	puts("• minimum CpN count = 5,");
	puts("• minimum differentially methylated CpN = 1,");
	puts("• a maximum p-value of 0.01");
	puts("• a minimum Percent difference of 10%.");
	puts("• If the Gap, \"-G\" option, is altered, -G<gap> will be in the output. If there are any allowed consecutive skips, this will be");
	puts("indicated in the filename by \" S2\" if there are 2 allowed skips, for example.");
//---------------running --------------------------
	printf("\n%c[1mRunning Defiant%c[0m\n\n", 27, 27);
	puts("Simplest possible case:");
	puts("1. ./defiant -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt\n");
	puts("Label each sample.");
	puts("2. ./defiant -L control,case -l my_name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt\n");
	puts("3. Annotation:");
	puts("./defiant -a refFlat.txt -L control,case -l my_name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt\n");
	puts("4. Generate figures with \"CpG\" key (requires an installation of the R programming language.)");
	puts("./defiant -x CpG -a refFlat.txt  -L control,case -l my_name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt\n");
}

void version (void) {
	puts("This is Defiant version 1.2");
	puts("Copyright David E. Condon, Ph.D. Wednesday, November 11th, 2026.");
}
#endif
