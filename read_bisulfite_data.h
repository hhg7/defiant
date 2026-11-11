#ifndef DEFIANT_READ_BISULFITE_DATA_H
#define DEFIANT_READ_BISULFITE_DATA_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <regex.h>
//requires safe_str_num.h and annotation.h to be included first

bool regex_match(const char *restrict string, const char *restrict pattern) {
//slightly modified from https://stackoverflow.com/questions/1631450/c-regular-expression-howto#1631458
  regex_t re;
  if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0) return false;
  int status = regexec(&re, string, 0, NULL, 0);
  regfree(&re);
  if (status != 0) return false;//doesn't match
  return true;//matches
}

bool is_string_integer (const char *restrict STRING) {
	if ((STRING == NULL) || (STRING[0] == '\0') || (STRING[0] == '\n')) {
		return false;//an empty string isn't an integer
	}
	for (size_t x = 0; STRING[x] != '\0'; x++) {
		if (STRING[x] == '\n') {
			break;
		}
		if (!isdigit((unsigned char)STRING[x])) {
			return false;
		}//if any digit isn't a digit, return false
	}
	return true;
}

bool is_string_float (const char *restrict STRING) {
	const size_t STRING_LENGTH = strlen(STRING);
	unsigned short int number_of_decimal_points = 0;
	for (size_t x = 0; x < STRING_LENGTH; x++) {
		if (STRING[x] == '\n') {
			break;
		}
		if ((x > 0) && (x < (STRING_LENGTH-1)) && (STRING[x] == '.')) {
			if (number_of_decimal_points == 1) {
				return false;
			}//only allowed 1 decimal point '.'
			number_of_decimal_points++;
			continue;
		}
		if (!isdigit((unsigned char)STRING[x])) {
			return false;
		}//if any digit isn't a digit, return false
	}
	if (number_of_decimal_points == 1) {
		return true;
	} else {
		return false;
	}
}

bool is_string_numeric (const char *restrict STRING) { // float or integer
	if (
	(is_string_integer(STRING) == true)
	||
	(is_string_float(STRING) == true)
	) {
		return true;
	}
	return false;
}

typedef struct {
	unsigned int nucleotide, coverage, methylated_C;
	char strand;//'+', '-', or 0 if the file type doesn't say
	char *restrict chromosome;
} LINE_STRUCT;

char get_strand (const char *restrict STRING) {
	if ((strcmp(STRING, "+") == 0) || (strcmp(STRING, "F") == 0)) {
		return '+';
	} else if ((strcmp(STRING, "-") == 0) || (strcmp(STRING, "R") == 0)) {
		return '-';
	}
	return 0;
}

typedef struct {
	unsigned short int type;
	bool header;
} FILE_TYPE_STRUCT;

bool bool_arrays_match (
const bool *restrict array1, const bool *restrict array2,
const unsigned short MIN_INDEX, const unsigned short MAX_INDEX
) {
	for (unsigned short i = MIN_INDEX; i < MAX_INDEX; i++) {
		if (array1[i] != array2[i]) {
			return false;
		}
	}
	return true;
}

/*void print_bool_array(
const bool *restrict bool_array,
const unsigned short MAX_INDEX
) {
	printf("%s", bool_array[0] ? "true" : "false");
	for (unsigned short i = 1; i < MAX_INDEX; i++) {
		printf(",%s", bool_array[i] ? "true" : "false");
	}
	puts("");
}*/

FILE_TYPE_STRUCT get_file_type (const char *restrict line) {//determines file type of input
//return 1: "chr1	768	768	0.4666	15"//Jon Schug's format
//return 2: "chr21.43008720	chr21	43008720	R	67	92.54	7.46"
//return 2 header "chrBase	chr	base	strand	coverage	freqC	freqT"//some program gives this output, I don't know which one, header
//return 3: "chr1	768	768	0.4666" without coverage information//found on Geobase
//return 4:"chr1	3391	6	14"//found on GeoBase
//return 5:<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>//Bismark coverage2cytosine format
//return 6 <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>//Bismark bismark2bedGraph format
//return 7::tid chr pos m-score conf
//return 8: chr1	3010874	3010876	'46/47'	978	+ (EPP)Epigenome Processing Pipeline
//return 9: bsmooth input: 10	60025	+	CG	25	26
//return 10:	chr1	10496	10497	79.69	64	+	10496	10497	180,60,0 (BisSNP), found in RnBeads
//check list of possible headers
	FILE_TYPE_STRUCT return_struct;
	return_struct.type = 0;
	if (strcmp("chrBase	chr	base	strand	coverage	freqC	freqT\n", line) == 0) {
		return_struct.header = true;
		return_struct.type = 2;
		return return_struct;
	} else if (strcmp(":tid chr pos m-score conf\n", line) == 0) {
		return_struct.header = true;
		return_struct.type = 7;
		return return_struct;
	}
	return_struct.header = false;
	char *restrict line_copy = strdup(line);
	if (line_copy == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	chomp(line_copy);//get rid of "\n" or "\r\n"
	char *restrict temporary_char_string,  *saveptr;//necessary for strtok_r
	temporary_char_string = strtok_r(line_copy, "\t ", &saveptr);
// ZZZ the following arrays record how the input line looks, it will then be matched against known formats
	bool cols_are_chr[]        = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_cg_chg_chh[] = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_cx[]         = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_cxx[]        = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_float[]      = {0,0,0,0,0,0,0,0,0,0,0}; // 0 = false, 1 = true
	bool cols_are_fraction[]   = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_int[]        = {0,0,0,0,0,0,0,0,0,0,0};
	bool col_strand_definit[]  = {0,0,0,0,0,0,0,0,0,0,0};
	bool cols_are_numeric[]		= {0,0,0,0,0,0,0,0,0,0,0};
//these are the definitions of each column
	
	unsigned short no_of_cols = 1;
	unsigned int intcol2 = 0, intcol3 = 0;

	while (temporary_char_string != NULL) {
		if (no_of_cols > 11) {
			printf("There are more columns than I can account for at %s line %u\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
//		printf("determining file type at %s line %u\n", __FILE__, __LINE__);
//		puts(temporary_char_string);
//		cols_are_chr[no_of_cols-1]      = regex_match(temporary_char_string, "chr");
		cols_are_int[no_of_cols-1]      = is_string_integer(temporary_char_string);
		cols_are_cxx[no_of_cols-1]      = regex_match(temporary_char_string, "C[ACGHNT]{2}");
		cols_are_float[no_of_cols-1]    = is_string_float(temporary_char_string);
		cols_are_fraction[no_of_cols-1] = regex_match(temporary_char_string, "[0-9]/[1-9]");
		cols_are_numeric[no_of_cols-1]	= is_string_numeric(temporary_char_string);
		if ((no_of_cols == 2) && (cols_are_int[1] == true)) {
			intcol2 = (unsigned)safe_strtoul(temporary_char_string);
		}
		if ((no_of_cols == 3) && (cols_are_int[2] == true)) {
			intcol3 = (unsigned)safe_strtoul(temporary_char_string);
		}
		if (
		(regex_match(temporary_char_string, "C[ACGT]") == true) && 
		(strlen(temporary_char_string) == 2)) {
			cols_are_cx[no_of_cols-1] = true;
		}
		if (
		(strlen(temporary_char_string) == 3) &&
		(regex_match(temporary_char_string, "C[ACGHNT]{2}") == true)
		) {
			cols_are_cg_chg_chh[no_of_cols-1] = true;
		}
		if ((regex_match(temporary_char_string, "[RF+-]")) && (strlen(temporary_char_string) == 1)) {
			col_strand_definit[no_of_cols-1] = true;
		}
		temporary_char_string = strtok_r(NULL, "\t ", &saveptr);
		if (temporary_char_string == NULL) {
			break;
		}
		no_of_cols++;
	}
	const bool all_cols_false[] = {0,0,0,0,0,0,0,0,0,0,0};
//only write the ones that aren't totally empty
//	const bool type1_chr_cols[]		= {1,0,0,0,0};
	const bool type1_float_cols[]		= {0,0,0,1,0}; // column 4 is either float or int
	const bool type1_int_colsA[]		= {0,1,1,0,1}; // column 4 is either float or int
	const bool type1_int_colsB[]		= {0,1,1,1,1}; // column 4 is either float or int
//chr1.762	chr1	762	R	100	17.64	82.36
//	const bool type2_chr_cols[]		= {1,1,0,0,0,0,0};// two chr columns
	const bool type2_float_cols[]		= {0,0,0,0,0,1,1};
	const bool type2_int_cols[]		= {0,0,1,0,1,0,0};
	const bool type2_strand_cols[]	= {0,0,0,1,0,0,0};
//chr1	762	763	0.1764
//	const bool type3_chr_cols[]		= {1,0,0,0};
//chr1	762	6	14
	const bool type4_int_cols[]		= {0,1,1,1};
//chr1	762	+	17	64	CG	CGA
//	const bool type5_chr_cols[]		= {1,0,0,0,0,0,0};
	const bool type5_int_cols[]		= {0,1,0,1,1,0,0};
	const bool type5_strand_cols[]	= {0,0,1,0,0,0,0};
	const bool type5_cg_cols[]			= {0,0,0,0,0,1,0};
	const bool type5_cxx_cols[]		= {0,0,0,0,0,0,1};
//chr1	762	763	0.223684	17	76
//	const bool type6_chr_cols[]		= {1,0,0,0,0,0};
	const bool type6_int_cols[]		= {0,1,1,0,1,1};
	const bool type6_numeric_cols[]	= {0,1,1,1,1,1};
//20 chr1 43924 0 0.5 0.5	e.g. rn5.Control-1_Hpa.AC56EJACXX.lane_8.hcount
//	const bool type7_chr_cols[]		= {0,1,0,0,0};
	const bool type7_numeric_cols[]	= {1,0,1,1,1};
	const bool type7_int_cols[]		= {1,0,1,1}; //only test to 4 cols for integers
// chr1	762	763	'17/76'	999	+
	const bool type8_int_cols[]		= {0,1,1};
	const bool type8_fraction_cols[]	= {0,0,0,1,0,0};
	const bool type8_strand_cols[]	= {0,0,0,0,0,1};
//X	762	+	CG	17	76
	const bool type9_strand_cols[]	= {0,0,1,0,0,0};
	const bool type9_cg_cols[]			= {0,0,0,1,0,0};
//chr1 10496 10497 79.69 64 + 10496 10497 180,60,0 0 0
//	const bool type10_chr_cols[]		= {1,0,0,0,0,0,0,0,0,0,0};
	const bool type10_int_cols[]		= {0,1,1};//,0,1,0,1,1,0,1,1};
	const bool type10_numeric_cols[]	= {0,1,1,1,1,0,1,1,0,1,1};
	const bool type10_strand_cols[]	= {0,0,0,0,0,1,0,0,0,0,0};
// chr1	C	3001631	CG	CG	1.0	5	5
//type 11 chr is the same as type5
	const bool type11_int_cols[]		= {0,0,1,0,0,0,1,1};
	const bool type11_float_cols[]	= {0,0,0,0,0,1,0,0};
/*//type12
1 6 - 0 2 CHH CAA (7 columns)
	const bool type12_strand_cols[]	= {0,0,1,0,0,0,0};
	const bool type12_cxx_cols[]		= {0,0,0,0,0,1,1};*/
	if (/* Input Type 1
Example:	chr1	762	763	0.1764	37
Column1:	chromosome
Column2:	nucleotide, an unsigned integer in [0,4294967295]
Column3:	ignored.
Column4:	methylation percent, a floating point in [0,1]
Column5:	coverage, an unsigned integer, an unsigned integer in [0,4294967295]
*/
	(no_of_cols == 5) &&
//	(bool_arrays_match(cols_are_chr, type1_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true) &&
	(
		(bool_arrays_match(cols_are_float, type1_float_cols, 1, no_of_cols) == true)
		||
		(bool_arrays_match(cols_are_float, all_cols_false, 1, no_of_cols) == true)
	) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(
		(bool_arrays_match(cols_are_int, type1_int_colsA, 1, no_of_cols) == true)
		||
		(bool_arrays_match(cols_are_int, type1_int_colsB, 1, no_of_cols) == true)
	) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 1, no_of_cols) == true)
	) {
		return_struct.type = 1;
	} else if (/*Input Type 2	known for MethylKit input
Example:	chr1.762	chr1	762	R	100	17.64	82.36
Column1:	unique name, this is ignored.
Column2:	chromosome
Column3:	nucleotide, an unsigned integer [0,4294967295]
Column4:	sense, this is ignored.
Column5:	coverage, an unsigned integer in [0,4294967295]
Column6:	methylation percent, a floating point in [0,100].
Column7:	cytosine percent, a floating point in [0,100].
*/
	(no_of_cols == 7) &&
//	(bool_arrays_match(cols_are_chr, type2_chr_cols, no_of_cols) == true) &&
//	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_float, type2_float_cols, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type2_int_cols, 2, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, type2_strand_cols, 2, no_of_cols) == true)
	) {
		return_struct.type = 2;
	} else if (/*Input Type 3
Example:	chr1	762	763	0.1764
Column1:	chromosome
Column2:	nucleotide, an unsigned integer [0,4294967295]
Column3:	nucleotide+1 (ignored)
Column4:	methylation percent, a floating point in [0,1].*/
	(no_of_cols == 4) &&
	((intcol2 + 1) == intcol3) &&
//	(bool_arrays_match(cols_are_chr, type3_chr_cols, no_of_cols) == true) && //type3 & 4 have same chr columns, no need to make new array identical to the first
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true) &&
	(cols_are_int[1] == true) &&
	(cols_are_int[2] == true) &&
	(cols_are_numeric[3] == true) && //can be float or integer
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 1, no_of_cols) == true)
	) {
		return_struct.type = 3;
	} else if (/*Input Type 4
Example:	chr1	762	6	14
Column1:	chromosome
Column2:	nucleotide, an unsigned integer in [0,4294967295]
Column3:	methylated C count, an unsigned integer in [0,4294967295]
Column4:	count unmethylated, an unsigned integer [0,4294967295]
*/
	(no_of_cols == 4) &&
//	(bool_arrays_match(cols_are_chr, type3_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_float, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type4_int_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 1, no_of_cols) == true)
) {
		return_struct.type = 4;
	} else if (/*Input Type 5
Example:	chr1	762	763	+	17	64	CG	CGA
Column1:	chromosome
Column2:	nucleotide/start position, an unsigned integer [0,4294967295]
Column3:	nucleotide/end position, ignored.
Column4:	strand
Column5:	methylated C count, an unsigned integer in [0,4294967295]
Column6:	count unmethylated, an unsigned integer in [0,4294967295]
Column7:	C-context, e.g. CG, CH, CHH.
Column8:	C-context, e.g. CGA, CGT, etc.*/
	(no_of_cols == 7) &&
//	(bool_arrays_match(cols_are_chr, type5_chr_cols, no_of_cols) == true) &&
//	(bool_arrays_match(cols_are_cg_chg_chh, type5_cg_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, type5_cg_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, type5_cxx_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_float, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type5_int_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, type5_strand_cols, 1, no_of_cols) == true)
) {
		return_struct.type = 5;
	} else if (
/*Input Type 6	Bismark coverage2cytosine format:
<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
Example:	chr1	762	763	0.223684	17	76
Column1:	chromosome, which is a string
Column2:	nucleotide/start position, an unsigned integer in [0,4294967295]
Column3:	nucleotide/end position, an unsigned integer in [0,4294967295]
Column4:	methylation percentage in [0,1]
Column5:	methylated C count, an unsigned integer in [0,4294967295]
Column6:	count unmethylated, an unsigned integer in [0,4294967295]*/
	(no_of_cols == 6) &&
//	(bool_arrays_match(cols_are_chr, type6_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type6_int_cols, 1, 3) == true) &&
	(bool_arrays_match(cols_are_numeric, type6_numeric_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 1, no_of_cols) == true)
	) {
		return_struct.type = 6;
	} else if (
/*Input Type 7	HELP-Tag data.  This can have a header.
Example:	20 chr1 43924 0 0.5 0.5	e.g. rn5.Control-1_Hpa.AC56EJACXX.lane_8.hcount
Column1:	ignored
Column2:	chromosome, a string
Column3:	position, an unsigned integer in [0,4294967295]
Column4:	methylation percent: a floating point number in [0,1]
Column5:	Conf. ignored.*/
	(no_of_cols == 5) &&
//	(bool_arrays_match(cols_are_chr, type7_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type7_int_cols, 2, 4) == true) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 2, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_numeric, type7_numeric_cols, 2, no_of_cols) == true)
	) {
		return_struct.type = 7;
	} else if (
/*Input Type 8	(EPP)Epigenome Processing Pipeline
Example:	chr1	762	763	'17/76'	999	+
Column1:	chromosome
Column2:	start nucleotide, an unsigned integer in [0,4294967295]
Column3:	end nucleotide, an unsigned integer in [0,4294967295]
Column4:	methylation percent as a fraction, two unsigned integers.  Coverage is given as the denominator. Everything after column 4 is ignored.*/
	(no_of_cols == 6) &&
//	(bool_arrays_match(cols_are_chr, type6_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_float, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, type8_fraction_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type8_int_cols, 1, 3) == true) &&
	(bool_arrays_match(col_strand_definit, type8_strand_cols, 1, no_of_cols) == true)
	) {
		return_struct.type = 8;
	} else if (
/*Input Type 9	Bsmooth Input
Example:	X	762	+	CG	17	76
Column1:	chromosome
Column2:	Nucleotide, an unsigned integer in [0,4294967295]
Column3:	strand sense, ignored
Column4:	context, ignored
Column5:	methylated C count, an unsigned integer in [0,4294967295]
Column6:	non-methylated C count, an unsigned integer in [0,4294967295]*/
	(no_of_cols == 6) &&
	(cols_are_int[1] == true) &&
	(cols_are_int[4] == true) &&
	(cols_are_int[5] == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, type9_cg_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, type9_strand_cols, 1, no_of_cols) == true)
	) {
		return_struct.type = 9;
	} else if (
/*Input Type 10	BisSNP (found in RnBeads)
Example:	chr1 10496 10497 79.69 64 + 10496 10497 180,60,0 0 0
Column1:	chromosome
Column2:	Nucleotide start, an unsigned integer in [0,4294967295]
Column3:	Nucleotide end, an unsigned integer in [0,4294967295]
Column4:	methylation value in [0:100].
Column5:	Coverage, an unsigned integer in [0,4294967295]. Everything after this column is ignored.*/
	(no_of_cols == 11) &&
//	(bool_arrays_match(cols_are_chr, type10_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_numeric, type10_numeric_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type10_int_cols, 1, 3) == true) &&
	(bool_arrays_match(col_strand_definit, type10_strand_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) == true)
	) {
		return_struct.type = 10;
	} else if (
/*Input Type 11	 BSSeeker https://github.com/BSSeeker/BSseeker2
Example:	chr1	C	3001631	CG	CG	1.0	5	5
Column1:	chromosome
Column2:	nucleotide on Watson (+) strand (ignored)
Column3:	Nucleotide end, an unsigned integer in [0,4294967295]
Column4:	context (CG/CHG/CHH) (ignored)
Column5:	dinucleotide-context (CA/CC/CG/CT) (ignored)
Column6:	methylation-level = #_of_C / (#_of_C + #_of_T)  (ignored)
Column7:	#_of_C (methylated C, the count of reads showing C here), an unsigned integer in [0,4294967295]
Column8:	#_of_C + #_of_T (all Cytosines, the count of reads showing C or T here), an unsigned integer in [0,4294967295]*/
	(no_of_cols == 8) &&
//	(bool_arrays_match(cols_are_chr, type5_chr_cols, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_float, type11_float_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, 1, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_int, type11_int_cols, 1, no_of_cols) == true) &&
	(bool_arrays_match(col_strand_definit, all_cols_false, 1, no_of_cols) == true)
	) {
		return_struct.type = 11;
/*	} else if (//<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
Input Type 12:	1 6 - 0 2 CHH CAA
Column1:	chromosome
Column2:	position
Column3:	strand +-
Column4:	count methylated
Column5:	count unmethylated
Column6: C context
Column7: trinucleotide context

	(no_of_cols == 7) &&
	(bool_arrays_match(col_strand_definit, type12_strand_cols, no_of_cols) == true) &&
	(cols_are_int[1] == true) &&
	(cols_are_int[3] == true) &&
	(cols_are_int[4] == true) &&
	(bool_arrays_match(cols_are_fraction, all_cols_false, no_of_cols) == true) &&
	(bool_arrays_match(cols_are_cg_chg_chh, type12_cxx_cols, no_of_cols) == true)
	) {
		return_struct.type = 12;*/
	} else {//everything below this prints details so I can debug
		free(line_copy); line_copy = NULL;
		fprintf(stderr, "could not identify file type for\n%s, which has %u columns.\n", line, no_of_cols);
		fprintf(stderr, "failed at %s line %u\n", __FILE__, __LINE__);
		for (unsigned short i = 0; i < no_of_cols; i++) {
			fprintf(stderr, "column %u:\n", i+1);
			fprintf(stderr, "\tcols_are_chr[%u] = %s\n", i, cols_are_chr[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_cg_chg_chh[%u] = %s\n", i, cols_are_cg_chg_chh[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_cx[%u] = %s\n", i, cols_are_cx[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_cxx[%u] = %s\n", i, cols_are_cxx[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_float[%u] = %s\n", i, cols_are_float[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_fraction[%u] = %s\n", i, cols_are_fraction[i] ? "true" : "false");
			fprintf(stderr, "\tcols_are_int[%u] = %s\n", i, cols_are_int[i] ? "true" : "false");
			fprintf(stderr, "\tcol_strand_definit[%u] = %s\n", i, col_strand_definit[i] ? "true" : "false");
		}
/*		puts("input type 2:");
		printf("cols_are_cx %s\n", bool_arrays_match(cols_are_cx, all_cols_false, 2, no_of_cols) ? "true" : "false");
		printf("cols_are_cxx %s\n", bool_arrays_match(cols_are_cxx, all_cols_false, 2, no_of_cols)  ? "true" : "false");
		printf("cols_are_float %s\n", bool_arrays_match(cols_are_float, type2_float_cols, 2, no_of_cols) ? "true" : "false");
		printf("cols_are_fraction %s\n", bool_arrays_match(cols_are_fraction, all_cols_false, 2, no_of_cols) ? "true" : "false");
		printf("cols_are_int %s\n", bool_arrays_match(cols_are_int, type2_int_cols, 2, no_of_cols) ? "true" : "false");
		printf("col_strand_definit %s\n", bool_arrays_match(col_strand_definit, type2_strand_cols, 2, no_of_cols) ? "true" : "false");*/
/*		puts("input type 5:");
		printf("cols_are_chr: %s\n", bool_arrays_match(cols_are_chr, type5_chr_cols, no_of_cols) ? "true" : "false");
//	(bool_arrays_match(cols_are_cg_chg_chh, type5_cg_cols, no_of_cols) == true) &&
		printf("cols_are_cx: %s\n", bool_arrays_match(cols_are_cx, type5_cg_cols, no_of_cols) ? "true" : "false");
		printf("cols_are_cxx: %s\n", bool_arrays_match(cols_are_cxx, type5_cxx_cols, no_of_cols) ? "true" : "false");
		printf("cols_are_float: %s\n", bool_arrays_match(cols_are_float, all_cols_false, no_of_cols) ? "true" : "false");
		printf("cols_are_fraction: %s\n", bool_arrays_match(cols_are_fraction, all_cols_false, no_of_cols) ? "true" : "false");
		printf("cols_are_int: %s\n", bool_arrays_match(cols_are_int, type5_int_cols, no_of_cols) ? "true" : "false");
		printf("col_strand_definit: %s\n", bool_arrays_match(col_strand_definit, type5_strand_cols, no_of_cols)  ? "true" : "false");*/
		puts("input type 8:");
		printf("cols_are_cg_chg_chh: %s\n", bool_arrays_match(cols_are_cg_chg_chh, all_cols_false, 1, no_of_cols) ? "true" : "false");
		printf("cols_are_cx: %s\n", bool_arrays_match(cols_are_cx, all_cols_false, 1, no_of_cols) ? "true" : "false");
		printf("cols_are_cxx: %s\n", bool_arrays_match(cols_are_cxx, all_cols_false, 1, no_of_cols) ? "true" : "false");
		printf("cols_are_float: %s\n", bool_arrays_match(cols_are_float, all_cols_false, 1, no_of_cols)  ? "true" : "false");
		printf("cols_are_fraction: %s\n", bool_arrays_match(cols_are_fraction, type8_fraction_cols, 1, no_of_cols)  ? "true" : "false");
		printf("cols_are_int: %s\n", bool_arrays_match(cols_are_int, type8_int_cols, 1, 3)  ? "true" : "false");
		printf("col_strand_definit: %s\n", bool_arrays_match(col_strand_definit, type8_strand_cols, 1, no_of_cols) ? "true" : "false");
		exit(EXIT_FAILURE);
	}
	free(line_copy); line_copy = NULL;
	return return_struct;
}

LINE_STRUCT read_data_line(const char *restrict line, const unsigned short int TYPE) {//, const bool FDEBUG) {//read line of data from file, necessary if the program reads multiple file types
/*	if (FDEBUG == true) {
		FILE_TYPE_STRUCT file_type_check = get_file_type(line);
		if (file_type_check.type != TYPE) {
			printf("line %s does not match file type @ %s line %u\n", line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}*/
	LINE_STRUCT return_struct;
	char *restrict bs_line_copy = strdup(line);
	if (bs_line_copy == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	char *saveptr, *restrict tmp_string;
	return_struct.nucleotide = 0;
	return_struct.methylated_C = 0;
	return_struct.strand = 0;
	return_struct.coverage = 0;
	return_struct.chromosome = NULL;
	chomp(bs_line_copy);//get rid of "\n" or "\r\n"
//	return_struct.sense = 0;
	if (TYPE == 1) {//bedgraph format
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		//safe_string_copy(return_struct.chromosome, tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);//get 2nd column as string
		tmp_string = strtok_r(NULL,"\t ",&saveptr);//throw away 3rd column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL,"\t ",&saveptr);//get 4th column as string
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		const double PERCENT = safe_strtod(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
 		return_struct.coverage = (unsigned)safe_strtoul(tmp_string);//5th column is coverage
 		return_struct.methylated_C = (unsigned)round(PERCENT*(double)return_struct.coverage);
	} else if (TYPE == 2) {//chrBase	chr	base	strand	coverage	freqC	freqT			methylSig input
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//like chr21
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//ignore first string of characters
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//like 43008720 = nucleotide
		if (tmp_string == NULL) {
			printf("I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);//export nucleotide to output struct
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//ignore sense
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.strand = get_strand(tmp_string);//F or R
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//will be coverage
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 6th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)safe_strtoul(tmp_string);//export nucleotide to output struct
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//will be percent
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 7th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)round(((safe_strtod(tmp_string)) / 100.0) * return_struct.coverage);
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
	} else if (TYPE == 3) {//chr10	119965	119966	0.03846
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//2nd column nucleotide
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);//get 2nd column as string
		tmp_string = strtok_r(NULL,"\t ",&saveptr);//discard 3rd column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL,"\t ",&saveptr);//get 4th column as string
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		const double PERCENT = (double)safe_strtod(tmp_string);//4th column is percent
		return_struct.coverage = 100000;//USHRT_MAX is higher coverage limit than I think a user could want
		return_struct.methylated_C = (unsigned)(round(PERCENT * return_struct.coverage));
	} else if (TYPE == 4) {//chr1	762	6	14
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//like chr21
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//like 43008720 = nucleotide
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);//export nucleotide to output struct
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//like 43008720 = nucleotide
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtoul(tmp_string);//export nucleotide to output struct
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//like 43008720 = nucleotide
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)(return_struct.methylated_C + safe_strtoul(tmp_string));//export nucleotide to output struct;
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
//		exit(EXIT_SUCCESS);
//		return_struct.methylated_C = (unsigned)(double)mC/return_struct.coverage;
//		printf("%u/%u = %f\n\n",return_struct.methylated_C,return_struct.coverage,return_struct.methylated_C);
//printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
	} else if (TYPE == 5) {//chr1	762	+	17	64	CG	CGA from Bismark
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column, chromosome "chr1"
		if (tmp_string == NULL) {
			fprintf(stderr, "Couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, position
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);//save 2nd column
		tmp_string = strtok_r(NULL,"\t ",&saveptr);//3rd column (strand +/-)
		if (tmp_string == NULL) {
			fprintf(stderr, "Couldn't get 3nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.strand = get_strand(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column, mC count
		if (tmp_string == NULL) {
			fprintf(stderr, "Couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th column, C count
		if (tmp_string == NULL) {
			fprintf(stderr, "Couldn't get 6th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = return_struct.methylated_C + (unsigned)safe_strtoul(tmp_string);
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);//DEBUGGING
	} else if (TYPE == 6) {//chr1	762	763	0.265625	17	76
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column is chromosome
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, position
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned int)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//ignore the 2nd position (3rd column)
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column, percent
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtod(tmp_string);//4th column is percent
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th column, methylated C count
		if (tmp_string == NULL) {
			printf("I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned int)safe_strtoul(tmp_string);
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 6th column, non-methylated C count
		if (tmp_string == NULL) {
			printf("I couldn't get 6th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)(return_struct.methylated_C + safe_strtoul(tmp_string));
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
	} else if (TYPE == 7) {//return 7::tid chr pos m-score conf (HELP-Tag)
		return_struct.coverage = 100000;
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column, throw it away
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, chromosome
		if (tmp_string == NULL) {
			printf("I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 3rd column, position + 1, & discard
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column, m-score (methylation percent)
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
//		const double PERCENT = strtod(tmp_string,NULL) / 100.0;
//		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th column, confidence
//		return_struct.methylated_C = return_struct.coverage * PERCENT;
		return_struct.methylated_C = (unsigned)round(safe_strtod(tmp_string) * return_struct.coverage / 100.0);
		if (return_struct.methylated_C > return_struct.coverage) {
			fprintf(stderr, "%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
			exit(EXIT_FAILURE);
		}
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %zu	coverage = %zu\n\n", line, return_struct.chromosome, return_struct.nucleotide,  return_struct.methylated_C, return_struct.coverage);//DEBUGGING
//the 5th column isn't relevant so I'm not going to waste time with strtok_r again
	} else if (TYPE == 8) {//"chr1	3010874	3010876	'10/14'	714	+"
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, position
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, position 2
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"'");
		char *restrict str2 = NULL, *subtoken;
		str2 = tmp_string;
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		const size_t FRACTION_LENGTH = strlen(str2);
		if ((FRACTION_LENGTH > 0) && (str2[FRACTION_LENGTH-1] == '\'')) {//only remove a closing quote, never a digit
			str2[FRACTION_LENGTH-1] = '\0';
		}
		str2 = strtok_r(str2, "/", &subtoken);
		if ((str2 == NULL) || (subtoken == NULL) || (subtoken[0] == '\0')) {
			fprintf(stderr, "I couldn't get a fraction like '10/14' from the 4th column of the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtoul(str2);
		return_struct.coverage = (unsigned)safe_strtoul(subtoken);
//		return_struct.methylated_C = (double)METHYLATED_COUNT/return_struct.coverage;
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);
	} else if (TYPE == 9) {//10	60025	+	CG	25	26
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, positioSn
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 3rd column, sense
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.strand = get_strand(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column, CpN type (not using it)
		if (tmp_string == NULL) {
			printf("I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th colum, number of methylated nucleotides
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtoul(tmp_string);//convert string to unsigned integer
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 6th column, coverage
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 6th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)safe_strtoul(tmp_string);
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage);//debugging
	} else if (TYPE == 10) {//return 10:	chr1	10496	10497	79.69	64	+	10496	10497	180,60,0 (BisSNP), found in RnBeads
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column, C of CpG
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 3rd column, nucleotide +1 (discard)
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column = methylation percent
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 4th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		const double PERCENT = safe_strtod(tmp_string) / 100.0;//make 4th column into double data type
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th column = coverage
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 5th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)safe_strtoul(tmp_string);
		return_struct.methylated_C = (unsigned)(round(PERCENT * return_struct.coverage));
//		printf("%schromosome = %s	nucleotide = %u	methylated_C = %u	coverage = %u, percent = %lf\n\n", line, return_struct.chromosome, return_struct.nucleotide, return_struct.methylated_C, return_struct.coverage, PERCENT);//debugging
	} else if (TYPE == 11) {//https://github.com/BSSeeker/BSseeker2
//chr1	C	3001631	CG	CG	1.0	5	5
		tmp_string = strtok_r(bs_line_copy, "\t ",&saveptr);//1st column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 1st column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		removeSubstring(tmp_string,"chr");
		return_struct.chromosome = strdup(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 2nd column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 2nd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 3rd column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 3rd column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.nucleotide = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 4th column
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 5th column
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 6th column
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 7th column
		if (tmp_string == NULL) {
			fprintf(stderr, "I couldn't get 7th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.methylated_C = (unsigned)safe_strtoul(tmp_string);
		tmp_string = strtok_r(NULL, "\t ",&saveptr);//get 8th column
		if (tmp_string == NULL) {
			printf("I couldn't get 8th column for the line \"%s\"	cf. %s line %u\n",line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		return_struct.coverage = (unsigned)safe_strtoul(tmp_string);
	} else {
		fprintf(stderr, "TYPE %u is not defined @ %s line %u\n\n", TYPE, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	free(bs_line_copy); bs_line_copy = NULL;
	if (return_struct.methylated_C > return_struct.coverage) {
		printf("%s has methylated_C = %u > coverage = %u, which isn't possible.\n", line, return_struct.methylated_C, return_struct.coverage);
		printf("Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return return_struct;
}
#endif
