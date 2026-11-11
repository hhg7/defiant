#ifndef DEFIANT_OPTION_RANGE_H
#define DEFIANT_OPTION_RANGE_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
//requires safe_str_num.h (safe_strtod) to be included first

typedef struct {
	double min, max, step;
} option_range_struct;

typedef struct {
	unsigned int coverage, d, skip, CpN;
	unsigned int number_of_DMRs;
	double p, percent;
} cutoff_struct;

option_range_struct get_option_range (const char *restrict STRING) {
//accepts "value", "min,max" or "min,max,step"
	option_range_struct return_struct;
	if (strstr(STRING,",") == NULL) {
		return_struct.step = 1.0;//never used, as min = max
		return_struct.min = safe_strtod(STRING);
		return_struct.max = return_struct.min;
		return return_struct;
	}
	char *restrict copy_string = strdup(STRING);
	if (copy_string == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	char *restrict string, *string_pointer = NULL;
	string = strtok_r(copy_string, ",", &string_pointer);
	char *restrict max_string = strtok_r(NULL, ",", &string_pointer);
	if ((string == NULL) || (max_string == NULL)) {//e.g. "5," or ",5"
		fprintf(stderr, "'%s' isn't a range. A range looks like \"min,max\" or \"min,max,step\".\n", STRING);
		exit(EXIT_FAILURE);
	}
	return_struct.min = safe_strtod(string);
	return_struct.max = safe_strtod(max_string);
	if (return_struct.max < return_struct.min) {
		const double TMP = return_struct.max;
		return_struct.max = return_struct.min;
		return_struct.min = TMP;
	}
	string = strtok_r(NULL, ",", &string_pointer);
	if (string == NULL) {//i.e., the string ends here, so there are only two values
		return_struct.step = return_struct.max - return_struct.min;
	} else {
		return_struct.step = safe_strtod(string);
	}
	free(copy_string); copy_string = NULL;
	if ((return_struct.min != return_struct.max) && (return_struct.step <= 0.0)) {//a step <= 0 would loop forever
		fprintf(stderr, "The step in the range '%s' must be > 0.\n", STRING);
		exit(EXIT_FAILURE);
	}
	return return_struct;
}

void check_uint_option_range (const option_range_struct RANGE, const char *restrict OPTION) {
//options like -c, -CpN, -d and -S are unsigned integers, so they need integer steps >= 1, and can't be negative
	if ((RANGE.min < 0.0) || (RANGE.max > (double)UINT_MAX)) {
		fprintf(stderr, "'%s' must be between 0 and %u.\n", OPTION, UINT_MAX);
		exit(EXIT_FAILURE);
	}
	if ((RANGE.min != floor(RANGE.min)) || (RANGE.max != floor(RANGE.max))) {
		fprintf(stderr, "'%s' only accepts integers.\n", OPTION);
		exit(EXIT_FAILURE);
	}
	if ((RANGE.min != RANGE.max) && ((RANGE.step < 1.0) || (RANGE.step != floor(RANGE.step)))) {
		fprintf(stderr, "The step for '%s' must be an integer >= 1.\n", OPTION);
		exit(EXIT_FAILURE);
	}
}

bool doubles_equal (const double A, const double B) {
//p-values & percents are stepped in floating point, so 0.01+0.01+0.01+0.01+0.01 isn't necessarily == 0.05
	const double SCALE = fmax(1.0, fmax(fabs(A), fabs(B)));
	return fabs(A - B) <= 1e-9 * SCALE;
}

unsigned int uint_range_length (const option_range_struct RANGE) {//how many values are in min, min+step, ..., <= max
	const unsigned int MIN = (unsigned)RANGE.min, MAX = (unsigned)RANGE.max;
	const unsigned int STEP = RANGE.step >= 1.0 ? (unsigned)RANGE.step : 1;//a step of 0 would never end
	return (MAX - MIN) / STEP + 1;
}

unsigned int uint_range_value (const option_range_struct RANGE, const unsigned int INDEX) {
	const unsigned int STEP = RANGE.step >= 1.0 ? (unsigned)RANGE.step : 1;
	return (unsigned)RANGE.min + INDEX * STEP;
}

unsigned int double_range_length (const option_range_struct RANGE) {//how many values are in min, min+step, ..., <= max
	if (RANGE.step <= 0.0) {
		return 1;
	}
	return (unsigned)floor((RANGE.max - RANGE.min) / RANGE.step + 1e-9) + 1;
}

double double_range_value (const double DEFAULT, const option_range_struct RANGE, const unsigned int INDEX) {
//calculated from the index, not by repeated addition, so that rounding errors don't accumulate
	const double VALUE = RANGE.min + (double)INDEX * RANGE.step;
	if (doubles_equal(VALUE, DEFAULT) == true) {
		return DEFAULT;//exactly the default, so that == comparisons with the default work
	}
	if (doubles_equal(VALUE, RANGE.max) == true) {
		return RANGE.max;
	}
	return VALUE;
}

unsigned int how_many_uint_array_elements (const unsigned int DEFAULT, const option_range_struct RANGE, const unsigned int CPU) {
	if (RANGE.min == RANGE.max) {
		if ((unsigned)RANGE.min == DEFAULT) {
			return 1;
		}
		return CPU > 1 ? 2 : 1;//will return the value and the default
	}
	const unsigned int LENGTH = uint_range_length(RANGE);
	bool defaults_present = false;
	for (unsigned int i = 0; i < LENGTH; i++) {
		if (uint_range_value(RANGE, i) == DEFAULT) {
			defaults_present = true;
		}
	}
	if ((defaults_present == false) && (CPU > 1)) {
		return LENGTH + 1;
	}
	return LENGTH;
}

unsigned int how_many_double_array_elements (const double DEFAULT, const option_range_struct RANGE, const unsigned int CPU) {
	if (RANGE.min == RANGE.max) {
		if (RANGE.min == DEFAULT) {
			return 1;
		}
		return CPU > 1 ? 2 : 1;//will return the value and the default
	}
	const unsigned int LENGTH = double_range_length(RANGE);
	bool defaults_present = false;
	for (unsigned int i = 0; i < LENGTH; i++) {
		if (double_range_value(DEFAULT, RANGE, i) == DEFAULT) {
			defaults_present = true;
		}
	}
	if ((defaults_present == false) && (CPU > 1)) {
		return LENGTH + 1;
	}
	return LENGTH;
}

void * uint_array_alloc (const unsigned int DEFAULT, const option_range_struct RANGE, const unsigned int CPU) {
	const unsigned int NUMBER_OF_ARRAY_ELEMENTS = how_many_uint_array_elements(DEFAULT, RANGE, CPU);
	unsigned int *restrict array = malloc(sizeof(unsigned int) * NUMBER_OF_ARRAY_ELEMENTS);
	if (array == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if (RANGE.min == RANGE.max) {
		if (NUMBER_OF_ARRAY_ELEMENTS == 2) {
			array[0] = DEFAULT;
			array[1] = (unsigned)RANGE.max;
		} else {
			array[0] = (unsigned)RANGE.max;
		}
		return array;
	}
	const unsigned int LENGTH = uint_range_length(RANGE);
	for (unsigned int i = 0; i < LENGTH; i++) {
		array[i] = uint_range_value(RANGE, i);
	}
	if (NUMBER_OF_ARRAY_ELEMENTS > LENGTH) {
		array[LENGTH] = DEFAULT;
	}
	return array;
}

void * double_array_alloc (const double DEFAULT, const option_range_struct RANGE, const unsigned int CPU) {
	const unsigned int NUMBER_OF_ARRAY_ELEMENTS = how_many_double_array_elements(DEFAULT, RANGE, CPU);
	double *restrict array = malloc(sizeof(double) * NUMBER_OF_ARRAY_ELEMENTS);
	if (array == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	if (RANGE.min == RANGE.max) {
		if (NUMBER_OF_ARRAY_ELEMENTS == 2) {
			array[0] = DEFAULT;
			array[1] = RANGE.max;
		} else {
			array[0] = RANGE.max;
		}
		return array;
	}
	const unsigned int LENGTH = double_range_length(RANGE);
	for (unsigned int i = 0; i < LENGTH; i++) {
		array[i] = double_range_value(DEFAULT, RANGE, i);
	}
	if (NUMBER_OF_ARRAY_ELEMENTS > LENGTH) {
		array[LENGTH] = DEFAULT;
	}
	return array;
}

unsigned int how_many_runs (const option_range_struct coverage, const option_range_struct CpN, const option_range_struct d, const option_range_struct p, const option_range_struct percent, const option_range_struct skip, const cutoff_struct defaults, const unsigned short int ALLOWED_NON_DEFAULTS, const unsigned int CPU) {
	const unsigned int COVERAGE_ITERATIONS = how_many_uint_array_elements(defaults.coverage, coverage, CPU);
	const unsigned int CPN_ITERATIONS = how_many_uint_array_elements(defaults.CpN, CpN, CPU);
	const unsigned int D_ITERATIONS = how_many_uint_array_elements(defaults.d, d, CPU);
	const unsigned int P_ITERATIONS = how_many_double_array_elements(defaults.p, p, CPU);
	const unsigned int PERCENT_ITERATIONS = how_many_double_array_elements(defaults.percent, percent, CPU);
	const unsigned int SKIP_ITERATIONS = how_many_uint_array_elements(defaults.skip, skip, CPU);
	unsigned int *restrict coverage_array = uint_array_alloc(defaults.coverage, coverage, CPU);
	unsigned int *restrict cpn_array = uint_array_alloc(defaults.CpN, CpN, CPU);
	unsigned int *restrict d_array = uint_array_alloc(defaults.d, d, CPU);
	double *restrict p_array = double_array_alloc(defaults.p, p, CPU);
	double *restrict percent_array = double_array_alloc(defaults.percent, percent, CPU);
	unsigned int *restrict skip_array = uint_array_alloc(defaults.skip, skip, CPU);
	unsigned int runs = 0;
	for (unsigned int c = 0; c < COVERAGE_ITERATIONS; c++) {
		for (unsigned int CpG = 0; CpG < CPN_ITERATIONS; CpG++) {
			for (unsigned int diff = 0; diff < D_ITERATIONS; diff++) {
				for (unsigned int pvalue = 0; pvalue < P_ITERATIONS; pvalue ++) {
					for (unsigned int Percent = 0; Percent < PERCENT_ITERATIONS; Percent ++) {
						for (unsigned int skips = 0; skips < SKIP_ITERATIONS; skips ++) {
							unsigned short int non_default_values = 0;
							if (coverage_array[c] != defaults.coverage) {
								non_default_values++;
							}
							if (cpn_array[CpG] != defaults.CpN) {
								non_default_values++;
							}
							if (d_array[diff] != defaults.d) {
								non_default_values++;
							}
							if (p_array[pvalue] != defaults.p) {
								non_default_values++;
							}
							if (percent_array[Percent] != defaults.percent) {
								non_default_values++;
							}
							if (skip_array[skips] != defaults.skip) {
								non_default_values++;
							}
//							printf("c%u	CpN%u	p%.2f	P%.2f	S%u\n", c, CpG, pvalue, Percent, skips);
							if (non_default_values <= ALLOWED_NON_DEFAULTS) {
								runs++;
							}
						}
					}
				}
			}
		}
	}
	free(coverage_array);	coverage_array = NULL;
	free(cpn_array);			cpn_array = NULL;
	free(d_array);				d_array = NULL;
	free(p_array);				p_array = NULL;
	free(percent_array);		percent_array = NULL;
	free(skip_array);			skip_array = NULL;
	return runs;
}

cutoff_struct * make_cutoff_struct (const option_range_struct coverage, const option_range_struct CpN, const option_range_struct d, const option_range_struct p, const option_range_struct percent, const option_range_struct skip,  const cutoff_struct defaults, const unsigned short int ALLOWED_NON_DEFAULTS, const unsigned int CPU) {
	const unsigned int RUNS = how_many_runs(coverage, CpN, d, p, percent, skip, defaults, ALLOWED_NON_DEFAULTS, CPU);
	if (RUNS == 0) {
		return NULL;
	}
	cutoff_struct *restrict array = malloc(RUNS * sizeof(cutoff_struct));
	const unsigned int COVERAGE_ITERATIONS = how_many_uint_array_elements(defaults.coverage, coverage, CPU);
	const unsigned int CPN_ITERATIONS = how_many_uint_array_elements(defaults.CpN, CpN, CPU);
	const unsigned int D_ITERATIONS = how_many_uint_array_elements(defaults.d, d, CPU);
	const unsigned int P_ITERATIONS = how_many_double_array_elements(defaults.p, p, CPU);
	const unsigned int PERCENT_ITERATIONS = how_many_double_array_elements(defaults.percent, percent, CPU);
	const unsigned int SKIP_ITERATIONS = how_many_uint_array_elements(defaults.skip, skip, CPU);
	unsigned int *restrict coverage_array = uint_array_alloc(defaults.coverage, coverage, CPU);
	unsigned int *restrict cpn_array = uint_array_alloc(defaults.CpN, CpN, CPU);
	unsigned int *restrict d_array = uint_array_alloc(defaults.d, d, CPU);
	double *restrict p_array = double_array_alloc(defaults.p, p, CPU);
	double *restrict percent_array = double_array_alloc(defaults.percent, percent, CPU);
	unsigned int *restrict skip_array = uint_array_alloc(defaults.skip, skip, CPU);
	unsigned int run = 0;//first run is all defaults
	for (unsigned int c = 0; c < COVERAGE_ITERATIONS; c++) {
		for (unsigned int CpG = 0; CpG < CPN_ITERATIONS; CpG++) {
	//		printf("	CpN = %u, %u/%u\n", cpn_array[CpG], CpG, CPN_ITERATIONS);
			for (unsigned int diff = 0; diff < D_ITERATIONS; diff++) {
//				printf("		d = %u, %u/%u\n", d_array[diff], diff, D_ITERATIONS);
				for (unsigned int pvalue = 0; pvalue < P_ITERATIONS; pvalue ++) {
//					printf("			p = %lf, %u/%u\n", p_array[pvalue], pvalue, P_ITERATIONS);
					for (unsigned int Percent = 0; Percent < PERCENT_ITERATIONS; Percent ++) {
//						printf("				%% = %lf,	 %u/%u\n", percent_array[Percent], Percent, PERCENT_ITERATIONS);
						for (unsigned int skips = 0; skips < SKIP_ITERATIONS; skips++) {
//							printf("					skips = %u, %u/%u\n", skip_array[skips], skips, SKIP_ITERATIONS);
							unsigned short int non_default_values = 0;
							if (coverage_array[c] != defaults.coverage) {
								non_default_values++;
							}
							if (cpn_array[CpG] != defaults.CpN) {
								non_default_values++;
							}
							if (d_array[diff] != defaults.d) {
								non_default_values++;
							}
							if (p_array[pvalue] != defaults.p) {
								non_default_values++;
							}
							if (percent_array[Percent] != defaults.percent) {
								non_default_values++;
							}
							if (skip_array[skips] != defaults.skip) {
								non_default_values++;
							}
							if (non_default_values <= ALLOWED_NON_DEFAULTS) {
								if (run >= RUNS) {
									fprintf(stderr, "You are about to segfault at @ %s line %u\n", __FILE__, __LINE__);
									exit(EXIT_FAILURE);
								}
								array[run].coverage			= coverage_array[c];
								array[run].CpN					= cpn_array[CpG];
								array[run].d					= d_array[diff];
								array[run].p					= p_array[pvalue];
								array[run].percent			= percent_array[Percent];
								array[run].skip				= skip_array[skips];
								array[run].number_of_DMRs	= 0;
//								printf("run %u	c = %u	CpN = %u	d = %u	p = %.2f	P = %lf	S = %u\n", run, array[run].coverage, array[run].CpN, array[run].d, array[run].p, array[run].percent, array[run].skip);
								run++;
							}
						}
					}
				}
			}
		}
	}
	free(coverage_array); coverage_array = NULL;
	free(cpn_array); cpn_array = NULL;
	free(d_array);	d_array = NULL;
	free(p_array); p_array = NULL;
	free(percent_array); percent_array = NULL;
	free(skip_array); skip_array = NULL;
	return array;
}
#endif
