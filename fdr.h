#ifndef DEFIANT_FDR_H
#define DEFIANT_FDR_H
#include <stdio.h>//printf
#include <stdlib.h>//qsort
#include <math.h>//fabs
#include <stdbool.h>//bool data type
#include <strings.h>//strcasecmp
#include <assert.h>//assert, necessary for random integer selection

unsigned int * seq_len(const unsigned int START, const unsigned int END) {
//named after R function of same name, but simpler function
	unsigned start = (unsigned)START;
	unsigned end = (unsigned)END;
	if (START == END) {
		unsigned int *restrict sequence = malloc( (end+1) * sizeof(unsigned int));
		if (sequence == NULL) {
			fprintf(stderr, "malloc failed at %s line %u\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		for (unsigned i = 0; i < end; i++) {
			sequence[i] = i+1;
		}
		return sequence;
	}
	if (START > END) {
		end = (unsigned)START;
		start = (unsigned)END;
	}
	const unsigned LENGTH = end - start ;
	unsigned int *restrict sequence = malloc( (1+LENGTH) * sizeof(unsigned int));
	if (sequence == NULL) {
		fprintf(stderr, "malloc failed at %s line %u\n", __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	if (START < END) {
		for (unsigned index = 0; index <= LENGTH; index++) {
			sequence[index] = start + index;
		}
	} else {
		for (unsigned index = 0; index <= LENGTH; index++) {
			sequence[index] = end - index;
		}
	}
	return sequence;
}

//modified from https://phoxis.org/2012/07/12/get-sorted-index-orderting-of-an-array/

double *restrict base_arr = NULL;

static int compar_increase (const void *restrict a, const void *restrict b) {
	int aa = *((int *restrict ) a), bb = *((int *restrict) b);
	if (base_arr[aa] < base_arr[bb]) {
		return 1;
	} else if (base_arr[aa] == base_arr[bb]) {
		return 0;
	} else {
		return -1;
	}
}

static int compar_decrease (const void *restrict a, const void *restrict b) {
	int aa = *((int *restrict ) a), bb = *((int *restrict) b);
	if (base_arr[aa] < base_arr[bb]) {
		return -1;
	} else if (base_arr[aa] == base_arr[bb]) {
		return 0;
	} else {
		return 1;
	}
}

unsigned int * order (const double *restrict ARRAY, const unsigned int SIZE, const bool DECREASING) {
//this has the same name as the same R function
	unsigned int *restrict idx = malloc(SIZE * sizeof(unsigned int));
	if (idx == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	base_arr = malloc(sizeof(double) * SIZE);
	if (base_arr == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	for (unsigned int i = 0; i < SIZE; i++) {
		base_arr[i] = ARRAY[i];
		idx[i] = i;
	}
	if (DECREASING == false) {
		qsort(idx, SIZE, sizeof(unsigned int), compar_decrease);
	} else if (DECREASING == true) {
		qsort(idx, SIZE, sizeof(unsigned int), compar_increase);
	}
	free(base_arr); base_arr = NULL;
	return idx;
}

double * cummin(const double *restrict ARRAY, const unsigned int NO_OF_ARRAY_ELEMENTS) {
//this takes the same name of the R function which it copies
//this requires a free() afterward where it is used
	if (NO_OF_ARRAY_ELEMENTS < 1) {
		fprintf(stderr, "cummin function requires at least one element.\n");
		fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	double *restrict output = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
	if (output == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	double cumulative_min = ARRAY[0];
	for (unsigned int i = 0; i < NO_OF_ARRAY_ELEMENTS; i++) {
		if (ARRAY[i] < cumulative_min) {
			cumulative_min = ARRAY[i];
		}
		output[i] = cumulative_min;
	}
	return output;
}

double * cummax(const double *restrict ARRAY, const unsigned int NO_OF_ARRAY_ELEMENTS) {
//this takes the same name of the R function which it copies
//this requires a free() afterward where it is used
	if (NO_OF_ARRAY_ELEMENTS < 1) {
		fprintf(stderr, "function requires at least one element.\n");
		fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	double *restrict output = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
	if (output == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	double cumulative_max = ARRAY[0];
	for (unsigned int i = 0; i < NO_OF_ARRAY_ELEMENTS; i++) {
		if (ARRAY[i] > cumulative_max) {
			cumulative_max = ARRAY[i];
		}
		output[i] = cumulative_max;
	}
	return output;
}

double * pminx(const double *restrict ARRAY, const unsigned int NO_OF_ARRAY_ELEMENTS, const double X) {
//named after the R function pmin
	if (NO_OF_ARRAY_ELEMENTS < 1) {
		fprintf(stderr, "pmin requires at least one element.\n");
		fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	double *restrict pmin_array = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
	if (pmin_array == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
		if (ARRAY[index] < X) {
			pmin_array[index] = ARRAY[index];
		} else {
			pmin_array[index] = X;
		}
	}
	return pmin_array;
}

void double_say (const double *restrict ARRAY, const size_t NO_OF_ARRAY_ELEMENTS) {
	printf("[1] %e", ARRAY[0]);
	for (unsigned int i = 1; i < NO_OF_ARRAY_ELEMENTS; i++) {
		printf(" %.10f", ARRAY[i]);
		if (((i+1) % 5) == 0) {
			printf("\n[%u]", i+1);
		}
	}
	puts("\n");
}

/*void uint_say (const unsigned int *restrict ARRAY, const size_t NO_OF_ARRAY_ELEMENTS) {
//for debugging
	printf("%u", ARRAY[0]);
	for (size_t i = 1; i < NO_OF_ARRAY_ELEMENTS; i++) {
		printf(",%u", ARRAY[i]);
	}
	puts("\n");
}*/

double * uint2double (const unsigned int *restrict ARRAY, const unsigned int NO_OF_ARRAY_ELEMENTS) {
	double *restrict doubleArray = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
	if (doubleArray == NULL) {
		fprintf(stderr, "Failure to malloc at %s line %u.\n", __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
		doubleArray[index] = (double)ARRAY[index];
	}
	return doubleArray;
}

double min2 (const double N1, const double N2) {
	if (N1 < N2) {
		return N1;
	} else {
		return N2;
	}
}

double * p_adjust (const double *restrict PVALUES, const unsigned int NO_OF_ARRAY_ELEMENTS, const char *restrict STRING) {
//this function is a translation of R's p.adjust "BH" method
// i is always i[index] = NO_OF_ARRAY_ELEMENTS - index - 1
	if (NO_OF_ARRAY_ELEMENTS < 1) {
		fprintf(stderr, "You may wish to run without the -fdr option\n");
		fprintf(stderr, "p_adjust requires at least one element.\n");
		fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	short int TYPE = -1;
	if (STRING == NULL) {
		TYPE = 0;
	} else if (strcasecmp(STRING, "BH") == 0) {
		TYPE = 0;
	} else if (strcasecmp(STRING, "fdr") == 0) {
		TYPE = 0;
	} else if (strcasecmp(STRING, "by") == 0) {
		TYPE = 1;
	} else if (strcasecmp(STRING, "Bonferroni") == 0) {
		TYPE = 2;
	} else if (strcasecmp(STRING, "hochberg") == 0) {
		TYPE = 3;
	} else if (strcasecmp(STRING, "holm") == 0) {
		TYPE = 4;
	} else if (strcasecmp(STRING, "hommel") == 0) {
		TYPE = 5;
	} else {
		fprintf(stderr, "%s doesn't match any accepted FDR methods.\n", STRING);
		fprintf(stderr, "Failed at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
	if (TYPE == 2) {//Bonferroni method
		double *restrict bonferroni = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		if (bonferroni == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			const double BONFERRONI = PVALUES[index] * NO_OF_ARRAY_ELEMENTS;
			if (BONFERRONI >= 1.0) {
				bonferroni[index] = 1.0;
			} else if ((0.0 <= BONFERRONI) && (BONFERRONI < 1.0)) {
				bonferroni[index] = BONFERRONI;
			} else {
				fprintf(stderr, "%g is outside of the interval I planned.\n", BONFERRONI);
				fprintf(stderr, "Failure at %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		return bonferroni;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
	} else if (TYPE == 4) {//Holm method
/*these values are computed separately from BH, BY, and Hochberg because they are
computed differently*/
		unsigned int *restrict o  = order(PVALUES, NO_OF_ARRAY_ELEMENTS, false);
//sorted in reverse of methods 0-3
		double *restrict o2double = uint2double(o, NO_OF_ARRAY_ELEMENTS);
		double *restrict cummax_input = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		for (unsigned index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			cummax_input[index] = (NO_OF_ARRAY_ELEMENTS - index ) * (double)PVALUES[o[index]];
//			printf("cummax_input[%zu] = %e\n", index, cummax_input[index]);
		}
		free(o); o = NULL;
		unsigned int *restrict ro = order(o2double, NO_OF_ARRAY_ELEMENTS, false);
		free(o2double); o2double = NULL;

		double *restrict cummax_output = cummax(cummax_input, NO_OF_ARRAY_ELEMENTS);
		free(cummax_input); cummax_input = NULL;

		double *restrict pmin = pminx(cummax_output, NO_OF_ARRAY_ELEMENTS, 1);
		free(cummax_output); cummax_output = NULL;
		double *restrict qvalues = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			qvalues[index] = pmin[ro[index]];
		}
		free(pmin); pmin = NULL;
		free(ro); ro = NULL;
		return qvalues;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
	} else if (TYPE == 5) {//Hommel method
//i <- seq_len(n)
//o <- order(p)
		unsigned int *restrict o = order(PVALUES, NO_OF_ARRAY_ELEMENTS, false);//false is R's default
//p <- p[o]
		double *restrict p = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		if (p == NULL) {
			printf("failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		for (size_t index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			p[index] = PVALUES[o[index]];
		}
//ro <- order(o)
		double *restrict o2double = uint2double(o, NO_OF_ARRAY_ELEMENTS);
		free(o); o = NULL;
		unsigned int *restrict ro = order(o2double, NO_OF_ARRAY_ELEMENTS, false);
		free(o2double); o2double = NULL;
//		puts("ro");
//q <- pa <- rep.int(min(n * p/i), n)
		double *restrict q   = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		if (q == NULL) {
			printf("failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		double *restrict pa  = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
		if (pa == NULL) {
			printf("failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		double min = (double)NO_OF_ARRAY_ELEMENTS * p[0];
		for (size_t index = 1; index < NO_OF_ARRAY_ELEMENTS; index++) {
			const double TEMP = (double)NO_OF_ARRAY_ELEMENTS * p[index] / (double)(1+index);
			if (TEMP < min) {
				min = TEMP;
			}
		}
		for (size_t index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			pa[index] = min;
			 q[index] = min;
		}
//		puts("q & pa");
//\\		double_say(q, NO_OF_ARRAY_ELEMENTS);
/*for (j in (n - 1):2) {
            ij <- seq_len(n - j + 1)
            i2 <- (n - j + 2):n
            q1 <- min(j * p[i2]/(2:j))
            q[ij] <- pmin(j * p[ij], q1)
            q[i2] <- q[n - j + 1]
            pa <- pmax(pa, q)
        }
*/
		for (size_t j = (NO_OF_ARRAY_ELEMENTS-1); j >= 2; j--) {
//			printf("j = %zu\n", j);
			unsigned int *restrict ij = seq_len(0,(unsigned)(NO_OF_ARRAY_ELEMENTS - j));
			const size_t I2_LENGTH = j - 1;
			unsigned int *restrict i2 = malloc(I2_LENGTH * sizeof(unsigned int));
			for (size_t i = 0; i < I2_LENGTH; i++) {
				i2[i] = (unsigned)(NO_OF_ARRAY_ELEMENTS-j+2+i-1);
//R's indices are 1-based, C's are 0-based, I added the -1
			}

			double q1 = (double)j * p[i2[0]] / 2.0;
			for (size_t i = 1; i < I2_LENGTH; i++) {//loop through 2:j
				const double TEMP_Q1 = (double)j * p[i2[i]] / (double)(2 + i);
				if (TEMP_Q1 < q1) {
					q1 = TEMP_Q1;
				}
			}

			for (size_t i = 0; i < (NO_OF_ARRAY_ELEMENTS - j + 1); i++) {//q[ij] <- pmin(j * p[ij], q1)
				q[ij[i]] = min2( (double)j*p[ij[i]], q1);
			}
			free(ij); ij = NULL;

			for (size_t i = 0; i < I2_LENGTH; i++) {//q[i2] <- q[n - j + 1]
				q[i2[i]] = q[NO_OF_ARRAY_ELEMENTS - j];//subtract 1 because of starting index difference
			}
			free(i2); i2 = NULL;

			for (size_t i = 0; i < NO_OF_ARRAY_ELEMENTS; i++) {//pa <- pmax(pa, q)
				if (pa[i] < q[i]) {
					pa[i] = q[i];
				}
			}
//			printf("j = %zu, pa = \n", j);
//				double_say(pa, N);
		}//end j loop
		free(p); p = NULL;
		for (size_t index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			q[index] = pa[ro[index]];//Hommel q-values
		}
//now free memory
		free(ro); ro = NULL;
		free(pa); pa = NULL;
		return q;
	}
//The methods are similarly computed and thus can be combined for clarity
	unsigned int *restrict o = order(PVALUES, NO_OF_ARRAY_ELEMENTS, true);
	if (o == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
//	double_say(PVALUES, NO_OF_ARRAY_ELEMENTS);
	double *restrict o_double = uint2double(o, NO_OF_ARRAY_ELEMENTS);
	for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
		if ((PVALUES[index] < 0) || (PVALUES[index] > 1)) {
			fprintf(stderr, "array[%u] = %lf, which is outside the interval [0,1]\n", index, PVALUES[index]);
			fprintf(stderr, "died at %s line %u\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	unsigned int *restrict ro = order(o_double, NO_OF_ARRAY_ELEMENTS, false);
	if (ro == NULL) {
			fprintf(stderr, "failed to malloc at %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
	}
	free(o_double); o_double = NULL;
	double *restrict cummin_input = malloc(sizeof(double) * NO_OF_ARRAY_ELEMENTS);
	if (TYPE == 0) {//BH method
		for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			const double NI = (double)NO_OF_ARRAY_ELEMENTS / (double)(NO_OF_ARRAY_ELEMENTS - index);// n/i simplified 
			cummin_input[index] = NI * PVALUES[o[index]];//PVALUES[o[index]] is p[o]
		}
	} else if (TYPE == 1) {//BY method
		double q = 1.0;
		for (unsigned int index = 2; index < (1+NO_OF_ARRAY_ELEMENTS); index++) {
			q +=  1.0/(double)index;
		}
		for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
			const double NI = (double)NO_OF_ARRAY_ELEMENTS / (double)(NO_OF_ARRAY_ELEMENTS - index);// n/i simplified 
			cummin_input[index] = q * NI * PVALUES[o[index]];//PVALUES[o[index]] is p[o]
		}
	} else if (TYPE == 3) {//Hochberg method
		for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
// pmin(1, cummin((n - i + 1L) * p[o]))[ro]
			cummin_input[index] = (double)(index + 1) * PVALUES[o[index]];
		}
	}
	free(o); o = NULL;
	double *restrict cummin_array = NULL;
	cummin_array = cummin(cummin_input, NO_OF_ARRAY_ELEMENTS);
	free(cummin_input); cummin_input = NULL;//I don't need this anymore
	double *restrict pmin = pminx(cummin_array, NO_OF_ARRAY_ELEMENTS, 1);
	free(cummin_array); cummin_array = NULL;
	double *restrict q_array = malloc(NO_OF_ARRAY_ELEMENTS*sizeof(double));
	for (unsigned int index = 0; index < NO_OF_ARRAY_ELEMENTS; index++) {
		q_array[index] = pmin[ro[index]];
	}

	free(ro); ro = NULL;
	free(pmin); pmin = NULL;
	return q_array;
}
#endif
