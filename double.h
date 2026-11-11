#ifndef DEFIANT_DOUBLE_H
#define DEFIANT_DOUBLE_H
#include <stdlib.h>
#include <stdio.h>

void * add_double (double *restrict array, const double NUMBER, const unsigned int INDEX) {
	array = realloc(array, sizeof(double)*((size_t)1+INDEX));
	if (array == NULL) {
		fprintf(stderr, "array failed to realloc @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	array[INDEX] = NUMBER;
	return array;
}

void * double_clear (double *restrict array) {
	free(array);//realloc(array, 0) is implementation-defined, and undefined as of C23
	return NULL;
}

void print_double_CpN_list(FILE *restrict fh, const double *restrict array, const unsigned int INDEX) {
	if ((array == NULL) || (INDEX == 0)) {
		fprintf(fh, "\tn/a");
		return;
	}
	fprintf(fh, "\t%g", array[0]);
	for (unsigned int i = 1; i < INDEX; i++) {
		fprintf(fh, ",%g", array[i]);
	}
}
#endif
