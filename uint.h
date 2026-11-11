#ifndef DEFIANT_UINT_H
#define DEFIANT_UINT_H
#include <stdlib.h>
#include <stdio.h>

void * add_uint (unsigned int *restrict array, const unsigned int NUMBER, const unsigned int INDEX) {
	array = realloc(array, sizeof(unsigned int)*((size_t)1+INDEX));
	if (array == NULL) {
		fprintf(stderr, "array failed to realloc @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	array[INDEX] = NUMBER;
	return array;
}

void * uint_clear (unsigned int *restrict array) {
	free(array);//realloc(array, 0) is implementation-defined, and undefined as of C23
	return NULL;
}

void print_CpN_list(FILE *restrict fh, const unsigned int *restrict array, const unsigned int INDEX) {
	if ((array == NULL) || (INDEX == 0)) {//e.g. "-d 0": a DMR without any differential CpN
		fprintf(fh, "\tn/a");
		return;
	}
	fprintf(fh, "\t%u", array[0]);
	for (unsigned int i = 1; i < INDEX; i++) {
		fprintf(fh, ",%u", array[i]);
	}
}
#endif
