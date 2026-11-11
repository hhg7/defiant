#ifndef DEFIANT_SAFE_STR_NUM_H
#define DEFIANT_SAFE_STR_NUM_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>//isfinite
#include <ctype.h>//isspace

double safe_strtod(const char *restrict STRING) {
	if ((STRING == NULL) || (STRING[0] == '\0')) {
		fprintf(stderr, "An empty string can't be converted to a double precision floating point @ %s line %u.\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	char *endptr;
	errno = 0;
	const double VALUE = strtod(STRING, &endptr);
	if ((endptr == STRING) || (*endptr != '\0') || (errno == ERANGE) || (!isfinite(VALUE))) {//"nan" & "inf" convert, but can't be used
		fprintf(stderr, "The string '%s' failed to convert to a double precision floating point @ %s line %u.\n", STRING, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return VALUE;
}//slightly modified from mousio on http://stackoverflow.com/questions/5601537/converting-a-string-to-double/5601559#5601559

unsigned long int safe_strtoul (const char *restrict STRING) {
	if ((STRING == NULL) || (STRING[0] == '\0')) {
		fprintf(stderr, "An empty string can't be converted to an unsigned integer @ %s line %u.\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	size_t first = 0;
	while (isspace((unsigned char)STRING[first])) {
		first++;
	}
	if (STRING[first] == '-') {//strtoul silently wraps negative numbers around to huge positive numbers
		fprintf(stderr, "Negative numbers like '%s' can't be fed to unsigned integers @ %s line %u.\n", STRING, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	char *endptr;
	errno = 0;
	const unsigned long int VALUE = strtoul(STRING, &endptr, 10);//base 10: base 0 would read "010" as octal 8
	if ((endptr == STRING) || (*endptr != '\0') || (errno == ERANGE)) {
		fprintf(stderr, "The string '%s' failed to convert to a unsigned integer @ %s line %u.\n", STRING, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return VALUE;
}
#endif
