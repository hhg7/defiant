#ifndef DEFIANT_MY_STRING_H
#define DEFIANT_MY_STRING_H
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

unsigned short int count_occurences_of_string (const char *restrict BIG_STRING, const char *restrict SMALL_STRING) {//count (non-overlapping) occurences of a small string in a big string
	const size_t SMALL_LENGTH = strlen(SMALL_STRING);
	if (SMALL_LENGTH == 0) {
		return 0;
	}
	unsigned short int count = 0;
	const char *position = BIG_STRING;
	while ((position = strstr(position, SMALL_STRING)) != NULL) {
		count++;
		position += SMALL_LENGTH;
	}
	return count;
}

unsigned short int count_tokens (const char *restrict STRING, const char DELIMITER) {
//the number of non-empty fields, i.e. the number of tokens that strtok/strtok_r will return, so "a,,b," has 2 tokens
	unsigned short int count = 0;
	bool inside_token = false;
	for (size_t i = 0; STRING[i] != '\0'; i++) {
		if (STRING[i] == DELIMITER) {
			inside_token = false;
		} else if (inside_token == false) {
			inside_token = true;
			count++;
		}
	}
	return count;
}
#endif
