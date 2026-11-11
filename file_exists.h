#ifndef DEFIANT_FILE_EXISTS_H
#define DEFIANT_FILE_EXISTS_H
#include <stdio.h>
#include <unistd.h>//access
#include <stdbool.h>//true and false definitions

bool file_exists (const char *restrict fname) {
	if (access( fname, F_OK ) != -1 ) {
		return true;// file exists
	} else {
		return false;// file doesn't exist
	}
}
#endif
