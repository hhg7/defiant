#ifndef DEFIANT_ANNOTATION_H
#define DEFIANT_ANNOTATION_H
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
//requires safe_str_num.h (safe_strtoul) to be included first, and _GNU_SOURCE to be defined before any system header

void * add_name(char **restrict array, size_t *restrict size, unsigned int *restrict f_number_of_list_genes, const char *restrict NAME, const bool UNIQUE, const bool KEEP_NA) {//this function works, when freed inside main()
	if ((KEEP_NA == false) && (strcmp(NAME,"n/a") == 0)) {
		return array;
	}//don't bother with empty strings
//	printf("@line %u, adding %s to gene list array, with size %zu and %u elements.\n", __LINE__, NAME, *size, *f_number_of_list_genes);//debugging
	if (UNIQUE == true) {
		for (int name = (int)*f_number_of_list_genes-1; name >= 0; name--) {//if I count a gene twice, it is most likely from the front 
	//		printf("array[%u] = %s with %s\n",name, array[name], NAME);
			if (strcmp(array[name], NAME) == 0) {
				return array;//if this name is already in the list, don't add it
			}
		}
	}
	const size_t NAME_SIZE = strlen(NAME) + 1;//was "strlen(NAME)*+2*sizeof(char*)", which is 0 for an empty name
	*size = ((size_t)*f_number_of_list_genes + 1) * sizeof(char*);//one pointer per name
	array = realloc(array, *size);
	if (array == NULL) {
		fprintf(stderr, "realloc of array failed.\n");
		exit(EXIT_FAILURE);
	}
	array[*f_number_of_list_genes] = malloc(NAME_SIZE);
	if (array[*f_number_of_list_genes] == NULL) {
		fprintf(stderr, "alloc of array[%u] failed.\n", *f_number_of_list_genes);
		exit(EXIT_FAILURE);
	}
	strcpy(array[*f_number_of_list_genes], NAME);
	*f_number_of_list_genes += 1;
	return array;
}

void * clear_names(char **restrict array, size_t *restrict size, unsigned int *restrict number_of_annotated_genes) {//empty the list of names
	for (unsigned int name = 0; name < *number_of_annotated_genes; name++) {
		free(array[name]); array[name] = NULL;//this is missing something
		if (array[name] != NULL) {
			printf("free of array[%u] failed.\n", name);
		}
	}
	*number_of_annotated_genes = 0;
	*size = 1;
	array = realloc(array, *size);
	if (array == NULL) {
		fprintf(stderr, "realloc of array failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return array;
}

void safe_string_copy ( char *restrict DESTINATION, const char *restrict SOURCE) {
//this doesn't check sizes: only use it when DESTINATION is known to be big enough for SOURCE, otherwise use bounded_string_copy
	strcpy(DESTINATION,SOURCE);//already ensures NULL terminator
}

void bounded_string_copy (char *restrict DESTINATION, const size_t DESTINATION_SIZE, const char *restrict SOURCE) {
//copy into a fixed-size buffer, and stop instead of writing past its end
	const size_t LENGTH = strlen(SOURCE);
	if (LENGTH >= DESTINATION_SIZE) {
		fprintf(stderr, "\"%s\" is %zu characters long, but at most %zu characters fit here.\n", SOURCE, LENGTH, DESTINATION_SIZE - 1);
		exit(EXIT_FAILURE);
	}
	memcpy(DESTINATION, SOURCE, LENGTH + 1);
}

void bounded_string_append (char *restrict destination, const size_t DESTINATION_SIZE, const char *restrict NEW) {
//like strcat, but a string that would be too long ends with "..." instead of overflowing the buffer
	const size_t USED = strlen(destination);
	const size_t LENGTH = strlen(NEW);
	if (USED + LENGTH < DESTINATION_SIZE) {
		memcpy(destination + USED, NEW, LENGTH + 1);
	} else if ((DESTINATION_SIZE >= 4) && (USED >= 3) && (strcmp(destination + USED - 3, "...") != 0)) {
		const size_t END = USED + 3 < DESTINATION_SIZE ? USED : DESTINATION_SIZE - 4;
		memcpy(destination + END, "...", 4);
	}
}

void chomp (char *restrict string) {//remove the line ending, "\n" or "\r\n", like perl's chomp
	size_t length = strlen(string);
	while ((length > 0) && ((string[length-1] == '\n') || (string[length-1] == '\r'))) {
		string[length-1] = '\0';
		length--;
	}
}
void removeSubstring(char *restrict s, const char *restrict toremove) {//remove a small string from a larger string
//http://stackoverflow.com/questions/4833347/removing-substring-from-a-string user411313
	while( s==strstr(s,toremove) ) {
		memmove(s,s+strlen(toremove),1+strlen(s+strlen(toremove)));
	}
}

unsigned short int get_annotation_type (const char *restrict line) {//determines file type of input, e.g. refFlat or Gencode GTF
//return 1: refFlat
//return 2: Gencode chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
	if (line[0] == '#') {
		return 0;//comment, assumes comments don't exist in refFlat (I don't see them)
	}
	if (strstr(line, "	Ensembl	chromosome	") != NULL) {
		return 0;
	}
//	printf("@ %s line %u\n", __FILE__, __LINE__);
	char *restrict line_copy = strdup(line);
	if (line_copy == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	chomp(line_copy);//get rid of newline
	char *restrict temporary_char_string,  *saveptr;//necessary for strtok_r
	temporary_char_string = strtok_r(line_copy, "\t ", &saveptr);
	if (temporary_char_string == NULL) {//blank line, keep looking
		free(line_copy); line_copy = NULL;
		return 0;
	}
	unsigned short int number_of_columns = 1, strand_definition = 0;
//	bool column1_is_chr = false;

//	if ((number_of_columns == 1) && (strstr(temporary_char_string,"chr") != NULL)) {
//		column1_is_chr = true;
//	}
	while (temporary_char_string != NULL) {
		if ((temporary_char_string[0] == '-') || (temporary_char_string[0] == '+')) {
			strand_definition = number_of_columns;
		}
		temporary_char_string = strtok_r(NULL, "\t ", &saveptr);
		if (temporary_char_string == NULL) {
			break;
		}
		number_of_columns++;
	}
	free(line_copy); line_copy = NULL;
//	printf("strand definition = %u\n", strand_definition);
	if (strand_definition == 0) {
		fprintf(stderr, "strand_definition isn't defined @ %s line %u\n", __FILE__, __LINE__);
		fprintf(stderr, "This happened for %s\n", line);
		exit(EXIT_FAILURE);
	}
	if (strand_definition == 4) {
		return 1;//refFlat
	} else if (strand_definition == 7) {
		return 2;//GTF & GFF
	} else {
		fprintf(stderr, "Couldn't recognize annotation file format with strand_definition = %u @ %s line %u\n", strand_definition, __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
}

bool get_gtf_attribute (const char *restrict line, const char *restrict KEY, char *restrict value, const size_t VALUE_SIZE) {
//find an attribute in the 9th column of a GTF line (gene_name "Xkr4";) or of a GFF3 line (gene_name=Xkr4;)
//returns false if the line doesn't have that attribute
	const size_t KEY_LENGTH = strlen(KEY);
	const char *position = line;
	while ((position = strstr(position, KEY)) != NULL) {
		const bool STARTS_FIELD = (position == line) || (position[-1] == ' ') || (position[-1] == '\t') || (position[-1] == ';');//don't match "gene_name" inside "other_gene_name"
		const char *after = position + KEY_LENGTH;
		position = after;
		if ((STARTS_FIELD == false) || ((*after != ' ') && (*after != '='))) {
			continue;
		}
		while ((*after == ' ') || (*after == '=')) {
			after++;
		}
		const char *end;
		if (*after == '"') {//quoted, GTF style
			after++;
			end = strchr(after, '"');
			if (end == NULL) {
				continue;
			}
		} else {//unquoted, GFF3 style, ends at the next ';' or the end of the line
			end = after;
			while ((*end != '\0') && (*end != ';') && (*end != '\n') && (*end != '\r') && (*end != '\t')) {
				end++;
			}
		}
		const size_t LENGTH = (size_t)(end - after);
		if ((LENGTH == 0) || (LENGTH >= VALUE_SIZE)) {
			continue;
		}
		memcpy(value, after, LENGTH);
		value[LENGTH] = '\0';
		return true;
	}
	return false;
}

typedef struct {
	char *restrict gene_name, *restrict chromosome;//, transcript[48];//1,2, M, X, Y etc.
	unsigned int gene_start, gene_end;//, code_start, code_end;
//	unsigned short int number_of_exons;
	short int sense;//-1 or +1
} ANNOTATION_LINE_STRUCT;

ANNOTATION_LINE_STRUCT read_gene_annotation_line (const char *restrict line, const unsigned short int TYPE) {
	ANNOTATION_LINE_STRUCT gene;
	gene.gene_name = NULL;//the caller frees gene_name & chromosome, also when this line is ignored
	gene.chromosome = NULL;
	gene.gene_start = 0;
	gene.gene_end = 0;
	gene.sense = 0;//0 tells the caller to ignore this line
	char *restrict copy_line = strdup(line);
	if (copy_line == NULL) {
		fprintf(stderr, "strdup failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	chomp(copy_line);
	if (TYPE == 1) {//refFlat format
		char *restrict column, *column_pointer;
		column = strtok_r(copy_line,"\t ",&column_pointer);//http://stackoverflow.com/questions/15961253/c-correct-usage-of-strtok-r
		if (column == NULL) {//blank line
			free(copy_line); copy_line = NULL;
			return gene;
		}
		gene.gene_name = strdup(column);
//		printf("@line %u	gene %s ID'd\n", __LINE__, gene.gene_name);
//2nd column: transcript
		column = strtok_r(NULL,"\t ",&column_pointer);
		if (column == NULL) {//error checking
			fprintf(stderr, "2nd column of annotation file is NULL for line %s @ %s line %u.\n",copy_line, __FILE__, __LINE__);//DEBUGGING
			exit(EXIT_FAILURE);
		}
//				strncpy(gene.transcript,column,sizeof(gene.transcript));
//3rd column: chromosome
		column = strtok_r(NULL,"\t ",&column_pointer);
		if (column == NULL) {//error checking
			fprintf(stderr, "3rd column of annotation file is NULL for line %s @ %s line %u.\n",copy_line,__FILE__,__LINE__);//DEBUGGING
			exit(EXIT_FAILURE);
		}
		removeSubstring(column,"chr");
		gene.chromosome = strdup(column);
//4th column: +/- sense
		column = strtok_r(NULL,"\t ",&column_pointer);
		if (column == NULL) {//error checking
			fprintf(stderr, "4th column of annotation file is NULL for line %s @ %s line %u.\n",copy_line,__FILE__,__LINE__);//DEBUGGING
			exit(EXIT_FAILURE);
		}
		gene.sense = 1;//1 if + sense, -1 if - sense
		if (strcmp(column,"-") == 0) {
			gene.sense = -1;
		}
//5th column: gene start
		column = strtok_r(NULL,"\t ",&column_pointer);
		if (column == NULL) {//error checking
			fprintf(stderr, "5th column of annotation file is NULL for line %s @ %s line %u.\n",copy_line,__FILE__,__LINE__);//DEBUGGING
			exit(EXIT_FAILURE);
		}
		gene.gene_start = (unsigned int)safe_strtoul(column);
//6th column: gene end
		column = strtok_r(NULL,"\t ",&column_pointer);
		if (column == NULL) {//error checking
			fprintf(stderr, "6th column of annotation file is NULL for line %s @ %s line %u.\n",copy_line,__FILE__,__LINE__);//DEBUGGING
			exit(EXIT_FAILURE);
		}
		gene.gene_end = (unsigned int)safe_strtoul(column);
	} else if (TYPE == 2) {//GTF format (GFF3)
//		printf("@ %s line %u\n", __FILE__, __LINE__);
		if (line[0] == '#') {//comment; a '#' elsewhere on the line is part of an attribute
			free(copy_line); copy_line = NULL;
			return gene;
		}
//the gene's name is read from the whole line before strtok_r cuts it apart
//gene_name is preferred (Ensembl, Gencode), then gene_id (e.g. UCSC's mm9 GTF), then GFF3's Name
		const size_t GENE_NAME_SIZE = strlen(copy_line) + 1;//an attribute can't be longer than the line that it's on
		gene.gene_name = malloc(GENE_NAME_SIZE);
		if (gene.gene_name == NULL) {
			fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if ((get_gtf_attribute(copy_line, "gene_name", gene.gene_name, GENE_NAME_SIZE) == false) &&
			(get_gtf_attribute(copy_line, "gene_id",   gene.gene_name, GENE_NAME_SIZE) == false) &&
			(get_gtf_attribute(copy_line, "gene",      gene.gene_name, GENE_NAME_SIZE) == false) &&
			(get_gtf_attribute(copy_line, "Name",      gene.gene_name, GENE_NAME_SIZE) == false)) {
			free(copy_line); copy_line = NULL;//e.g. a GFF3 "chromosome" or "region" line, which isn't a gene
			return gene;
		}
		char *restrict column, *column_pointer;
		column = strtok_r(copy_line,"\t ",&column_pointer);//1st column, chromosome
		if (column == NULL) {//blank line
			free(copy_line); copy_line = NULL;
			return gene;
		}
		removeSubstring(column,"chr");
		gene.chromosome = strdup(column);
		column = strtok_r(NULL,"\t ",&column_pointer);//2nd column, "source", discard
		column = strtok_r(NULL,"\t ",&column_pointer);//3rd column, "type", discard
		if (column == NULL) {//error checking
			fprintf(stderr, "3rd column of annotation file is NULL for line %s @ %s line %u.\n", line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		column = strtok_r(NULL,"\t ",&column_pointer);//4th column, "start"
		if (column == NULL) {//error checking
			fprintf(stderr, "4th column of annotation file is NULL for line %s @ %s line %u.\n", line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		gene.gene_start = (unsigned int)safe_strtoul(column);
		column = strtok_r(NULL,"\t ",&column_pointer);//5th column, "end"
		if (column == NULL) {//error checking
			fprintf(stderr, "5th column of annotation file is NULL for line %s @ %s line %u.\n", line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		gene.gene_end = (unsigned)safe_strtoul(column);
		column = strtok_r(NULL,"\t ",&column_pointer);//6th col, score, discard
		column = strtok_r(NULL,"\t ",&column_pointer);//7th col, sense
		if (column == NULL) {//error checking
			fprintf(stderr, "7th column of annotation file is NULL for line %s @ %s line %u.\n", line, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if (strcmp(column, "-") == 0) {
			gene.sense = -1;
		} else {
			gene.sense = 1;//1 if + sense, -1 if - sense
		}
	} else {
		free(copy_line); copy_line = NULL;
		fprintf(stderr, "I got a type = %u, which is not a recognized type. This is @ %s line %u\n", TYPE,__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	free(copy_line); copy_line = NULL;
	return gene;
}
#endif
