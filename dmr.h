#ifndef DEFIANT_DMR_H
#define DEFIANT_DMR_H
//requires annotation.h and fdr.h to be included first
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <float.h>

typedef struct {//store nucleotide information
	unsigned int nucleotide, methylated_C, coverage;
} BISULFITE_DATA_STRUCT;

typedef struct {
	size_t methylated_C, coverage;
} DMR_struct;

void * malloc_nucleotides (DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES) {
	array = malloc((size_t)(NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1]) * sizeof(DMR_struct));
	if (array == NULL) {
		fprintf(stderr, "alloc failed @ %s line %u.\n", __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	for (unsigned short int set = 0; set < 2; set++) {
		array[set] = malloc( NO_OF_REPLICATES[set]*sizeof(DMR_struct));
		if (array[set] == NULL) {
			fprintf(stderr, "alloc failed @ %s line %u.\n", __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			array[set][replicate] = malloc(sizeof(DMR_struct));
			if (array[set][replicate] == NULL) {
				fprintf(stderr, "alloc failed @ %s line %u.\n", __FILE__, __LINE__);
				perror("");
				exit(EXIT_FAILURE);
			}
		}
	}
	return array;
}

void * free_nucleotides (DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES) {
	array = realloc(array, (size_t)(NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1]) * sizeof(DMR_struct));
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			free(array[set][replicate]); 	array[set][replicate] = NULL;
		}
		array[set] = realloc(array[set], NO_OF_REPLICATES[set]*sizeof(DMR_struct));
	}
	return array;
}
void * final_free_nucleotides (DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES) {
	array = realloc(array, (size_t)(NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1]) * sizeof(DMR_struct));
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			free(array[set][replicate]); 	array[set][replicate] = NULL;
		}
		free(array[set]); array[set] = NULL;
	}
	free(array); array = NULL;
	return array;
}

void * add_nucleotides (DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES, const unsigned int INDEX, DMR_struct **restrict single_CpN) {
//allocates extra memory for an additional CpN, and then writes to RAM
//array (2 sets) and array[set] (one pointer per replicate) were sized by malloc_nucleotides, only array[set][replicate] grows
	if (array == NULL) {
		fprintf(stderr, "add_nucleotides needs an array from malloc_nucleotides @ %s line %u.\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			array[set][replicate]							= realloc(array[set][replicate], ((size_t)1+INDEX)*sizeof(DMR_struct));
			if (array[set][replicate] == NULL) {
				fprintf(stderr, "realloc failed @ %s line %u.\n", __FILE__, __LINE__);
				perror("");
				exit(EXIT_FAILURE);
			}
			array[set][replicate][INDEX].methylated_C	= single_CpN[set][replicate].methylated_C;
			array[set][replicate][INDEX].coverage		= single_CpN[set][replicate].coverage;
//			printf("array[%u][%u][%u] = %zu/%zu @ %s line %u\n", set, replicate, INDEX, array[set][replicate][INDEX].methylated_C, array[set][replicate][INDEX].coverage, __FILE__, __LINE__);
		}
	}
	return array;
}

void * clear_nucleotides (DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES) {
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			free(array[set][replicate]);	array[set][replicate] = NULL;
			array[set][replicate] = realloc(array[set][replicate], NO_OF_REPLICATES[set]*sizeof(DMR_struct));
		}
		array[set] = realloc(array[set], NO_OF_REPLICATES[set]*sizeof(DMR_struct) );
	}
	array = realloc(array, (size_t)(NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1]) * sizeof(DMR_struct));
	return array;
}

void print_methylation_percentages (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, FILE *restrict answers_filehandle) {
//this assumes coverages are all > 0
	for (unsigned short int set_loop = 0; set_loop < 2; set_loop++) {
		fprintf(answers_filehandle,"\t[%.1f",100.0*(double)dmr[set_loop][0].methylated_C / (double)dmr[set_loop][0].coverage);
		for (unsigned short int r = 1; r < NO_OF_REPLICATES[set_loop]; r++) {
			fprintf(answers_filehandle,",%.1f",100.0*(double)dmr[set_loop][r].methylated_C / (double)dmr[set_loop][r].coverage);
		}
		fprintf(answers_filehandle,"]");
	}
}

void Rprint_methylation_percentages (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, FILE *restrict r, const unsigned short int SET) {
	bool first = true;
	for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[SET]; replicate++) {
		if (dmr[SET][replicate].coverage > 0) {
			if (first == true) {//first will not have comma
				fprintf(r,"%zu/%zu",dmr[SET][replicate].methylated_C, dmr[SET][replicate].coverage);
				first = false;
			} else {
				fprintf(r,",%zu/%zu", dmr[SET][replicate].methylated_C, dmr[SET][replicate].coverage);
			}
		}
	}
}

double single_set_DMR_mean (DMR_struct **restrict DMR, const unsigned short int *restrict NO_OF_REPLICATES, const unsigned short int SET) {
	size_t methylated_C_sum = 0, fcoverage_sum = 0;
	for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[SET]; replicate++) {
//		if (DMR[SET][replicate].methylated_C > DMR[SET][replicate].coverage) {
//			printf("@line %u	methylated_C = %zu > coverage = %zu, not possible.\n", __LINE__, DMR[SET][replicate].methylated_C, DMR[SET][replicate].coverage);
//			exit(EXIT_FAILURE);
//		}
		if (DMR[SET][replicate].coverage > 0) {
			methylated_C_sum += DMR[SET][replicate].methylated_C;
			fcoverage_sum += DMR[SET][replicate].coverage;
		}
	}
	if (fcoverage_sum == 0) {
		fprintf(stderr, "@ file %s line %u you tried to divide by 0 when taking a set mean.\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	return (double) methylated_C_sum / (double)fcoverage_sum;
}

/*const unsigned int coverage_sum (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, const unsigned short int SET) {
	unsigned int sum_coverage = 0;
	for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[SET]; replicate++) {
		if (dmr[SET][replicate].coverage > 0) {
			sum_coverage++;
		}
	}
	return sum_coverage;
}*/

double Pvalue (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, const bool PAIRED) {//calculate a p-value based on an array
	size_t methylated_C_1 = 0, methylated_C_2 = 0, coverage1 = 0, coverage2 = 0;
	for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[0]; replicate++) {
/*		if ((dmr[0][replicate].coverage == 0) || (dmr[0][replicate].methylated_C > dmr[0][replicate].coverage)) {
			printf("got 5mC/C = %zu / %zu @ %s line %u\n", dmr[0][replicate].methylated_C , dmr[0][replicate].coverage, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}*/
		methylated_C_1 += dmr[0][replicate].methylated_C;
		coverage1 += dmr[0][replicate].coverage;
	}
//	printf("@ %s line %u sum = %zu/%zu\n", __FILE__, __LINE__, methylated_C_1, coverage1);
	for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[1]; replicate++) {
/*		if ((dmr[1][replicate].coverage == 0) || (dmr[1][replicate].methylated_C > dmr[1][replicate].coverage)) {
			printf("got 5mC/C = %zu / %zu @ %s line %u\n", dmr[1][replicate].methylated_C , dmr[1][replicate].coverage, __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}*/
		methylated_C_2 += dmr[1][replicate].methylated_C;
		coverage2 += dmr[1][replicate].coverage;
	}
	if ((coverage1 == 0) || (coverage2 == 0)) {
		return 1;
	}
//	printf("@ %s line %u sum = %zu/%zu\n", __FILE__, __LINE__, methylated_C_2, coverage2);
	if ((NO_OF_REPLICATES[0] == 1) || (NO_OF_REPLICATES[1] == 1)) {//Fisher's exact test
//		printf("A = %u	B = %u	C = %u	D = %u, coverage1 = %u,	coverage2 = %u\n", A, B, C, D, coverage1, coverage2);
		const double P = exp(lgamma((double)(methylated_C_1+methylated_C_2+1)) + lgamma((double)(coverage1 - methylated_C_1+coverage2 - methylated_C_2+1)) + lgamma((double)(coverage1+1)) + lgamma((double)(coverage2+1)) - lgamma((double)(methylated_C_1+1)) - lgamma((double)(methylated_C_2+1)) - lgamma((double)(coverage1 - methylated_C_1+1)) - lgamma((double)(coverage2 - methylated_C_2+1)) - lgamma((double)(coverage1 + coverage2 +1)));
		return P > 1.0 ? 1.0 : P;//rounding can give 1.0000000000000002, and p_adjust stops if p > 1
	}
	const double FMEAN1 = (double)methylated_C_1 / (double)coverage1, FMEAN2 = (double)methylated_C_2 / (double)coverage2;
	if (FMEAN1 == FMEAN2) {
		return 1.0;//if the means are equal, the p-value is 1, leave the function
	}
	double unbiased_sample_variance1 = 0.0, unbiased_sample_variance2 = 0.0;
	double T_STATISTIC = DBL_MAX;
	double DEGREES_OF_FREEDOM = DBL_MAX;
	if (PAIRED == false) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[0]; replicate++) {//1st part of added unbiased_sample_variance
			unbiased_sample_variance1 += ((double)dmr[0][replicate].methylated_C/(double)dmr[0][replicate].coverage-FMEAN1)*((double)dmr[0][replicate].methylated_C/(double)dmr[0][replicate].coverage-FMEAN1)*(double)dmr[0][replicate].coverage;
		}
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[1]; replicate++) {
			unbiased_sample_variance2 += ((double)dmr[1][replicate].methylated_C/(double)dmr[1][replicate].coverage-FMEAN2)*((double)dmr[1][replicate].methylated_C/(double)dmr[1][replicate].coverage-FMEAN2)*(double)dmr[1][replicate].coverage;
		}
	//	printf("unbiased_sample_variance1 = %lf\tunbiased_sample_variance2 = %lf\n",unbiased_sample_variance1,unbiased_sample_variance2);//DEBUGGING
		unbiased_sample_variance1 = unbiased_sample_variance1/(double)(coverage1-1);
		unbiased_sample_variance2 = unbiased_sample_variance2/(double)(coverage2-1);
		
		DEGREES_OF_FREEDOM = pow((unbiased_sample_variance1/NO_OF_REPLICATES[0]+unbiased_sample_variance2/NO_OF_REPLICATES[1]),2.0)//numerator
		 /
		(
			(unbiased_sample_variance1*unbiased_sample_variance1)/(NO_OF_REPLICATES[0]*NO_OF_REPLICATES[0]*(NO_OF_REPLICATES[0]-1))+
			(unbiased_sample_variance2*unbiased_sample_variance2)/(NO_OF_REPLICATES[1]*NO_OF_REPLICATES[1]*(NO_OF_REPLICATES[1]-1))
		);
		T_STATISTIC = (FMEAN1-FMEAN2)/sqrt(unbiased_sample_variance1/NO_OF_REPLICATES[0]+unbiased_sample_variance2/NO_OF_REPLICATES[1]);
		
	} else {//PAIRED = true
		double X_mean = 0.0;
		double D_mean = 0.0;
		double xd = 0.0;
		
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[0]; r++) {//for each replicate
			X_mean += (1.0/(double)dmr[0][r].coverage + 1.0/(double)dmr[1][r].coverage);
			const double PERCENT1 = (double)(dmr[0][r].methylated_C) / (double)dmr[0][r].coverage;
			const double PERCENT2 = (double)dmr[1][r].methylated_C / (double)dmr[1][r].coverage;
			D_mean += (PERCENT2 - PERCENT1);
			xd += (PERCENT2 - PERCENT1) * (double)( dmr[0][r].coverage + dmr[1][r].coverage);
		}
		if (xd == 0.0) {
			return 1.0;
		}
		X_mean /= NO_OF_REPLICATES[0];
		D_mean /= NO_OF_REPLICATES[0];
//		printf("X_mean = %lf\n", X_mean);
//		printf("D_mean = %lf\n", D_mean);
		
		double b_numerator = 0.0;
		double b_denominator = 0.0;
		double m_mean = 0.0;
		size_t weight_sum = 0;
		
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[0]; r++) {//for each replicate
			weight_sum += ( dmr[0][r].coverage + dmr[1][r].coverage );
			const double Xi = 1.0 / ((double)dmr[0][r].coverage) + 1.0 / ((double)dmr[1][r].coverage);
			const double PERCENT1 = (double)(dmr[0][r].methylated_C) / (double)dmr[0][r].coverage;
			const double PERCENT2 = (double)dmr[1][r].methylated_C / (double)dmr[1][r].coverage;
			const double di = PERCENT2 - PERCENT1;
			const double mi = (di - D_mean) * (di - D_mean);
			m_mean += mi;
			b_numerator += (Xi - X_mean) * mi;
			b_denominator += (Xi - X_mean) * (Xi - X_mean);
		}
		xd /= (double)weight_sum;
		m_mean /= NO_OF_REPLICATES[0];
//		printf("m_mean = %lf\n", m_mean);
		double var = 0.0;
		
		if (b_denominator != 0) {
			const double b_hat = (double) b_numerator / b_denominator;
			const double a_hat = m_mean - b_hat * X_mean;
			
			for (unsigned short int r = 0; r < NO_OF_REPLICATES[0]; r++) {//for each replicate
				const double Xi = 1.0 / ((double)dmr[0][r].coverage) + 1.0 / ((double)dmr[1][r].coverage);
				var += (pow((double)(dmr[0][r].coverage + dmr[1][r].coverage), 2.0)) * (Xi*b_hat + a_hat);
			}
			T_STATISTIC = (FMEAN2- FMEAN1) / (var / sqrt((double)weight_sum) );
		} else {
			double sd = 0.0;
//			printf("xd = %lf\n", xd);
			for (unsigned short int r = 0; r < NO_OF_REPLICATES[0]; r++) {//for each replicate
				const double PERCENT1 = (double)dmr[0][r].methylated_C / (double)dmr[0][r].coverage;
				const double PERCENT2 = (double)dmr[1][r].methylated_C / (double)dmr[1][r].coverage;
//				printf("pct1 = %lf	pct2 = %lf\n", PERCENT1, PERCENT2);
				sd += pow(xd - (PERCENT2 - PERCENT1), 2.0);
//				printf("r = %u	sd = %lf\n", r, sd);
			}
//			printf("sd = %lf\n", sd);
			sd = sqrt(sd / (NO_OF_REPLICATES[0] - 1));
			if (sd == 0.0) {
				return 1.0;
			}
			T_STATISTIC = xd / (sd / sqrt(NO_OF_REPLICATES[0]));
		}
		DEGREES_OF_FREEDOM = NO_OF_REPLICATES[0] - 1;
	}
//	printf("t = %lf\n", T_STATISTIC);
	const double a = DEGREES_OF_FREEDOM/2;
	double value = DEGREES_OF_FREEDOM/(T_STATISTIC*T_STATISTIC+DEGREES_OF_FREEDOM);
//	printf("%e	%e	%e\n", T_STATISTIC, DEGREES_OF_FREEDOM, value);
	if ((isinf(value) != 0) || (isnan(value) != 0)) {
		return 1.0;
	}
	if ((isinf(value) != 0) || (isnan(value) != 0)) {
		return 1.0;
	}

/*  Purpose:

    BETAIN computes the incomplete Beta function ratio.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 November 2010

  Author:

    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
    C version by John Burkardt.

  Reference:

    KL Majumder, GP Bhattacharjee,
    Algorithm AS 63:
    The incomplete Beta Integral,
    Applied Statistics,
    Volume 22, Number 3, 1973, pages 409-411.

  Parameters:
https://www.jstor.org/stable/2346797?seq=1#page_scan_tab_contents
    Input, double X, the argument, between 0 and 1.

    Input, double P, Q, the parameters, which
    must be positive.

    Input, double BETA, the logarithm of the complete
    beta function.

    Output, int *IFAULT, error flag.
    0, no error.
    nonzero, an error occurred.

    Output, double BETAIN, the value of the incomplete
    Beta function ratio.
*/
	const double beta = lgamma(a)+0.57236494292470009-lgamma(a+0.5);
	const double acu = 0.1E-14;
  double ai;
  double cx;
  int indx;
  int ns;
  double pp;
  double psq;
  double qq;
  double rx;
  double temp;
  double term;
  double xx;

//  ifault = 0;
//Check the input arguments.
  if ( (a <= 0.0)) {// || (0.5 <= 0.0 )){
//    *ifault = 1;
//    return value;
  }
  if ( value < 0.0 || 1.0 < value )
  {
//    *ifault = 2;
    return value;
  }
/*
  Special cases.
*/
  if ( value == 0.0 || value == 1.0 )   {
    return value;
  }
  psq = a + 0.5;
  cx = 1.0 - value;

  if ( a < psq * value )
  {
    xx = cx;
    cx = value;
    pp = 0.5;
    qq = a;
    indx = 1;
  }
  else
  {
    xx = value;
    pp = a;
    qq = 0.5;
    indx = 0;
  }

  term = 1.0;
  ai = 1.0;
  value = 1.0;
  ns = ( int ) ( qq + cx * psq );
/*
  Use the Soper reduction formula.
*/
  rx = xx / cx;
  temp = qq - ai;
  if ( ns == 0 )
  {
    rx = xx;
  }

  for ( ; ; )
  {
    term = term * temp * rx / ( pp + ai );
    value = value + term;;
    temp = fabs ( term );

    if ( temp <= acu && temp <= acu * value )
    {
      value = value * exp ( pp * log ( xx ) 
      + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

      if ( indx )
      {
        value = 1.0 - value;
      }
      break;
    }

    ai = ai + 1.0;
    ns = ns - 1;

    if ( 0 <= ns )
    {
      temp = qq - ai;
      if ( ns == 0 )
      {
        rx = xx;
      }
    }
    else
    {
      temp = psq;
      psq = psq + 1.0;
    }
  }
  return value;
}


void zero_DMR_struct (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES) {
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			dmr[set][replicate].methylated_C = 0;
			dmr[set][replicate].coverage = 0;
		}
	}
}

typedef struct {//this allows me to keep RAM use down by saving chromosome changes instead of chromosome for each line read
	unsigned int index;
	char *restrict chromosome;
} CHROMOSOME_CHANGES_STRUCT;

typedef struct {
	char *restrict gene_name;
	unsigned int gene_start, gene_end;
	short int sense;//-1 or +1
} GENES_STRUCT;

void print_DMR_struct ( DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES) {
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
			printf("dmr[%u][%u].methylated_C = %zu	dmr[%u][%u].coverage = %zu\n", set, replicate, dmr[set][replicate].methylated_C, set, replicate, dmr[set][replicate].coverage);
		}
	}
}/*debugging*/


double mean_percent_difference (DMR_struct **restrict DMR, const unsigned short int *restrict NO_OF_REPLICATES) {
	size_t methylated_C_sum[] = {0,0};
	size_t coverage_sum[] = {0,0};
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
/*			if (DMR[set][replicate].methylated_C > DMR[set][replicate].coverage) {
				printf("@line %u	methylated_C = %zu > coverage = %zu, not possible.\n", __LINE__, DMR[set][replicate].methylated_C, DMR[set][replicate].coverage);
				exit(EXIT_FAILURE);
			}*/
			methylated_C_sum[set] += DMR[set][replicate].methylated_C;
			coverage_sum[set] += DMR[set][replicate].coverage;
		}
	}
	if ((coverage_sum[0] == 0) || (coverage_sum[1] == 0)) {
		return 0;//actually undefined
	}
	return (double)((double)methylated_C_sum[1]/(double)coverage_sum[1] - (double)methylated_C_sum[0]/(double)coverage_sum[0]);
}

unsigned int promoter_cutoff = 10000;//used in assigning
unsigned int number_of_chromosomes_in_GENES = 0;//declared outside main so I can use it in functions
GENES_STRUCT *restrict GENES = 0;//declare struct to hold gene data
CHROMOSOME_CHANGES_STRUCT *restrict gene_chromosome_changes = NULL;//declared outside main so I can use it in functions
unsigned int GENES_chromosome_change_index = 0, number_of_annotated_genes = 0;
bool gene_annotation_defined = false, printing_CpN_list = false;

#define GENE_STRING_SIZE 4096//size of the strings that identify_gene writes gene names into

bool list_contains (const char *restrict LIST, const char *restrict NAME, const char *restrict DELIMITER) {
//is NAME one of the entries of the delimited LIST?  strstr alone would find "Abc1" inside "Abc10"
	const size_t NAME_LENGTH = strlen(NAME), DELIMITER_LENGTH = strlen(DELIMITER);
	const char *position = LIST;
	while ((position = strstr(position, NAME)) != NULL) {
		const bool STARTS = (position == LIST) || ((size_t)(position - LIST) >= DELIMITER_LENGTH && strncmp(position - DELIMITER_LENGTH, DELIMITER, DELIMITER_LENGTH) == 0);
		const bool ENDS = (position[NAME_LENGTH] == '\0') || (strncmp(position + NAME_LENGTH, DELIMITER, DELIMITER_LENGTH) == 0);
		if ((STARTS == true) && (ENDS == true)) {
			return true;
		}
		position++;
	}
	return false;
}

void append_string (char *restrict destination, const size_t DESTINATION_SIZE, const char *restrict NEW) {
	if (strcmp(destination,"n/a") == 0) {//if the variable list is empty,
		destination[0] = '\0';
		bounded_string_append(destination, DESTINATION_SIZE, NEW);
	} else if (list_contains(destination, NEW, ",") == false) {
		bounded_string_append(destination, DESTINATION_SIZE, ",");
		bounded_string_append(destination, DESTINATION_SIZE, NEW);
	}
}

void identify_gene (const unsigned int START, const unsigned int END, char *restrict gene_string, char *restrict intergene_string, char *restrict list_gene_string, const size_t STRING_SIZE, unsigned int *restrict min_gene_search_index, unsigned int *restrict max_gene_search_index, const unsigned int first_nucleotide) {//given nucleotide and globally accessed lists, how can I identify this nucleotide's genomic region?
//gene_string, intergene_string & list_gene_string must each have room for STRING_SIZE characters
//GENES[*min_gene_search_index] to GENES[*max_gene_search_index] are this chromosome's genes, sorted by gene_start
	if (number_of_annotated_genes == 0) {
		return;
	}
	const unsigned int LAST_GENE = number_of_annotated_genes - 1;//the last valid index of GENES
	if (*min_gene_search_index > *max_gene_search_index) {
		*min_gene_search_index = 0;
		*max_gene_search_index = LAST_GENE;
	}
	if (*max_gene_search_index > LAST_GENE) {
		*max_gene_search_index = LAST_GENE;
	}
	if (END < first_nucleotide) {//is this nucleotide in the beginning of the chromosome?
		snprintf(intergene_string, STRING_SIZE, "Chr_Start-%u-%s",GENES[*min_gene_search_index].gene_start - END, GENES[*min_gene_search_index].gene_name);
		if ((GENES[*min_gene_search_index].sense == 1) && ((GENES[*min_gene_search_index].gene_start - END) < promoter_cutoff)) {
			bounded_string_copy(list_gene_string, STRING_SIZE, GENES[*min_gene_search_index].gene_name);
		}
	} else if (START > GENES[*max_gene_search_index].gene_end) {//is this nucleotide in the end of the chromosome?
		snprintf(intergene_string, STRING_SIZE, "%s-%u-End_Chr",GENES[*max_gene_search_index].gene_name,START-GENES[*max_gene_search_index].gene_end);
	}
	const unsigned int FIRST_GENE_SEARCHED = *min_gene_search_index;
	bool search_start_moved = false;
	for (unsigned int i = FIRST_GENE_SEARCHED; i <= *max_gene_search_index; i++) {
		if ((GENES[i].gene_start <= START) && (END <= GENES[i].gene_end)) {//is this nucleotide in this gene?
			append_string(gene_string, STRING_SIZE, GENES[i].gene_name);
			if (strcmp(list_gene_string, "n/a") == 0) {
				bounded_string_copy(list_gene_string, STRING_SIZE, GENES[i].gene_name);
			} else {
				bounded_string_append(list_gene_string, STRING_SIZE, "; ");//don't overwrite gene name
				bounded_string_append(list_gene_string, STRING_SIZE, GENES[i].gene_name);
			}
//the next DMR is further along the chromosome, so its search can start at the FIRST gene that this DMR matched.
//Starting at the last match would skip a long gene that contains both this DMR and the next one, when shorter genes start inside it
			if (search_start_moved == false) {
				*min_gene_search_index = i;
				search_start_moved = true;
			}
		} else if ((i < *max_gene_search_index) && (GENES[i].gene_end < START) && (END < GENES[i+1].gene_start)) {//intergenic
			const unsigned int intergenic_start_nucleotide	= START - GENES[i].gene_end;
			const unsigned int intergenic_end_nucleotide		= GENES[i+1].gene_start - END;
			snprintf(intergene_string, STRING_SIZE, "%s-%u,%u-%s",GENES[i].gene_name, intergenic_start_nucleotide, intergenic_end_nucleotide,GENES[i+1].gene_name);
			if ((GENES[i].sense == -1) && (intergenic_start_nucleotide <= promoter_cutoff)) {//the gene is transcribed 3'-5', and within promoter cutoff
				bounded_string_copy(list_gene_string, STRING_SIZE, GENES[i].gene_name);
			}
			if ((GENES[i+1].sense == 1) && (intergenic_end_nucleotide <= promoter_cutoff)) {//the gene is transcribed 5'-3', and within promoter cutoff
				bounded_string_copy(list_gene_string, STRING_SIZE, GENES[i+1].gene_name);
			}
			if (search_start_moved == false) {
				*min_gene_search_index = i;
				search_start_moved = true;
			}
		} else if (END < GENES[i].gene_start) {//I've passed what this gene can possibly be
			return;
		}
	}
//a DMR that straddles a gene's border is neither inside a gene nor between genes, and keeps "n/a"
}

void add_2nd_dmr_to_1st_dmr(DMR_struct **restrict dmr1, DMR_struct **restrict dmr2, const unsigned short int *restrict NO_OF_REPLICATES) {
//this function assumes the array dimensions are identical
	for (unsigned short int set= 0; set < 2; set++) {
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set]; r++) {
/*			if (dmr1[set][r].methylated_C > dmr1[set][r].coverage) {
				printf("@line %u, methylated_C = %zu > coverage = %zu, which isn't possible.\n",__LINE__,dmr1[set][r].methylated_C,dmr1[set][r].coverage);
				exit(EXIT_FAILURE);
			}
			if (dmr2[set][r].methylated_C > dmr2[set][r].coverage) {
				printf("@line %u, methylated_C = %zu > coverage = %zu, which isn't possible.\n",__LINE__,dmr1[set][r].methylated_C,dmr1[set][r].coverage);
				exit(EXIT_FAILURE);
			}*/
			dmr1[set][r].methylated_C += dmr2[set][r].methylated_C;
			dmr1[set][r].coverage += dmr2[set][r].coverage;
		}
	}
}


void subtract_2nd_dmr_from_1st_dmr(DMR_struct **restrict dmr1, DMR_struct **restrict dmr2, const unsigned short int *restrict NO_OF_REPLICATES) {
//this function assumes the array dimensions are identical
	for (unsigned short int set = 0; set < 2; set++) {
		for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
//			if (dmr1[set][replicate].methylated_C < dmr2[set][replicate].methylated_C) {
//				printf("%zu - %zu	Вниманиe! subtracting a large unsigned integer from a small unsigned integer will give an unsigned integer a negative value @line %u. :,(\n", dmr1[set][replicate].methylated_C, dmr2[set][replicate].methylated_C, LINE);
//				exit(EXIT_FAILURE);
//			}
//			if (dmr1[set][replicate].methylated_C > dmr1[set][replicate].coverage) {
//				printf("@line %u, methylated_C = %zu > coverage = %zu, which isn't possible.\n",__LINE__,dmr1[set][replicate].methylated_C,dmr1[set][replicate].coverage);
//				exit(EXIT_FAILURE);
//			}
			dmr1[set][replicate].methylated_C -= dmr2[set][replicate].methylated_C;
			dmr1[set][replicate].coverage -= dmr2[set][replicate].coverage;
		}
	}
}

void get_dmr_from_sum(DMR_struct ***restrict dmrs, DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, const size_t INDEX) {
//this function assumes the array dimensions are identical
	for (unsigned short int set= 0; set < 2; set++) {
//		printf("@ %s line %u doing set %u from 0 to < %u\n", __FILE__, __LINE__, set, NO_OF_REPLICATES[set]);
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set]; r++) {
			dmr[set][r].methylated_C = dmrs[set][r][INDEX].methylated_C;
			dmr[set][r].coverage     = dmrs[set][r][INDEX].coverage;
//			printf("dmr[%u][%u]5mC/Cyt = %zu/%zu\n", set, r, dmr[set][r].methylated_C, dmr[set][r].coverage);
		}
	}
}

void print_DMRs (const char *restrict FILENAME,
	DMR_struct ***restrict all_DMRs,
	const unsigned short SET1,
	const unsigned int *restrict DMR_chromosome_indices,
	CHROMOSOME_CHANGES_STRUCT ***restrict bs_data_chromosome_changes,
	const unsigned int *restrict DMR_starts,
	const unsigned int *restrict DMR_ends,
	const unsigned int NO_OF_DMRs,
	const unsigned short int *restrict NO_OF_REPLICATES,
	const char *restrict FDR_TYPE,
	const char *restrict label_set1,
	const char *restrict label_set2,
	const bool gene_annotation_defined,
	const unsigned PROMOTER_CUTOFF,
	const unsigned SKIPS,
	const bool printing_CpN_list,
	const unsigned *restrict nCpN,
	const unsigned *restrict diffCpN,
	const unsigned *restrict CpN_skips_in_this_DMR,
	char **restrict GENE_LIST,
	const bool PAIRED
	) {
	if (NO_OF_DMRs == 0) {
		return;
	}
	double *restrict pvalues = malloc(sizeof(double) * NO_OF_DMRs);
	if (pvalues == NULL) {
		fprintf(stderr, "Failed to malloc @ %s line %u\n", __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	DMR_struct **restrict single_DMR = malloc(sizeof(DMR_struct) * (unsigned)(NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1]));
	if (single_DMR == NULL) {
		fprintf(stderr, "Failed to malloc @ %s line %u\n", __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	for (unsigned short int set = 0; set < 2; set++) {
		single_DMR[set] = malloc(sizeof(DMR_struct) * NO_OF_REPLICATES[set]);
		if (single_DMR[set] == NULL) {
			fprintf(stderr, "Failed to malloc single_DMR[%u] @ %s line %u\n", set, __FILE__, __LINE__);
			perror("");
			exit(EXIT_FAILURE);
		}
	}
	for (unsigned int dmr = 0; dmr < NO_OF_DMRs; dmr++) {
		for (unsigned short set = 0; set < 2; set++) {
			for (unsigned short r = 0; r < NO_OF_REPLICATES[set]; r++) {
//				printf("@ %s line %u all_DMRs[%u][%u][%u]: methC = %zu, coverage = %zu\n", __FILE__, __LINE__, set, r, dmr, all_DMRs[set][r][dmr].methylated_C, all_DMRs[set][r][dmr].coverage);
				if (all_DMRs[set][r][dmr].coverage == 0) {
					printf("got coverage = 0\n");
					exit(EXIT_FAILURE);
				}
			}
		}/**/
		get_dmr_from_sum(all_DMRs,
		single_DMR,
		NO_OF_REPLICATES,
		dmr);
//void get_dmr_from_sum(DMR_struct ***restrict dmrs, DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, const size_t INDEX) {
		const double PVALUE = Pvalue(single_DMR, 
		NO_OF_REPLICATES,
		PAIRED);
		pvalues[dmr] = PVALUE;
	}
	double *restrict qvalues = p_adjust(pvalues, NO_OF_DMRs, FDR_TYPE);
/*
make filename for output
*/
	char *restrict filename = malloc((strlen(FILENAME)+strlen("_q."))+sizeof(char*));
	strcpy(filename, FILENAME);
	filename[strlen(filename)-4] = '\0';
	strcat(filename, "_q.tsv");
	FILE *restrict fh = fopen(filename, "w");
	if (fh == NULL) {
		fprintf(stderr, "failed to open %s for writing @ %s line %u\n", filename, __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	printf("Writing answers to %s\n", filename);
	free(filename); filename = NULL;
/*
Start Output Header
*/
	fprintf(fh,"Chromosome\tStart\tEnd\t#mCpN\t#Diff.CpN");
	fprintf(fh,"\tMean_Difference\t%s\t%s",label_set1, label_set2);
	if (SKIPS > 0) {
		fprintf(fh, "\tSkips");
	}
	fprintf(fh, "\tp(DMR)");
	if (FDR_TYPE == NULL) {
		fprintf(fh, "\tq(DMR, BH)");
	} else {
		fprintf(fh, "\tq(DMR, %s)", FDR_TYPE);
	}
	if (gene_annotation_defined == true) {//part of the answer file header
		fprintf(fh,"\tInside_Genes\tBetween_Genes\tGene_Promoter_Cutoff_%u_Nucleotides", PROMOTER_CUTOFF);
	}
	(void)printing_CpN_list;//the CpN list is only in the main answers file, so it doesn't get a header here
	fprintf(fh, "\n");
/*
End Output Header
*/
	for (unsigned int dmr = 0; dmr < NO_OF_DMRs; dmr++) {
//		printf("%s line %u: DMR_chromosome_indices[%u] = %u\n", __FILE__, __LINE__, dmr, DMR_chromosome_indices[dmr]);
		fprintf(fh, "%s", bs_data_chromosome_changes[SET1][0][DMR_chromosome_indices[dmr]].chromosome);
		fprintf(fh, "\t%u\t%u", DMR_starts[dmr], DMR_ends[dmr]);
		fprintf(fh, "\t%u\t%u", nCpN[dmr], diffCpN[dmr]);
		get_dmr_from_sum(all_DMRs, single_DMR, NO_OF_REPLICATES, dmr);
		fprintf(fh, "\t%.1f", 100.0*mean_percent_difference(single_DMR, NO_OF_REPLICATES));
//void print_methylation_percentages (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, FILE *restrict answers_filehandle) {
		print_methylation_percentages(single_DMR, NO_OF_REPLICATES, fh);
		if (SKIPS > 0) {
			fprintf(fh,"\t%u", CpN_skips_in_this_DMR[dmr]);
		}
		fprintf(fh, "\t%g\t%g", pvalues[dmr], qvalues[dmr]);
		if (gene_annotation_defined == true) {
			fprintf(fh, "\t%s\t%s\t%s", GENE_LIST[dmr*3], GENE_LIST[dmr*3+1], GENE_LIST[dmr*3+2]);
		}
		fprintf(fh, "\n");
	}
	for (unsigned short int set = 0; set < 2; set++) {
		free(single_DMR[set]); single_DMR[set] = NULL;
	}
	free(single_DMR); single_DMR = NULL;
	free(pvalues); pvalues = NULL;
	free(qvalues); qvalues = NULL;
	fclose(fh);
}

void print_methylation_heatmap (DMR_struct **restrict dmr, const unsigned short int *restrict NO_OF_REPLICATES, FILE *restrict tsv_filehandle, const size_t heatmap_row) {
//this assumes coverages are all > 0
	size_t mC_count = 0;
	size_t coverage_count = 0;
	for (unsigned short int set_loop = 0; set_loop < 2; set_loop++) {
//		fprintf(answers_filehandle,"\t[%.1f",100.0*(double)dmr[set_loop][0].methylated_C / (double)dmr[set_loop][0].coverage);
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set_loop]; r++) {
			mC_count += dmr[set_loop][r].methylated_C;
			coverage_count += dmr[set_loop][r].coverage;
//			fprintf(answers_filehandle,",%.1f",100.0*(double)dmr[set_loop][r].methylated_C / (double)dmr[set_loop][r].coverage);
		}
//		fprintf(answers_filehandle,"]");
	}
	for (unsigned short int set_loop = 0; set_loop < 2; set_loop++) {
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set_loop]; r++) {
			if (dmr[set_loop][r].coverage == 0) {
				return;//0/0 would write "nan" into the heatmap's table
			}
		}
	}
	const double MEAN = (double) mC_count / (double) coverage_count;
	double stdev = 0.0;
	for (unsigned short int set_loop = 0; set_loop < 2; set_loop++) {
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set_loop]; r++) {
			stdev += pow((MEAN - (double)(dmr[set_loop][r].methylated_C) / (double)dmr[set_loop][r].coverage), 2.0);
		}
	}
	if (stdev == 0.0) {
		return;//prevents division by 0
	}
	stdev = sqrt(stdev / (NO_OF_REPLICATES[0] + NO_OF_REPLICATES[1] - 1) );
	fprintf(tsv_filehandle, "%zu", heatmap_row);
	for (unsigned short int set_loop = 0; set_loop < 2; set_loop++) {
		for (unsigned short int r = 0; r < NO_OF_REPLICATES[set_loop]; r++) {
			const double x = (double)(dmr[set_loop][r].methylated_C) / (double)dmr[set_loop][r].coverage;
			fprintf(tsv_filehandle, "\t%lf", (x - MEAN) / stdev);
		}
	}
	fprintf(tsv_filehandle, "\n");
}

void make_heatmap_eps (
const char *restrict HEATMAP_TSV,
const char *restrict SET1_LABEL,
const char *restrict SET2_LABEL,
const char *restrict OPTIONS_STRING
) {
	const size_t LENGTH = strlen("heatmap_") 
	+ strlen(SET1_LABEL)
	+ strlen("_vs_")
	+ strlen(SET2_LABEL)
	+ strlen("_")
	+ strlen(".eps")
	+ strlen(OPTIONS_STRING)
	+ 1;//'\0'; the "_" and the '\0' used to be missing, so the .eps filename overflowed by 2 bytes
	char *restrict heatmap_r_filename = malloc(LENGTH);
	if (heatmap_r_filename == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	safe_string_copy(heatmap_r_filename, "heatmap_");
	strcat(heatmap_r_filename, SET1_LABEL);
	strcat(heatmap_r_filename, "_vs_");
	strcat(heatmap_r_filename, SET2_LABEL);
	strcat(heatmap_r_filename, "_");
	strcat(heatmap_r_filename, OPTIONS_STRING);
	strcat(heatmap_r_filename, ".R");
	
	FILE *restrict fh_r = fopen(heatmap_r_filename, "w");
	if (fh_r == NULL) {
		fprintf(stderr, "Failed to write %s at %s line %u\n", heatmap_r_filename, __FILE__, __LINE__);
		perror("");
		exit(EXIT_FAILURE);
	}
	char *restrict heatmap_eps_filename = malloc(LENGTH);
	if (heatmap_eps_filename == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	safe_string_copy(heatmap_eps_filename, "heatmap_");
	strcat(heatmap_eps_filename, SET1_LABEL);
	strcat(heatmap_eps_filename, "_vs_");
	strcat(heatmap_eps_filename, SET2_LABEL);
	strcat(heatmap_eps_filename, "_");
	strcat(heatmap_eps_filename, OPTIONS_STRING);
	strcat(heatmap_eps_filename, ".eps");
		
	fprintf(fh_r, "library(gplots)\n");
	fprintf(fh_r, "data <- read.table(file = '%s', sep='\\t', header=T)\n", HEATMAP_TSV);
	fprintf(fh_r, "row.names(data) <- data$Gene.Name\n");
	fprintf(fh_r, "data <- data[,2:dim(data)[2]]\n");
	fprintf(fh_r, "distancet <- as.dist(1-cor(t(as.matrix(data))))\n");
	fprintf(fh_r, "heirarchical_cluster <- hclust(distancet, method='ward.D2')\n");
	fprintf(fh_r, "dendcomplete <- as.dendrogram(heirarchical_cluster)\n");
	fprintf(fh_r, "setEPS()\n");
	fprintf(fh_r, "postscript('%s')\n", heatmap_eps_filename);
	free(heatmap_eps_filename); heatmap_eps_filename = NULL;
	fprintf(fh_r, "heatmap.2(as.matrix(data), col=greenred(75), scale='row', key=T, keysize=1, density.info='none', dendrogram=c('col'), trace='none', cexCol=1.2, margins=c(10,5),srtCol=45,Rowv=T,labRow=F, main = '%s vs. %s')\n", SET1_LABEL, SET2_LABEL);
	fprintf(fh_r, "dev.off()\n");
	fclose(fh_r);
	char *restrict command = malloc(LENGTH + strlen("Rscript ") * 2);
	if (command == NULL) {
		fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	safe_string_copy(command, "Rscript ");
	strcat(command, heatmap_r_filename);
	free(heatmap_r_filename); heatmap_r_filename = NULL;
	if (system(command) != 0) {
		fprintf(stderr, "%s FAILED.\n", command);
		perror("");
		free(command); command = NULL;
		exit(EXIT_FAILURE);
	}
	free(command); command = NULL;
}

/*unsigned int get_next_CpN (const unsigned int CURRENT_CPN, const size_t **restrict pmin, const size_t **restrict pmax, BISULFITE_DATA_STRUCT ***restrict fdata, const unsigned short int *restrict FNO_OF_REPLICATES, const unsigned short int *restrict FCOMPARED_SETS) {

	size_t **restrict fmin = malloc(sizeof(size_t)*(FNO_OF_REPLICATES[0]+FNO_OF_REPLICATES[1]));//this crap is done so that array values don't change
	size_t **restrict fmax = malloc(sizeof(size_t)*(FNO_OF_REPLICATES[0]+FNO_OF_REPLICATES[1]));
	for (unsigned short int fset = 0; fset < 2; fset++) {//do not alter the current indices, make a copy here
		fmin[fset] = malloc(sizeof(size_t)*(FNO_OF_REPLICATES[fset]));
		fmax[fset] = malloc(sizeof(size_t)*(FNO_OF_REPLICATES[fset]));
		for (unsigned short int fr = 0; fr < FNO_OF_REPLICATES[fset]; fr++) {
			fmin[fset][fr] = pmin[fset][fr];
			fmax[fset][fr] = pmax[fset][fr];
		}
	}

	unsigned int return_nucleotide = UINT_MAX;
	for (unsigned short int fset = 0; fset < 2; fset++) {//get all the CpN equal to each other, if possible
		for (unsigned short int fr = 0; fr < FNO_OF_REPLICATES[fset]; fr++) {
			if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide == CURRENT_CPN) {
				if (fmin[fset][fr] <= fmax[fset][fr]) {//prevent segfault, I'm not at the array max?
					if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]+1].nucleotide < return_nucleotide) {
//						printf("is %u the next CpG after %u?\n", fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide, CURRENT_CPN);
						return_nucleotide = fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]+1].nucleotide;
					}
				}
			} else if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide < CURRENT_CPN) {//if the 
				while (true) {
					fmin[fset][fr]++;
					if ((fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide > CURRENT_CPN) && (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide < return_nucleotide)) {
//						printf("is %u the next CpG after %u\n", fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide, CURRENT_CPN);
						return_nucleotide = fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide;
						break;
					}
					if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide > CURRENT_CPN) {
						break;
					}
					if (fmin[fset][fr] == fmax[fset][fr]) {
						break;
					}
				}
			} else if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide > CURRENT_CPN) {//keep going down until I go below
				while (true) {
					if (fmin[fset][fr] == 0) {
						break;
					}
					fmin[fset][fr]--;
					if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide <= CURRENT_CPN) {
						break;
					}
					if (fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide < return_nucleotide) {
//						printf("is %u the next CpG after %u\n", fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide, CURRENT_CPN);
						return_nucleotide = fdata[FCOMPARED_SETS[fset]][fr][fmin[fset][fr]].nucleotide;
					}
				}
			}
		}
	}

	if (return_nucleotide == UINT_MAX) {
		printf("Failed to find the next nucleotide @ %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	for (unsigned short int fset = 0; fset < 2; fset++) {//get all the CpN equal to each other, if possible
		free(fmin[fset]); fmin[fset] = NULL;
		free(fmax[fset]); fmax[fset] = NULL;
	}
	free(fmin); fmin = NULL;	free(fmax); fmax = NULL;
	return return_nucleotide;
}*/
#endif
