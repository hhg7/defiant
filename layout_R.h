#ifndef DEFIANT_LAYOUT_R_H
#define DEFIANT_LAYOUT_R_H
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//requires dmr.h (DMR_struct) to be included first

void R_number (FILE *restrict rscript, const char *restrict FORMAT, const double NUMBER) {
//C prints 0/0 as "nan" or "-nan", which R can't read; a replicate without coverage is NA in R
	if (isfinite(NUMBER)) {
		fprintf(rscript, FORMAT, NUMBER);
	} else {
		fprintf(rscript, "NA");
	}
}

//NO_OF_CPN must have similar CpN subtracted
void R_figure (FILE *restrict rscript, const unsigned int *restrict CpN_list, const unsigned int NO_OF_CPN, const char *restrict LEGEND1, const char *restrict LEGEND2, const char *restrict XAXIS_LABEL, const char *restrict LIST_GENE, const double *restrict P, const double *restrict PERCENT, const double *restrict SUM_P, const double *restrict SUM_PERCENT, DMR_struct ***restrict array, const unsigned short int *restrict NO_OF_REPLICATES, const unsigned int START, const unsigned int END) {
	if (NO_OF_CPN < 1) {
		printf("there must be at least one CpN to print into the image. This is at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
//clear everything
	fprintf(rscript, "g1 <- c()\ng2 <- c()\np <- c()\npercent <- c()\nsum_p <- c()\nsum_percent <- c()\ng1_means <- c()\ng2_means <- c()\nxaxis <- c()\n");
	for (unsigned int cpn = 0; cpn < NO_OF_CPN; cpn++) {
		if ((CpN_list[cpn] < START) || (CpN_list[cpn] > END)) {
			printf("DMR = [%u,%u] yet you put in %u\n", START, END, CpN_list[cpn]);
			exit(EXIT_FAILURE);
		}
		if (P[cpn] < 0) {
			printf("somehow got p = %lf @ %s %u\n", P[cpn], __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		if (SUM_P[cpn] < 0) {
			printf("somehow got sum_p = %lf @ %s %u\n", SUM_P[cpn], __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
		fprintf(rscript, "sum_p[%u] <- ", cpn+1);
		R_number(rscript, "%g", SUM_P[cpn]);
		fprintf(rscript, "\nxaxis[%u] <- %u\n", cpn+1, CpN_list[cpn]);
		fprintf(rscript, "p[%u] <- ", cpn+1);
		R_number(rscript, "%g", P[cpn]);
		fprintf(rscript, "\npercent[%u] <- ", cpn+1);
		R_number(rscript, "%g", PERCENT[cpn]);
		fprintf(rscript, "\nsum_percent[%u] <- ", cpn+1);
		R_number(rscript, "%g", SUM_PERCENT[cpn]);
		fprintf(rscript, "\n");
		for (unsigned short int set = 0; set < 2; set++) {
			double *restrict replicates = malloc(sizeof(double)*NO_OF_REPLICATES[set]);
			if (replicates == NULL) {
				fprintf(stderr, "malloc failed @ %s line %u\n", __FILE__, __LINE__);
				exit(EXIT_FAILURE);
			}
			size_t set_methylated_C = 0, set_coverage = 0;
			for (unsigned short int replicate = 0; replicate < NO_OF_REPLICATES[set]; replicate++) {
				replicates[replicate] = 100.0*(double)array[set][replicate][cpn].methylated_C/(double)array[set][replicate][cpn].coverage;
				set_methylated_C += array[set][replicate][cpn].methylated_C;
				set_coverage += array[set][replicate][cpn].coverage;
			}
			fprintf(rscript, "g%u[%u] <- list(c(", set+1, cpn+1);//write the list of methylation percentages for the boxplot
			R_number(rscript, "%lf", replicates[0]);
			for (unsigned short int replicate = 1; replicate < NO_OF_REPLICATES[set]; replicate++) {
				fprintf(rscript, ",");
				R_number(rscript, "%lf", replicates[replicate]);
			}
			free(replicates); replicates = NULL;
			fprintf(rscript, "))\n");
			fprintf(rscript, "g%u_means[%u] <- ", set+1, cpn+1);
			R_number(rscript, "%lf", 100.0*(double)set_methylated_C / (double)set_coverage);//write the weighted mean
			fprintf(rscript, "\n");
		}
	}
	const bool TITLE_PRESENT = strcmp(LIST_GENE,"n/a") == 0 ? false : true;
	if (TITLE_PRESENT == false) {//this image will have a title
		fprintf(rscript, "layout(matrix(c(1,1,2,2,3,3,4,4), nrow = 4, byrow = TRUE), heights = c(1,1,1,1.46))\n");
	} else {//this image will not have a title
		fprintf(rscript, "layout(matrix(c(0,0,1,1,2,2,3,3,4,4), nrow = 5, byrow = TRUE), heights = c(0.2,1,1,1,1.46))\n");
	}
	fprintf(rscript, "par(mar = c(bottom, left, top, right), las=2)\n");/*las=2 sets y-axis labels horizontal style */
	fprintf(rscript, "boxplot(g1, xaxt = 'n', range = 0, ylab = 'Methylation %%', ylim = c(0,100), col = 'white', cex.lab=1.5, cex.axis=cex_main, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(g1)+0.5), yaxt = 'n')\n");
	fprintf(rscript, "axis(2,cex.axis=0.75)\n");
	if (TITLE_PRESENT == true) {
		fprintf(rscript, "title('%s', outer = TRUE, line = -1.5)\n", LIST_GENE);
	}
	fprintf(rscript, "lines(g1_means, col='dark green', lwd = 3)\n");
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");
	fprintf(rscript, "legend('topright',inset = c(legend_space,0), c('%s', 'Weighted Mean'), col = c('black','dark green'), lwd = c(1,3))\n", LEGEND1);
	fprintf(rscript, "boxplot(g2, xaxt = 'n', range = 0, main = NULL, ylim = c(0,100), ylab = 'Methylation %%', col = 'gray',  cex.lab=1.5, cex.axis=cex_main, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(g2)+0.5), yaxt = 'n')\n");
	fprintf(rscript, "axis(2,cex.axis=0.75)\n");//resize y tick label to prevent overlap
	fprintf(rscript, "lines(g2_means, col='dark green', lwd = 3)\n");
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");
	fprintf(rscript, "legend('topright',inset = c(legend_space,0),c('%s', 'Weighted Mean'), col = c('black','dark green'), lwd = c(1,3))\n", LEGEND2);
	fprintf(rscript, "par(mar = c(bottom, left, top, right))#'mar’ A numerical vector of the form 'c(bottom, left, top, right)’\n");
	fprintf(rscript, "plot(p, xaxt='n', ylab = 'P', type = 'l', lty = 1, lwd = 3, cex.lab=1.5, cex.axis=1, cex.main=cex_main, cex.sub=1.5, , log = 'y', xlim = c(0.5,length(p)+0.5), ylim = c(min(p,sum_p), max(p, sum_p)), yaxt = 'n')\n");
	fprintf(rscript, "axis(2,cex.axis=0.75)\n");//resize y tick label to prevent overlap
	fprintf(rscript, "points(sum_p, xaxt='n', ylab = 'P', type = 'l', col = 'blue', lty = 2, lwd = 3)\n");
//	fprintf(rscript, "segments(%u, min(p)/100, %u, 6, col = 'red', lwd = 3, lty = 4)\n", DMR_START, DMR_START);
//	fprintf(rscript, "segments(%u, min(p)/100, %u, 6, col = 'red', lwd = 3, lty = 4)\n", DMR_END, DMR_END);
	fprintf(rscript, "legend('topright', inset=c(legend_space,0), c('%s P', 'Moving P Mean'), col = c('black','blue'), lwd=c(3,3), lty=c(1,2))\n", XAXIS_LABEL);//'DMR Borders removed
	fprintf(rscript, "par(mar = c(bottom+5, left, top, right))#'mar’ A numerical vector of the form 'c(bottom, left, top, right)’\n");
	fprintf(rscript, "plot(percent, xaxt='n', ylab = 'Methylation %% Diff.', xlab = '%s', type = 'l', lty = 1, lwd = 3, cex.lab=1.5, cex.axis=1, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(percent)+0.5), ylim = c(min(percent, sum_percent), max(percent, sum_percent)), yaxt = 'n')\n", XAXIS_LABEL);
	fprintf(rscript, "axis(2,cex.axis=0.75)\n");//resize y tick label to prevent overlap
	fprintf(rscript, "points(sum_percent, xaxt='n', xlab = '%s', type = 'l', col = 'blue', lty = 2, lwd = 3)\n", XAXIS_LABEL);
//	fprintf(rscript, "segments(%u, -100, %u, 100, col = 'red', lwd = 3, lty = 4)\n", DMR_START, DMR_START);
//	fprintf(rscript, "segments(%u, -100, %u, 100, col = 'red', lwd = 3, lty = 4)\n", DMR_END, DMR_END);
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");

	fprintf(rscript, "legend('topright', inset=c(legend_space,0), c('%s Percent', 'Moving Mean %%'), col = c('black','blue'), lwd=c(3,3), lty=c(1,2))\n", XAXIS_LABEL);
	fprintf(rscript, "axis(las=1, 1, at=1:length(xaxis), xaxis)\n");
	fprintf(rscript, "dev.off()\n");
	fprintf(rscript, "g1 <- c()\ng2 <- c()\np <- c()\npercent <- c()\nsum_p <- c()\nsum_percent <- c()\ng1_means <- c()\ng2_means <- c()\nxaxis <- c()\n");
}
/*
void R_layout (FILE *restrict rscript, const unsigned int *restrict CpN_list, const unsigned int NO_OF_CPN, const unsigned int SIMILAR_CPN, const char *restrict LEGEND1, const char *restrict LEGEND2, const char *restrict XAXIS_LABEL, const char *restrict LIST_GENE) {
	if (NO_OF_CPN < 1) {
		printf("there must be at least one CpN to print into the image. This is at %s line %u\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}
	const bool TITLE_PRESENT = strcmp(LIST_GENE,"n/a") == 0 ? false : true;
	if (SIMILAR_CPN > 0) {//remove the trails
		fprintf(rscript, "g1 <- g1[1:(length(g1)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "g2 <- g2[1:(length(g2)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "g1_means <- g1_means[1:(length(g1_means)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "g2_means <- g2_means[1:(length(g2_means)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "p <- p[1:(length(p)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "percent <- percent[1:(length(percent)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "sum_p <- sum_p[1:(length(sum_p)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "sum_percent <- sum_percent[1:(length(sum_percent)-%u)]\n", SIMILAR_CPN);
		fprintf(rscript, "xaxis <- xaxis[1:(length(xaxis)-%u)]\n", SIMILAR_CPN);
	}
	if (TITLE_PRESENT == false) {//this image will have a title
		fprintf(rscript, "layout(matrix(c(1,1,2,2,3,3,4,4), nrow = 4, byrow = TRUE), heights = c(1,1,1,1.4))\n");
	} else {//this image will not have a title
		fprintf(rscript, "layout(matrix(c(0,0,1,1,2,2,3,3,4,4), nrow = 5, byrow = TRUE), heights = c(0.2,1,1,1,1.4))\n");
	}
	fprintf(rscript, "par(mar = c(bottom, left, top, right))\n");
	fprintf(rscript, "boxplot(g1, xaxt = 'n', range = 0, ylab = 'Methylation %%', ylim = c(0,100), col = 'white', cex.lab=1.5, cex.axis=cex_main, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(g1)+0.5))\n");
	if (TITLE_PRESENT == true) {
		fprintf(rscript, "title('%s', outer = TRUE, line = -1.5)\n", LIST_GENE);
	}
	fprintf(rscript, "lines(g1_means, col='dark green', lwd = 3)\n");
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");
	fprintf(rscript, "legend('topright',inset = c(legend_space,0), c('%s', 'Weighted Mean'), col = c('black','dark green'), lwd = c(1,3))\n", LEGEND1);
	fprintf(rscript, "boxplot(g2, xaxt = 'n', range = 0, main = NULL, ylim = c(0,100), ylab = 'Methylation %%', col = 'gray',  cex.lab=1.5, cex.axis=cex_main, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(g2)+0.5))\n");
	fprintf(rscript, "lines(g2_means, col='dark green', lwd = 3)\n");
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");
	fprintf(rscript, "legend('topright',inset = c(legend_space,0),c('%s', 'Weighted Mean'), col = c('black','dark green'), lwd = c(1,3))\n", LEGEND2);
	fprintf(rscript, "par(mar = c(bottom, left, top, right))#'mar’ A numerical vector of the form 'c(bottom, left, top, right)’\n");
	fprintf(rscript, "plot(p, xaxt='n', ylab = 'P', type = 'l', lty = 1, lwd = 3, cex.lab=1.5, cex.axis=1, cex.main=cex_main, cex.sub=1.5, , log = 'y', xlim = c(0.5,length(p)+0.5), ylim = c(min(p,sum_p), max(p, sum_p)))\n");
	fprintf(rscript, "points(sum_p, xaxt='n', ylab = 'P', type = 'l', col = 'blue', lty = 2, lwd = 3)\n");
//	fprintf(rscript, "segments(%u, min(p)/100, %u, 6, col = 'red', lwd = 3, lty = 4)\n", DMR_START, DMR_START);
//	fprintf(rscript, "segments(%u, min(p)/100, %u, 6, col = 'red', lwd = 3, lty = 4)\n", DMR_END, DMR_END);
	fprintf(rscript, "legend('topright', inset=c(legend_space,0), c('%s P', 'Moving P Mean'), col = c('black','blue'), lwd=c(3,3), lty=c(1,2))\n", XAXIS_LABEL);//'DMR Borders removed
	fprintf(rscript, "par(mar = c(bottom+5, left, top, right))#'mar’ A numerical vector of the form 'c(bottom, left, top, right)’\n");
	fprintf(rscript, "plot(percent, xaxt='n', ylab = 'Methylation %% Diff.', xlab = '%s', type = 'l', lty = 1, lwd = 3, cex.lab=1.5, cex.axis=1, cex.main=cex_main, cex.sub=1.5, xlim = c(0.5,length(percent)+0.5), ylim = c(min(percent, sum_percent), max(percent, sum_percent)))\n", XAXIS_LABEL);
	fprintf(rscript, "points(sum_percent, xaxt='n', xlab = '%s', type = 'l', col = 'blue', lty = 2, lwd = 3)\n", XAXIS_LABEL);
//	fprintf(rscript, "segments(%u, -100, %u, 100, col = 'red', lwd = 3, lty = 4)\n", DMR_START, DMR_START);
//	fprintf(rscript, "segments(%u, -100, %u, 100, col = 'red', lwd = 3, lty = 4)\n", DMR_END, DMR_END);
	fprintf(rscript, "par(xpd=TRUE, las=2)\n");

	fprintf(rscript, "legend('topright', inset=c(legend_space,0), c('%s Percent', 'Moving Mean %%'), col = c('black','blue'), lwd=c(3,3), lty=c(1,2))\n", XAXIS_LABEL);
	fprintf(rscript, "axis( 1, at=1:length(xaxis), xaxis)\n");
	fprintf(rscript, "dev.off()\n");
	fprintf(rscript, "g1 <- c()\ng2 <- c()\np <- c()\npercent <- c()\nsum_p <- c()\nsum_percent <- c()\ng1_means <- c()\ng2_means <- c()\nxaxis <- c()\n");
}*/
#endif
