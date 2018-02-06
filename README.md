# defiant
**D**ifferential methylation: **E**asy, **F**ast, **I**dentification and **AN**no**T**ation

by David E. Condon and Kyoung-Jae Won, University of Pennsylvania, 2015-2017.  Email Dave at dec986@gmail.com with questions/complaints/suggestions.  This has been published in BMC Bioinformatics: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2037-1

# installation
Download the file defiant.zip and unzip the file like this on a Linux command line: `unzip defiant.zip`
which opens a file called install.sh  Do this on the command line: `./install.sh` which will compile an executable "defiant" for your system.
WARNING: If you're using Mac, use `install_mac.sh`. _I strongly advise against using Mac, however.  I have found Mac to be significantly less reliable and usable than Linux.  Also, the multi-processing capability is absent with the Mac version as GCC and omp.h cannot work under Mac._

# Synopsis
./defiant [OPTIONS]... [FILES]...

where [FILES] is a list of files with spaces separating groups and commas separating replicates within groups. Options are like "-a refFlat.gtf" etc. (cf. Options).

Minimal example: ./defiant -i control1,control2 case1,case2

All files in should be specified after "-i", where commas separate replicates in the same group.

# Options

-a	Specify annotation file, e.g. "-a mm10.gtf"

-b	Output DMRs in bed file. This option does not take an argument.

-c	minimum coverage, e.g. "-c 10".  This option accepts positive integers and can be parallelized to test multiple options.

-CpN	minimum CpN/CpG/CH/CHH in a DMR, e.g. "-CpN 10".  This option accepts positive integers and can parallelized.  "CpN" is case insensitive.

-cpu	Set number of CPU when running multiple options, e.g. "-cpu 4". "CPU" is case insensitive and accepts integers > 0.

-d	Minimum differential nucleotide count in a DMR, e.g. "-d 3".  This option can be parallelized.

-debug Turn on debugging mode.  This slows down the execution significantly, but can help diagnose problems if they arise.  This option does not accept any arguments.

-D	Maximum non-default options in a parallel run, e.g. "-D 4"

-debug	Turn on debugging mode.  This slows down the execution significantly, but can help diagnose problems if they arise.  This option does not accept any arguments.

-E	print statistics for every CpN. This option does not take an argument. This slows Defiant down significantly.

-f	make EPS figures for each DMR. Warning: requires R installation. This option does not take an argument, and will slow defiant's execution.

-fdr Calculate FDR-adjusted q-value for each CpN.  'FDR' is case insensitive.  This option can take case-insensitive arguments 'fdr' or 'bh' for Benjamini-Hothberg method, 'Bonferroni', 'Hochberg', 'Hommel', 'Holm', or 'BY' for Benjamini & Yekutieli.  If no argument is given, 'Holm' is assumed.  This function is a translation of R's 'p.adjust'.  I recommend against using this as for genome-scale CpG measurements, almost everything will be q = 1 and no DMRs will be obtained in any case.  This option will substantially increase RAM use and slow execution.  'Hommel' is so slow I strongly recommend against it.

-G	Maximum allowed gap between CpN, e.g. "-G 1000"

-h	Print this help menu; "-h" is case insensitive. This option does not take an argument, and defiant exits after this option is read.

-i	Start reading input files.  This is the only required argument.  All further entries to the command line are assumed to be files.

-l	Set output file(s) label, e.g. "-l new"

-L	give labels for each set in a comma-delimited string, e.g. "-L case,control"

-N	list CpG Nucleotides in the DMR in output file. This option does not take an argument.

-o	overwrite existing files if present.  This option does not take an argument.

-p	Maximum p-value, which is 0<=p<=1.  This option can be parallelized to test multiple options.  Default 0.05.

-P	Minimum Percent methylation difference (0 <= P <= 100). This option can be parallelized to test multiple options (default 10%).

-q	Promoter cutoff for gene assignment of intergenic DMRs (default 10,000 nucleotides).  This option accepts positive integers, e.g. "-q 15000".

-r	Minimum nucleotide range, which accepts a non-negative integer.  Default range is 0 nucleotides.

-R	include "Random" chromosomes.  This option does not accept an argument.

-s	Maximum allowed consecutive similar CpN, default is 5 CpN.  This accepts non-negative integers, e.g. "-s 3".

-S	Allow some number of consecutive skips of low coverage, default is 0.  This accepts positive integers, e.g. "-S 1".

-U	Include "Un" chromosomes (default is to ignore them).  This option does not accept an argument.

-v	Print a p-value for each DMR. This option accepts the same arguments that the '-FDR' option does.

-x	x-axis & legend labels in figures.  "-x" activates "-f" option and requires an R installation. This doesn't accept an argument.

# Parallelization

As each experiment is different, a different set of parameters may be appropriate for each experiment. You may not know these parameters ahead of time. Thus, defiant has been set to easily test multiple parameters in parallel via a shared-memory model. Parameters underlined and bold faced in Table 1, e.g. p, can be written like a C-style for loop and delimited with commas: ./defiant -p < min >, < max >, < step > which would increment the p-value from min to max in steps of step.
For example, 

./defiant -c 5,15,5

which would run minimum coverage from 5 to 15 in steps of 5.
If you only wish to run two parameters, you can simply write a comma in between the two parameters you wish to vary, e.g. 

./defiant -p 0.01,0.05

will run p = 0.01 and then p = 0.05.

The data is read off of the hard drive and into memory, which will then be shared among all the CPU. This is done to make 3D graphs, i.e. x vs. y with the 3rd dimension in color. However, the -D option can be used to vary all parameters as a nested for loop. I strongly recommend not to use the "-f" or "-x" options with multiple runs. One of defiant’s advantages is speed and low resource use, using both "-f" and "-x" options will make the runs take much much longer and potentially create a lot of files which will make I/O on your computer very slow. All DMR counts are then saved to a table, which will end in something like dmr_count.tsv  This file can be processed with plot_results.pl which will produce a tidy report on how each parameter influences the number of DMRs found.  The perl script plot_results.pl depends on LaTeX and GNUPlot.

# Input Formats
Input data saves each CpN as a chromosome, a nucleotide number, a count unmethylated, and a count methylated.
Each CpN is shown as a single line in all input formats.
The number 4294967295 appears frequently here, and is the maximum unsigned integer value on a 64-bit machine.
Defiant is set up to automatically identify and read the following input formats:

## Input Type 1
Example:	chr1	762	763	0.1764	37

Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer in [0,4294967295].
Column3:	ignored.
Column4:	methylation percent, a floating point in [0,1].
Column5:	coverage, an unsigned integer, an unsigned integer in [0,4294967295]

## Input Type 2	known for MethylKit input
Example:	chr1.762	chr1	762	R	10000	17.64	82.36

Column1:	unique name, this is ignored.
Column2:	chromosome, which is a string.
Column3:	nucleotide, an unsigned integer [0,4294967295].
Column4:	sense, this is ignored.
Column5:	coverage, an unsigned integer in [0,4294967295].
Column6:	methylation percent, a floating point in [0,100].
Column7:	cytosine percent, a floating point in [0,100].

## Input Type 3
Example:	chr1	762	763	0.1764

Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer [0,4294967295].
Column3:	ignored.
Column4:	methylation percent, a floating point in [0,1].

## Input Type 4
Example:	chr1	762	6	14

Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer in [0,4294967295].
Column3:	methylated C count, an unsigned integer in [0,4294967295].
Column4:	C count, an unsigned integer [0,4294967295].

## Input Type 5	Bismark coverage2cytosine format:
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>//Bismark coverage2cytosine format
Example:	chr1	762	763	+	17	64	CG	CGA

Column1:	chromosome, which is a string.
Column2:	nucleotide/start position, an unsigned integer [0,4294967295].
Column3:	strand.
Column4:	methylated C count, an unsigned integer in [0,4294967295].
Column5:	C count, an unsigned integer in [0,4294967295].
Column6:	C-context, e.g. CG, CH, CHH.
Column7:	C-context, e.g. CGA, CGT, etc.

## Input Type 6	Bismark coverage2cytosine format:
<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
Example:	chr1	762	763	0.265625	17	76

Column1:	chromosome, which is a string.
Column2:	nucleotide/start position, an unsigned integer in [0,4294967295].
Column3:	nucleotide/end position, an unsigned integer in [0,4294967295].
Column4:	methylation percentage, which is calculated by Defiant.
Column5:	methylated C count, an unsigned integer in [0,4294967295].
Column6:	C count, an unsigned integer in [0,4294967295].

## Input Type 7	HELP-Tag data.  This can have a header.
Example:	1 chr1	762	763	0.2656	0.1776

Column1:	ignored
Column2:	chromosome, a string.
Column3:	position, an unsigned integer in [0,4294967295].
Column4:	methylation percent: a floating point number in [0,1].
Column5:	Conf. ignored.

## Input Type 8	(EPP)Epigenome Processing Pipeline
Example:	chr1	762	763	'17/76'	999	+

Column1:	chromosome, which is a string.
Column2:	start nucleotide, an unsigned integer in [0,4294967295].
Column3:	end nucleotide, an unsigned integer in [0,4294967295].
Column4:	methylation percent as a fraction, two unsigned integers.  Coverage is given as the denominator. Everything after this column is ignored.

## Input Type 9	Bsmooth Input
Example:	X	762	+	CG	17	76

Column1:	chromosome, which is a string.
Column2:	Nucleotide, an unsigned integer in [0,4294967295].
Column3:	strand sense, ignored.
Column4:	context, ignored.
Column5:	methylated C count, an unsigned integer in [0,4294967295].
Column6:	C count, an unsigned integer in [0,4294967295].

## Input Type 10	BisSNP (found in RnBeads)
Example:	X	762	763	17.76	82.24	10000	762	763	180,60,0

Column1:	chromosome, which is a string.Column2:	Nucleotide start, an unsigned integer in [0,4294967295].
Column3:	Nucleotide end, an unsigned integer in [0,4294967295].
Column4:	methylation value in [0:100].
Column5:	Coverage, an unsigned integer in [0,4294967295]. Everything after this column is ignored.

# Output

The file names are formatted according to the options set at the command line. For example, consider the output file:

control_vs_case_c10_CpN5_d1_p0.01_P10.tsv

• the start of the filename indicates that the "control" group was compared against the "case" group, as indicated by the "-L" option on the command line
• minimum coverage of 10,
• minimum CpN count = 5,
• minimum differentially methylated CpN = 1,
• a maximum p-value of 0.01
• a minimum Percent difference of 10%.
• If the Gap, "-G" option, is altered, -G<gap> will be in the output. If there are any allowed consecutive skips, this will be
indicated in the filename by " S2" if there are 2 allowed skips, for example.

# Running Defiant

## Simplest possible case:

1. ./defiant -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

## Label each sample.
2. ./defiant -L control,case -l my_name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

## Annotation:
./defiant -a refFlat.txt -L control,case -l my ̇name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

## Generate figures with "CpG" key (requires an installation of the R programming language.)
./defiant -x CpG -a refFlat.txt  -L control,case -l my ̇name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt
