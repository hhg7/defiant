
    ____  ___________________    _   ________
   / __ \/ ____/ ____/  _/   |  / | / /_  __/
  / / / / __/ / /_   / // /| | /  |/ / / /   
 / /_/ / /___/ __/ _/ // ___ |/ /|  / / /    
/_____/_____/_/   /___/_/  |_/_/ |_/ /_/

Differential methylation: Easy, Fast, Identification and ANnoTation

 by David E. Condon, University of Pennsylvania, 2015-2017.


Help for DEFIANT
[1mSynopsis[0m
./defiant [OPTIONS]... [FILES]...

where [FILES] is a list of files with spaces separating groups and commas separating replicates within groups. Options are like "-a refFlat.gtf" etc. (cf. Options).

Minimal example: ./defiant -i control1,control2 case1,case2

All files in should be specified after "-i", where commas separate replicates in the same group.

[1mOptions[0m

[1m-a[0m	Specify annotation file, e.g. "-a mm10.gtf"

[1m-b[0m	Output DMRs in bed file. This option does not take an argument.

[1m-c[0m	minimum coverage, e.g. "-c 10".  This option accepts positive integers and can be parallelized to test multiple options.

[1m-CpN[0m	minimum CpN/CpG/CH/CHH in a DMR, e.g. "-CpN 10".  This option accepts positive integers and can parallelized.  "CpN" is case insensitive.

[1m-cpu[0m	Set number of CPU when running multiple options, e.g. "-cpu 4". "CPU" is case insensitive and accepts integers > 0.

[1m-d[0m	Minimum differential nucleotide count in a DMR, e.g. "-d 3".  This option can be parallelized.

[1m-D[0m	Maximum non-default options in a parallel run, e.g. "-D 4"

[1m-debug[0m	Turn on debugging mode.  This slows down the execution significantly, but can help diagnose problems if they arise.  This option does not accept any arguments.

[1m-E[0m	print statistics for every CpN. This option does not take an argument. This slows Defiant down significantly.

[1m-f[0m	make EPS figures for each DMR. Warning: requires R installation. This option does not take an argument, and will slow defiant's execution.

[1m-G[0m	Maximum allowed gap between CpN, e.g. "-G 1000"

[1m-h[0m	Print this help menu; "-h" is case insensitive. This option does not take an argument, and defiant exits after this option is read.

[1m-i[0m	Start reading input files.  This is the only required argument.  All further entries to the command line are assumed to be files.

[1m-I[0m	Print isolation of nucleotides in output files. This option does not take an argument.

[1m-l[0m	Set output file(s) label, e.g. "-l new"

[1m-L[0m	give labels for each set in a comma-delimited string, e.g. "-L case,control"

[1m-N[0m	list CpG Nucleotides in the DMR in output file. This option does not take an argument.

[1m-o[0m	overwrite existing files if present.  This option does not take an argument.

[1m-p[0m	Maximum p-value, which is 0<=p<=1.  This option can be parallelized to test multiple options.  Default 0.05.

[1m-P[0m	Minimum Percent methylation difference (0 <= P <= 100). This option can be parallelized to test multiple options (default 10%).

[1m-q[0m	Promoter cutoff for gene assignment of intergenic DMRs (default 10,000 nucleotides).  This option accepts positive integers, e.g. "-q 15000".

[1m-r[0m	Minimum nucleotide range, which accepts a non-negative integer.  Default range is 0 nucleotides.

[1m-R[0m	include "Random" chromosomes.  This option does not accept an argument.

[1m-s[0m	Maximum allowed consecutive similar CpN, default is 5 CpN.  This accepts non-negative integers, e.g. "-s 3".

[1m-S[0m	Allow some number of consecutive skips of low coverage, default is 0.  This accepts positive integers, e.g. "-S 1".

[1m-U[0m	Include "Un" chromosomes (default is to ignore them).  This option does not accept an argument.

[1m-x[0m	x-axis & legend labels in figures.  "-x" activates "-f" option and requires an R installation. This doesn't accept an argument.

[1mParallelization[0m

As each experiment is different, a different set of parameters may be appropriate for each experiment. You may not know these parameters ahead of time. Thus, defiant has been set to easily test multiple parameters in parallel via a shared-memory model. Parameters underlined and bold faced in Table 1, e.g. p, can be written like a C-style for loop and delimited with commas: ./defiant -p < min >, < max >, < step > which would increment the p-value from min to max in steps of step.
For example, 

./defiant -c 5,15,5

which would run minimum coverage from 5 to 15 in steps of 5.
If you only wish to run two parameters, you can simply write a comma in between the two parameters you wish to vary, e.g. 

./defiant -p 0.01,0.05

will run p = 0.01 and then p = 0.05.

The data is read off of the hard drive and into memory, which will then be shared among all the CPU. This is done to make 3D graphs, i.e. x vs. y with the 3rd dimension in color. However, the -D option can be used to vary all parameters as a nested for loop. I strongly recommend not to use the "-f" or "-x" options with multiple runs. One of defiant’s advantages is speed and low resource use, using both "-f" and "-x" options will make the runs take much much longer and potentially create a lot of files which will make I/O on your computer very slow. All DMR counts are then saved to a table, which will end in something like dmr_count.tsv

[1mInput Formats[0m
Defiant is set up to automatically identify and read the following input formats:


[1mInput Type 1[0m
[1mExample[0m:	chr1	762	763	0.1764	37
Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer in [0,4294967295].
Column3:	ignored.
Column4:	methylation percent, a floating point in [0,1].
Column5:	coverage, an unsigned integer, an unsigned integer in [0,4294967295]

[1mInput Type 2[0m	known for MethylKit input
[1mExample[0m:	chr1.762	chr1	762	R	10000	17.64	82.36
Column1:	unique name, this is ignored.
Column2:	chromosome, which is a string.
Column3:	nucleotide, an unsigned integer [0,4294967295].
Column4:	sense, this is ignored.
Column5:	coverage, an unsigned integer in [0,4294967295].
Column6:	methylation percent, a floating point in [0,100].
Column7:	cytosine percent, a floating point in [0,100].

[1mInput Type 3[0m
[1mExample[0m:	chr1	762	763	0.1764
Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer [0,4294967295].
Column3:	ignored.
Column4:	methylation percent, a floating point in [0,1].

[1mInput Type 4[0m
[1mExample[0m:	chr1	762	6	14
Column1:	chromosome, which is a string.
Column2:	nucleotide, an unsigned integer in [0,4294967295].
Column3:	methylated C count, an unsigned integer in [0,4294967295].
Column4:	C count, an unsigned integer [0,4294967295].

[1mInput Type 5[0m	Bismark coverage2cytosine format:
<chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>//Bismark coverage2cytosine format
[1mExample[0m:	chr1	762	763	+	17	64	CG	CGA
Column1:	chromosome, which is a string.
Column2:	nucleotide/start position, an unsigned integer [0,4294967295].
Column3:	strand.
Column4:	methylated C count, an unsigned integer in [0,4294967295].
Column5:	C count, an unsigned integer in [0,4294967295].
Column6:	C-context, e.g. CG, CH, CHH.
Column7:	C-context, e.g. CGA, CGT, etc.

[1mInput Type 6[0m	Bismark coverage2cytosine format:
<chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
[1mExample[0m:	chr1	762	763	0.265625	17	76
Column1:	chromosome, which is a string.
Column2:	nucleotide/start position, an unsigned integer in [0,4294967295].
Column3:	nucleotide/end position, an unsigned integer in [0,4294967295].
Column4:	methylation percentage, which is calculated by Defiant.
Column5:	methylated C count, an unsigned integer in [0,4294967295].
Column6:	C count, an unsigned integer in [0,4294967295].

[1mInput Type 7[0m	HELP-Tag data.  This can have a header.
[1mExample[0m:	1 chr1	762	763	0.2656	0.1776
Column1:	ignored
Column2:	chromosome, a string.
Column3:	position, an unsigned integer in [0,4294967295].
Column4:	methylation percent: a floating point number in [0,1].
Column5:	Conf. ignored.

[1mInput Type 8[0m	(EPP)Epigenome Processing Pipeline
[1mExample[0m:	chr1	762	763	'17/76'	999	+
Column1:	chromosome, which is a string.
Column2:	start nucleotide, an unsigned integer in [0,4294967295].
Column3:	end nucleotide, an unsigned integer in [0,4294967295].
Column4:	methylation percent as a fraction, two unsigned integers.  Coverage is given as the denominator. Everything after this column is ignored.

[1mInput Type 9[0m	Bsmooth Input
[1mExample[0m:	X	762	+	CG	17	76
Column1:	chromosome, which is a string.
Column2:	Nucleotide, an unsigned integer in [0,4294967295].
Column3:	strand sense, ignored.
Column4:	context, ignored.
Column5:	methylated C count, an unsigned integer in [0,4294967295].
Column6:	C count, an unsigned integer in [0,4294967295].

[1mInput Type 10[0m	BisSNP (found in RnBeads)
[1mExample[0m:	X	762	763	17.76	82.24	10000	762	763	180,60,0
Column1:	chromosome, which is a string.Column2:	Nucleotide start, an unsigned integer in [0,4294967295].
Column3:	Nucleotide end, an unsigned integer in [0,4294967295].
Column4:	methylation value in [0:100].
Column5:	Coverage, an unsigned integer in [0,4294967295]. Everything after this column is ignored.

[1mOutput[0m

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

[1mRunning Defiant[0m

Simplest possible case:
1. ./defiant -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

Label each sample.
2. ./defiant -L control,case -l my_name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

3. Annotation:
./defiant -a refFlat.txt -L control,case -l my ̇name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

4. Generate figures with "CpG" key (requires an installation of the R programming language.)
./defiant -x CpG -a refFlat.txt  -L control,case -l my ̇name -i control1.txt,control2.txt case1.txt,case2.txt,case3.txt

