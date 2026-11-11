#!/usr/bin/env perl

use strict;
use warnings::unused;
use warnings FATAL => 'all';
use feature 'say';
use autodie ':all';
use Const::Fast;
use Text::Diff;
use Storable;
const my $annotation_file => 'head1000.Mus_musculus.GRCm38.94.gtf';

if (not -e $annotation_file) {
	say "$annotation_file doesn't exist.";
	die;
}
open my $out, '>', "$0.cmd";
sub execute {
	my $command = shift;
	say "Executing: $command";
	say $out $command;
	if (system($command) != 0) {
		die;
	}
	0
}

execute("gcc -O4 -o new new.c -Wall -pedantic -std=gnu11 -lm -Wextra -Wconversion");
const my $new => './new';

const my $dir => '/home/con/DNA_Methylation';
const my @file_types => (2,4..6,8,10,11);
foreach my $type ('RRBS', 'WGBS') {
	foreach my $x (1 ,2) {
		foreach my $y (1..4) {
			my $file_list = '';
			my @translated_files;
			foreach my $group ('B','F') {
				foreach my $replicate (1..10) {
					my $n;
					if ($replicate == 10) {
						$n = 10;
					} elsif ($replicate < 10) {
						$n = "0$replicate";
					} else {
						print "$replicate doesn't fit my pigenholes.\n";
					}
					my $file = $group . "L$n" . "_artif.chr10.$type" . "_background.$x" . "_DMR.class.$y.bsmooth";
					unless (-e $file) {
						print "$file doesn't exist.\n";
						die;
					}
					if ($replicate == 1) {
						$file_list .= $file;
					} else {
						$file_list .= ",$file";
					}
					my @fh;
					open my $original, '<', $file;
					foreach my $file_type (@file_types) {
						my $new_name = $group . "L$n" . "_artif.chr10.$type" . "_background.$x" . "_DMR.class.$y.file_type$file_type.bsmooth";
						say "now writing $new_name";
						push @translated_files, $new_name;
						open $fh[$file_type], '>', $new_name;
					}
					while (<$original>) {
						my @line = split;
						shift @line;#chromosome is always 10
						my $nuc   = shift @line;
						my $sense        = shift @line;# + or -
						shift @line;#CG
						my $methC        = shift @line;
						my $c_count      = shift @line;
						my $coverage     = $methC + $c_count;
						my $meth_percent = ($methC / $coverage);
						my $meth_percent100 = 100 * $meth_percent;
						my $c_percent    = ($c_count / $coverage);
						my $c_percent100 = 100 * $c_percent;
						my $fr           = $sense eq '-' ? 'R' : 'F';# file type 2 names the strands F & R
						my $nuc1         = $nuc+1;
# file type 2
						say {$fh[2]} "chr10.$nuc\tchr10\t$nuc	$fr	$coverage	$meth_percent100	$c_percent100";# both percents are in [0,100]; $c_percent was in [0,1]
#						say {$fh[3]} "chr10\t$nuc	$nuc1	$meth_percent";#3 has no coverage, and will come with different answers
						say {$fh[4]} "chr10\t$nuc	$methC	$c_count";
						say {$fh[5]} "chr10\t$nuc	$sense	$methC	$c_count	CG	CGA";# the real strand: 2 lines of the same CpG are only merged if they are + then -
						say {$fh[6]} "chr10\t$nuc	$nuc1	$meth_percent	$methC	$c_count";
#						say {$fh[7]} "chr10\t$nuc	$nuc1	$meth_percent";#has no coverage
						say {$fh[8]} "chr10\t$nuc	$nuc1	'$methC/$coverage'	9	+";
						say {$fh[10]} "10\t$nuc	$nuc1	$meth_percent100	$coverage\tX\tY\t1,0,0\t0\t0";
						say {$fh[11]} join ("\t", 'chr10', 'C', $nuc, 'CG', 'CG', $meth_percent, $methC, $coverage);
					}
				}
				$file_list .= ' ';
			}
			const my $defiant_output => "b_vs_f_defiant$type.$x.$y" . '_c10_CpN5_d1_p0.05_P10.tsv';
			if (not -e $defiant_output) {
				say "$defiant_output doesn't exist";
				die;
			}
			foreach my $file_type (@file_types) {
				if ($file_type == 2) {
					next;
				}
				my $file_type_list = $file_list;
				const my $label => "new.file_type$file_type.$type.$x.$y";
				$file_type_list =~ s/\.bsmooth/.file_type$file_type.bsmooth/g;
				execute("/usr/bin/time $new -E -L b,f -l $label -i $file_type_list");
				my $output = "b_vs_f_$label" . '_c10_CpN5_d1_p0.05_P10.tsv';
				my $diff = diff $output, $defiant_output;
				if ($diff =~ m/^$/) {
					say "$defiant_output & $output are the same:	$label succeeded.";
				} else {
#					unlink @translated_files;#they use too much disk space
					say "$defiant_output & $output are different:	$label failed.";
					say "diff $defiant_output $output";
#					die;
				}
			}
			my $label = "new.file_type9.$type.$x.$y";
			execute("/usr/bin/time $new -L b,f -l $label -i $file_list");
			my $new_output = "b_vs_f_$label" . '_c10_CpN5_d1_p0.05_P10.tsv';
			my $diff = diff $new_output, $defiant_output;
			if ($diff =~ m/^$/) {
				say "\n$label succeeded.\n";
			} else {
				say "$label failed.";
				die;
			}
			unlink @translated_files;#they use too much disk space
			foreach my $fdr ('bh','by','hochberg','holm', 'hommel') {
				my $fdr_label = "new$type.$x.$y.$fdr";
				execute("/usr/bin/time $new -v $fdr -L b,f -l $fdr_label -i $file_list");
				my $no_fdr_file = "b_vs_f_$fdr_label" . '_c10_CpN5_d1_p0.05_P10.tsv';
				my $fdrs_found = 1;
				open my $p, '<', $no_fdr_file;
				while (<$p>) {
					if (/^No_DMR_was_found_with_current_parameters/) {
						$fdrs_found = 0;
						last
					}
					if (/^[0-9XYM]{1,2}\s/) {# [1-9] never matched chromosome 10
						last # there are DMRs
					}
				}
				close $p;
				if ($fdrs_found == 0) {
					next #there won't be a '_q.tsv' file
				}
				my $fdr_new_output = "b_vs_f_$fdr_label" . '_c10_CpN5_d1_p0.05_P10_q.tsv';
				my $fdr_defiant_output = $fdr_new_output;
				$fdr_defiant_output =~ s/new/defiant/;
				my $x2 = diff $fdr_new_output, $fdr_defiant_output;
				if ($x2 =~ m/^$/) {
					say "\n$fdr_label succeeded.\n";
				} else {
					say "$fdr_label failed.";
					say $x2;
					die;
				}
			}
		}
	}
}
#test multiple groups
#const my $g1 => 'BL01_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL02_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL03_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL04_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL05_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL06_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL07_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL08_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL09_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,BL10_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth';
#const my $g2 => 'FL01_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL02_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL03_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL04_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL05_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL06_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL07_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL08_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL09_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth,FL10_artif.chr10.RRBS_background.1_DMR.class.1.bsmooth';
#const my $g3 => 'BL01_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL02_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL03_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL04_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL05_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL06_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL07_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL08_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL09_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth,BL10_artif.chr10.RRBS_background.1_DMR.class.2.bsmooth';


#test multiple genes in same DMR
#execute("/usr/bin/time ./new -l some_name -v -L s1,s2,s3 -i $g1 $g2 $g3");
#my $x = `diff control_vs_FeDef_newdup_c10_CpN5_d1_p0.05_P10.tsv control_vs_FeDef_new_c10_CpN5_d1_p0.05_P10.tsv`;
#exit;
#const my $xref => retrieve('dup_diff.pldata');

#if ($x ne $$xref) {
# 	say "test of multiple genes in 1 DMR failed.";
# 	say "\$x = $x";
 #	die;
 #}
#test long chr names

#execute("./new -l long_chr -i $dir/Merged_9966/debug.long_chr.tab,$dir/Merged_9967/debug.long_chr.tab,$dir/Merged_9968/debug.long_chr.tab $dir/Merged_9969/debug.long_chr.tab,$dir/Merged_9970/debug.long_chr.tab,$dir/Merged_9971/debug.long_chr.tab");

#execute("gcc -g -O4 -o new new.c -Wall -pedantic -std=gnu11 -lm -Wextra -Wconversion");
#execute("valgrind -v --tool=memcheck --log-file=long_chr.log --max-stackframe=9951448 --leak-check=full --show-reachable=yes --track-origins=yes /home/con/DNA_Methylation/debug/new -debug -l long_chr -i $dir/Merged_9966/debug.long_chr.tab,$dir/Merged_9967/debug.long_chr.tab,$dir/Merged_9968/debug.long_chr.tab $dir/Merged_9969/debug.long_chr.tab,$dir/Merged_9970/debug.long_chr.tab,$dir/Merged_9971/debug.long_chr.tab");

execute("/usr/bin/time $new -a rn4.refFlat.txt -l new -L control,FeDef -i ../Merged_9966/bs_seeker-CG.tab,../Merged_9967/bs_seeker-CG.tab,../Merged_9968/bs_seeker-CG.tab ../Merged_9969/bs_seeker-CG.tab,../Merged_9970/bs_seeker-CG.tab,../Merged_9971/bs_seeker-CG.tab");
my $x = diff 'control_vs_FeDef_new_c10_CpN5_d1_p0.05_P10.tsv', 'control_vs_FeDef_defiant_c10_CpN5_d1_p0.05_P10.tsv';

if ($x =~ m/^$/) {
	say "rat data suceeded.\n";
} else {
	say "rat data failed.\n";
	say $x;
	die
}

foreach my $option ("-a $annotation_file", '-heatmap', '-l labeltest1234', '-L labeltestgroup1,labeltestgroup2', '-c 3,5,1','-v', '-b', '-c 4', '-CpN 4', '-d 2', '-E','-fdr bonferroni', '-fdr holm', '-fdr hochberg', '-fdr by','-G 1000', '-p 0.04', '-P 7', '-q 7000', '-R', '-s 1', '-U', '-N', '-v holm', '-v hochberg', '-v by', '-v hommel', '-v bonferroni', '-f', '-x CpG') {
   my $label = $option;
   $label    =~ s/-//g;#remove first '-', 'g' option so I can combine options
   $label    =~ s/\//_/g;#remove / as it may confuse the computer
   $label    =~ s/ /_/g;#remove all spaces
	if (-e "$label.log") {
		unlink "$label.log";#prevent reading a file from a previous run
	}
	execute("valgrind -v --tool=memcheck --log-file=$label.log --max-stackframe=9951448 --leak-check=full --show-reachable=yes --track-origins=yes /home/con/DNA_Methylation/debug/new $option -debug -i $dir/Merged_9966/debug.tab,$dir/Merged_9967/debug.tab,$dir/Merged_9968/debug.tab $dir/Merged_9969/debug.tab,$dir/Merged_9970/debug.tab,$dir/Merged_9971/debug.tab");
	my $check = 0;
	open my $log, '<', "$label.log";
	while (<$log>) {
		if (/==\d+==\s+ERROR SUMMARY: 0 errors from 0 contexts \(suppressed: 0 from 0\)/) {
			$check++;
		}
	}
	close $log;
	if ($check != 2) {
		say "$option failed.";
		die;
	}
	say "\n$option succeeded.\n";
	$label .= '_no_DMR';
	if (-e "$label.log") {
		unlink "$label.log";#prevent reading a file from a previous run
	}
	execute("valgrind -v --tool=memcheck --log-file=$label.log --max-stackframe=9951448 --leak-check=full --show-reachable=yes --track-origins=yes ./new $option -debug -i $dir/Merged_9966/debug_noDMR.tab,$dir/Merged_9967/debug_noDMR.tab,$dir/Merged_9968/debug_noDMR.tab $dir/Merged_9969/debug_noDMR.tab,$dir/Merged_9970/debug_noDMR.tab,$dir/Merged_9971/debug_noDMR.tab");
	$check = 0;
	open $log, '<', "$label.log";
	while (<$log>) {
		if (/==\d+==\s+ERROR SUMMARY: 0 errors from 0 contexts \(suppressed: 0 from 0\)/) {
			$check++;
		}
	}
	close $log;
	if ($check != 2) {
		say "$option failed.";
		die;
	}
	unlink "$label.log";
	say "\n$option succeeded.\n";
}
