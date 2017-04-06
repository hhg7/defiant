#!/usr/bin/env perl

use strict; use warnings; use Cwd;
#---------------------------
#Start Subroutines
#---------------------------
my $TOP_DIRECTORY = getcwd();
local $SIG{__WARN__} = sub {#kill the program if there are any warnings
	my $message = shift;
	my $fail_filename = "$TOP_DIRECTORY/$0.fail";
	open my $fh, '>', $fail_filename or die "Can't write $fail_filename: $!";
	printf $fh ("$message @ %s\n", getcwd());
	close $fh;
	die "$message\n";
};#http://perlmaven.com/how-to-capture-and-save-warnings-in-perl

sub execute {
	my $command = shift;
	print "Executing Command: $command\n";
	if (system($command) != 0) {
		my $fail_filename = "$TOP_DIRECTORY/$0.fail";
		open my $fh, '>', $fail_filename or die "Can't write $fail_filename: $!";
		print $fh "$command failed.\n";
		close $fh;
		print "$command failed.\n";
		die;
	}
}

sub make_2d_gnuplot {
	my $filename = shift;#should be something like 'cpn'
	my $title = shift;
	my $x_label = shift;
	my $png_file_list = shift;
	my $log_y_axis = shift;
	my $filename_base;
	my $Title = ucfirst $title;
	if ($filename =~ m/(\S+)\.tsv$/) {
		$filename_base = $1;#the filename must end in '.tsv'
	} else {
		print "$filename didn't end in '.tsv'\n";
		die;
	}
	$filename_base =~ s/\./_/g;
	$filename_base =~ s/_tsv$/.tsv/;
	open my $gnuplot, '>',  "$filename_base.gnuplot" or die "Can't write $filename_base.gnuplot: $!";
	print $gnuplot "#This file written by $TOP_DIRECTORY/$0\n";
	if ($log_y_axis eq 'yes') {
		print $gnuplot "set logscale y\n";
	}
	print $gnuplot "set xlabel '$x_label' font 'Arial, 16'
set ylabel 'Number of DMRs' font 'Arial, 16'
set title '$title' font 'Arial, 20'
set terminal eps
set xtics font 'Arial, 12'
set key off
set ytics font 'Arial, 12'
set output '$filename_base.eps'\n";

	print $gnuplot "plot '$filename' with lines lw 2 linecolor 'red'\n";
	close $gnuplot;
	execute("gnuplot $filename_base.gnuplot");
	push @{ $png_file_list }, "$filename_base.eps";
}

sub pm3d_gnuplot {
	my $filename = shift;
	my $x_label = shift;
	my $y_label = shift;
	my $image_file_list = shift;
	my $log_cb_axis = shift;
	my $filename_base;
	if ($filename =~ m/(\S+)\.tsv$/) {
		$filename_base = $1;
	} else {
		print "$filename didn't end in '.tsv'\n";
		die;
	}
	$filename_base =~ s/\./_/g;
	$filename_base =~ s/_tsv$/.tsv/;
	if (-z $filename) {#if the file is 0 size, don't plot
		print "$filename has 0 size, so I won't plot.\n";
		unlink $filename;
		return;
	}
	my $isolines = 0;
	open my $tsv_fh, '<', $filename or die "Can't read $filename: $!";
	while (<$tsv_fh>) {
		if (/^$/) {
			$isolines++;
		}
	}
	close $tsv_fh;
	if ($isolines < 2) {
		return;
	}
	my $Title = ucfirst $filename_base;
	$Title =~ s/$TOP_DIRECTORY//;
	$Title =~ s/\/Gnuplot\///;
	$Title =~ s/cpn/CpN/i;
	$Title =~ s/_/ /g;
	open my $gnuplot_fh, '>', "$filename_base.gnuplot" or die "Can't write $filename_base.gnuplot: $!";
	print $gnuplot_fh "#this file written by $TOP_DIRECTORY/$0\n";
	print $gnuplot_fh "set pm3d map
set pm3d corners2color c1\n";
	
	if ($log_cb_axis eq 'yes') {
		print $gnuplot_fh "set logscale zcb\n";
	}
	print $gnuplot_fh "set xlabel '$x_label'
set ylabel '$y_label'
set cblabel 'Number of DMRs'
set title '$Title'
set title font 'Arial, 20'
set xtics font 'Arial, 15'
set ytics font 'Arial, 15'
set key off
set lmargin at screen 0.1
set rmargin at screen 0.7
set xlabel font 'Arial, 18'
set ylabel font 'Arial, 14'
set terminal eps
set output '$filename_base.eps'
splot '$filename'\n";
	close $gnuplot_fh;
	execute("gnuplot $filename_base.gnuplot");
	push @{ $image_file_list }, "$filename_base.eps";
}
use List::Util qw(min max);
sub mean_stdev_array {
	my $array = shift;
	if (scalar @$array == 1) {
		return (0,0,0);
	}
	my ($mean,$standard_deviation) = (0,0);
	foreach my $element (@$array) {
		$mean += $element;
	}
	foreach my $e (@$array) {
		$mean += $e;
	}
	$mean /= scalar @$array;
	foreach my $element (@$array) {
		$standard_deviation += ($element-$mean)*($element-$mean);
	}
	$standard_deviation = sqrt($standard_deviation/(scalar @$array-1));
	my $range = max(@$array) - min(@$array);
	return($mean,$standard_deviation, $range);
}
#---------------------------
##End Subroutines
#---------------------------
#get options

my %default = ('coverage' => 10, 'CpN' => 5, 'd' => 1, 'p' => 0.05, 'percent' => '10.000000', 'skip' => 0);
my $directory = '.';
my $cpn = 'CpN';
my $twoD_graphs_only = 0;
use Getopt::Long;
GetOptions 
				('c:i' => \$default{coverage},#option '-c' is an integer, ':' indicates that it's optional
				 'cpn:i' => \$default{CpN},
				 'd:s' => \$default{d},
				 'directory:s' => \$directory,
				 'p:f' => \$default{p},
				 'P:f' => \$default{percent},
				 'S:i' => \$default{skip},
				 'x:s' => \$cpn,
				 '2d' => \$twoD_graphs_only
				);
#end getting options
my @dmr_count_files;

opendir my $dh, $directory or die "Can't opendir on $directory: $!";
while (my $file = readdir $dh) {
	if ($file =~ m/.*_dmr_count\.tsv$/) {
		print "$file\n";
		push @dmr_count_files, $file;
	}
}
closedir $dh;
mkdir 'Gnuplot';
foreach my $file (@dmr_count_files) {
	my $group;
	if ($file =~ m/(.+)_dmr_count\.tsv$/) {
		$group = $1;
	} else {
		print "Failed to get group name for $file\n";
		die;
	}
	print "$group\n";
	my %number_of_DMRs;
	open my $fh, '<', $file or die "Can't read $file: $!";
	while (<$fh>) {
		if (/^Coverage\s+/) {
			next;#header
		}
		chomp;#remove newline
		my @line = split;#make string into array
		if (scalar @line != 7) {
			print "$_ has an incorrect number of columns.\n";
		}
		$number_of_DMRs{$line[0]}{$line[1]}{$line[2]}{$line[3]}{$line[4]}{$line[5]} = $line[-1];
	}
	close $fh;
	my %filenames;
	my %filehandles;
	my @image_list;
	if ($twoD_graphs_only == 0) {#3D graphs
		my %dmr_counts;
		my %log_cb_axis;
		foreach my $comparison ('coverage_CpN','coverage_d','coverage_p','coverage_percent','coverage_skips',
		'CpN_d','CpN_p','CpN_percent','CpN_skips',
		'd_p','d_percent','d_skip','p_percent','p_skip','percent_skips') {
			$filenames{$comparison} = "$TOP_DIRECTORY/Gnuplot/$group" . "_$comparison.tsv";
			open $filehandles{$comparison}, '>', $filenames{$comparison} or die "Can't write $filenames{$comparison}: $!";
			$log_cb_axis{$comparison} = 'no';
		}

		foreach my $coverage (sort {$a <=> $b} keys %number_of_DMRs) {
#coverage vs. CpN
			foreach my $CpN (sort {$a <=> $b} keys %{$number_of_DMRs{$coverage} }) {
				if (!defined $number_of_DMRs{$coverage}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}) {
					print "\$number_of_DMRs{$coverage}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
					die;
				}
				print {$filehandles{'coverage_CpN'}} "$coverage\t$CpN\t$number_of_DMRs{$coverage}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'coverage_CpN'} }, $number_of_DMRs{$coverage}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}};
			}
			print {$filehandles{'coverage_CpN'}} "\n";
#coverage vs. d
			foreach my $d (sort {$a <=> $b} keys %{ $number_of_DMRs{$coverage}{$default{CpN}} }) {
				if (!defined $number_of_DMRs{$coverage}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}}) {
					print "\$number_of_DMRs{$coverage}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
					die;
				}
				print {$filehandles{'coverage_d'}} "$coverage	$d	$number_of_DMRs{$coverage}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'coverage_d'} }, $number_of_DMRs{$coverage}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}}
			}
			print {$filehandles{'coverage_d'}} "\n";
#coverage vs. p
			foreach my $p (sort {$a <=> $b} keys %{    $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}} }) {
				print {$filehandles{'coverage_p'}} "$coverage	$p	$number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'coverage_p'} }, $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}};
			}
			print {$filehandles{'coverage_p'}} "\n";
#coverage vs. percent
			foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}} }) {
				print {$filehandles{'coverage_percent'}} "$coverage	$percent\t$number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$default{skip}}\n";
				push @{ $dmr_counts{'coverage_percent'} }, $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$default{skip}};
			}
			print {$filehandles{'coverage_percent'}} "\n";
#coverage vs. skips
			foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}}) {
				print {$filehandles{'coverage_skips'}} "$coverage	$skip\t$number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$skip}\n";
				push @{ $dmr_counts{'coverage_skips'} }, $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$skip};
			}
			print {$filehandles{'coverage_skips'}} "\n";
		}
		close $filehandles{'coverage_CpN'};
		close $filehandles{'coverage_d'};
		close $filehandles{'coverage_p'};
		close $filehandles{'coverage_percent'};
		close $filehandles{'coverage_skips'};
		foreach my $graph (keys %dmr_counts) {
			my ($mean, $stddev, $range) = mean_stdev_array(\@{ $dmr_counts{$graph} });
			delete $dmr_counts{$graph};
			if ($stddev < $range/2) {#the data will be very spread out
				$log_cb_axis{$graph} = 'yes';
			}
		}
		undef %dmr_counts;
		pm3d_gnuplot($filenames{'coverage_CpN'}, 'Minimum Coverage (Reads)', "Minimum $cpn (Nucleotides)", \@image_list,	$log_cb_axis{'coverage_CpN'});
		pm3d_gnuplot($filenames{'coverage_d'}, 'Minimum Coverage (Reads)', "Minimum Diff. $cpn (Nucleotides)", \@image_list,	$log_cb_axis{'coverage_d'});
		pm3d_gnuplot($filenames{'coverage_p'}, 'Minimum Coverage (Reads)', 'Maximum p-value', \@image_list,	$log_cb_axis{'coverage_p'});
		pm3d_gnuplot($filenames{'coverage_percent'}, 'Minimum Coverage (Reads)', 'Minimum Percent Change', \@image_list,	$log_cb_axis{'coverage_percent'});
		pm3d_gnuplot($filenames{'coverage_skips'}, 'Minimum Coverage (Reads)', 'Maximum Consecutive Skips (Nucleotides)', \@image_list,	$log_cb_axis{'coverage_skips' });
		foreach my $CpN (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}} }) {
			foreach my $d (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$CpN} }) {
				if (!defined $number_of_DMRs{$default{coverage}}{$CpN}{$d}{$default{p}}{$default{percent}}{$default{skip}} ) {
					print "\$number_of_DMRs{$default{coverage}}{$CpN}{$d}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
					die;
				}
				print { $filehandles{'CpN_d'} } "$CpN	$d	$number_of_DMRs{$default{coverage}}{$CpN}{$d}{$default{p}}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'CpN_d'} }, $number_of_DMRs{$default{coverage}}{$CpN}{$d}{$default{p}}{$default{percent}}{$default{skip}};
			}
			print { $filehandles{'CpN_d'} } "\n";
			foreach my $p (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}} } ) {
				print { $filehandles{'CpN_p'} } "$CpN	$p\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'CpN_p'} }, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}};
			}
			print { $filehandles{'CpN_p'} } "\n";
#CpN vs. percent
			foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}} }) {
				print {$filehandles{'CpN_percent'}} "$CpN\t$percent	$number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$percent}{$default{skip}}\n";
				push @{ $dmr_counts{'CpN_percent'}}, $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$percent}{$default{skip}};
			}
			print {$filehandles{'CpN_percent'}} "\n";
#CpN vs. skips
			foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}}) {
				print {$filehandles{'CpN_skips'}} "$CpN	$skip\t$number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$skip}\n";
				push @{ $dmr_counts{'CpN_skips'} }, $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$skip};
			}
			print {$filehandles{'CpN_skips'}} "\n";
		}
		close $filehandles{'CpN_d'};
		close $filehandles{'CpN_p'};
		close $filehandles{'CpN_percent'};
		close $filehandles{'CpN_skips'};
		foreach my $graph (keys %dmr_counts) {
			my ($mean, $stddev, $range) = mean_stdev_array(\@{ $dmr_counts{$graph} });
			delete $dmr_counts{$graph};
			if ($stddev < $range/2) {#the data will be very spread out
				$log_cb_axis{$graph} = 'yes';
			}
		}
		undef %dmr_counts;
		pm3d_gnuplot($filenames{'CpN_d'}, "Minimum $cpn (Nucleotides)", "Minumum diff. $cpn", \@image_list,	$log_cb_axis{'CpN_d'});
		pm3d_gnuplot($filenames{'CpN_p'}, "Minimum $cpn (Nucleotides)", 'Maximum p-value', \@image_list,	$log_cb_axis{'CpN_p'});
		pm3d_gnuplot($filenames{'CpN_percent'}, "Minimum $cpn (Nucleotides)", 'Minimum Percent Change', \@image_list,	$log_cb_axis{'CpN_percent'});
		pm3d_gnuplot($filenames{'CpN_skips'}, "Minimum $cpn (Nucleotides)", 'Maximum Consecutive Skips (Nucleotides)', \@image_list,	$log_cb_axis{'CpN_skips'});
#d
		foreach my $d (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}} }) {
			foreach my $p (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d} }) {
				print {$filehandles{'d_p'} } "$d\t$p	$number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$p}{$default{percent}}{$default{skip}}\n";
				push @{ $dmr_counts{'d_p'} }, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$p}{$default{percent}}{$default{skip}};
			}
			print {$filehandles{'d_p'} } "\n";
#d vs. percent
			foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}} }) {
				print {$filehandles{'d_percent'}} "$d	$percent\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$percent}{$default{skip}}\n";
				push @{ $dmr_counts{'d_percent'} }, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$percent}{$default{skip}};
			}
			print {$filehandles{'d_percent'}} "\n";
#d vs. skip
			foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}} }) {
				print {$filehandles{'d_skip'}} "$d	$skip	$number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$skip}\n";
				push @{ $dmr_counts{'d_skip'}}, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$skip};
			}
			print {$filehandles{'d_skip'}} "\n";
		}
		close $filehandles{'d_p'};
		close $filehandles{'d_percent'};
		close $filehandles{'d_skip'};
		foreach my $graph (keys %dmr_counts) {
			my ($mean, $stddev, $range) = mean_stdev_array(\@{ $dmr_counts{$graph} });
			delete $dmr_counts{$graph};
			if ($stddev < $range/2) {#the data will be very spread out
				$log_cb_axis{$graph} = 'yes';
			}
		}
		undef %dmr_counts;
		pm3d_gnuplot($filenames{'d_p'}, "Minimum $cpn (Nucleotides)", "Maximum p-value", \@image_list,	$log_cb_axis{'d_p'});
		pm3d_gnuplot($filenames{'d_percent'}, "Minimum $cpn (Nucleotides)", "Minimum Percent Change", \@image_list,	$log_cb_axis{'d_percent'});
		pm3d_gnuplot($filenames{'d_skip'}, "Minimum $cpn (Nucleotides)", "Maximum Consecutive Skips (Nucleotides)", \@image_list,	$log_cb_axis{'d_skip'});
#p-value
		foreach my $p (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}} } ) {
			foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p} }) {
				print {$filehandles{'p_percent'}} "$p	$percent\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$percent}{$default{skip}}\n";
				push @{ $dmr_counts{'p_percent'} }, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$percent}{$default{skip}};
			}
			print {$filehandles{'p_percent'}} "\n";
#p vs. skips
			foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}}) {
				print {$filehandles{'p_skip'}} "$p	$skip\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$skip}\n";
				push @{ $dmr_counts{'p_skip'}}, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$skip};
			}
			print {$filehandles{'p_skip'}} "\n";
		}
		foreach my $graph (keys %dmr_counts) {
			my ($mean, $stddev, $range) = mean_stdev_array(\@{ $dmr_counts{$graph} });
			delete $dmr_counts{$graph};
			if ($stddev < $range/2) {#the data will be very spread out
				$log_cb_axis{$graph} = 'yes';
			}
		}
		undef %dmr_counts;
		close $filehandles{'p_percent'};
		close $filehandles{'p_skip'};
		pm3d_gnuplot($filenames{'p_percent'}, 'Maximum p-value', 'Minimum Percent Change', \@image_list,	$log_cb_axis{'p_percent'});
		pm3d_gnuplot($filenames{'p_skip'}, 'Maximum p-value', 'Maximum Consecutive Skips (Nucleotides)', \@image_list,	$log_cb_axis{'p_skip'});
#percent vs. Skips
		foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}} }) {
			foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$percent} }) {
				print {$filehandles{'percent_skips'}} "$percent\t$skip	$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$skip}\n";
				push @{ $dmr_counts{'percent_skips'}}, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$skip};
			}
			print {$filehandles{'percent_skips'}} "\n";
		}
		close $filehandles{'percent_skips'};
		foreach my $graph (keys %dmr_counts) {
			my ($mean, $stddev, $range) = mean_stdev_array(\@{ $dmr_counts{$graph} });
			delete $dmr_counts{$graph};
			if ($stddev < $range/2) {#the data will be very spread out
				$log_cb_axis{$graph} = 'yes';
			}
		}
		undef %dmr_counts;
		pm3d_gnuplot($filenames{'percent_skips'}, 'Minimum Percent Change','Maximum Consecutive Skips (Nucleotides)', \@image_list,	$log_cb_axis{'percent_skips'});
		undef %filenames;
		undef %filehandles;
	}
#----------------------
#START 2D PLOTS
#----------------------
	my %log_y_axis;
	foreach my $comparison ('coverage','CpN','d','p','percent','skip') {#all different than foreach loop above
		$filenames{$comparison} = $comparison;
		$filenames{$comparison} =~ s/\./_/g;#latex can't figure out file endings if there are dots in the name
		$filenames{$comparison} .= "_2d.tsv";
		$filenames{$comparison} =~ s/ /_/g;
		$filenames{$comparison} = "$TOP_DIRECTORY/Gnuplot/$group" . "_$filenames{$comparison}";
		print "$filenames{$comparison}\n";
		open $filehandles{$comparison}, '>', $filenames{$comparison} or die "Can't write $filenames{$comparison}: $!";
		$log_y_axis{$comparison} = 'no';
	}
#coverage
	my @number_of_DMR;
	foreach my $coverage (sort {$a <=> $b} keys %number_of_DMRs) {
		if (!defined $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}) {
			print "\$number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
			die;
		}
		print {$filehandles{'coverage'}} "$coverage	$number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}\n";
		push @number_of_DMR, $number_of_DMRs{$coverage}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}};
	}
	close $filehandles{'coverage'};
	my ($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{coverage} = 'yes';
	}
	undef @number_of_DMR;
#CpN
	foreach my $CpN (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}} }) {
		if (!defined $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}) {
		 	print "\$number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
			die;
		}
		print {$filehandles{'CpN'}} "$CpN	$number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}}\n";
		push @number_of_DMR, $number_of_DMRs{$default{coverage}}{$CpN}{$default{d}}{$default{p}}{$default{percent}}{$default{skip}};
	}
	close $filehandles{'CpN'};
	($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{CpN} = 'yes';
	}
	undef @number_of_DMR;
#d
	foreach my $d (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}} }) {
		if (!defined $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}}) {
			print "\$number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}} isn't defined.\n";
			die;
		}
		print {$filehandles{'d'}} "$d\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}}\n";
		push @number_of_DMR, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$d}{$default{p}}{$default{percent}}{$default{skip}};
	}
	close $filehandles{'d'};
	($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{d} = 'yes';
	}
	undef @number_of_DMR;
#p-value
	foreach my $p (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}} }) {
		print {$filehandles{'p'}} "$p\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}}\n";
		push @number_of_DMR, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$p}{$default{percent}}{$default{skip}};
	}
	close $filehandles{'p'};
	($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{p} = 'yes';
	}
#percent
	undef @number_of_DMR;
	foreach my $percent (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}} }) {
		print {$filehandles{'percent'}} "$percent\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$default{skip}}\n";
		push  @number_of_DMR, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$percent}{$default{skip}};
	}
	close $filehandles{'percent'};
	($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{percent} = 'yes';
	}
	undef @number_of_DMR;
	foreach my $skip (sort {$a <=> $b} keys %{ $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}}) {
		print {$filehandles{'skip'}} "$skip\t$number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$skip}\n";
		push @number_of_DMR, $number_of_DMRs{$default{coverage}}{$default{CpN}}{$default{d}}{$default{p}}{$default{percent}}{$skip};
	}
	close $filehandles{'skip'};
	($mean, $stddev, $range) = mean_stdev_array(\@number_of_DMR);
	if ($stddev < $range/2) {#the data will be very spread out
		$log_y_axis{skip} = 'yes';
	}
	undef @number_of_DMR;
	make_2d_gnuplot($filenames{'coverage'}, 'Coverage', 'Coverage (Reads)', \@image_list, $log_y_axis{coverage});
	make_2d_gnuplot($filenames{'CpN'}, $cpn, "Min. $cpn (Nucleotides)", \@image_list, $log_y_axis{CpN});
	make_2d_gnuplot($filenames{'d'}, "d", "Min. Differential $cpn", \@image_list, $log_y_axis{d});
	make_2d_gnuplot($filenames{'p'}, "Max. p-value", 'p-value', \@image_list, $log_y_axis{p});
	make_2d_gnuplot($filenames{'percent'}, "Minimum Percent Change", 'Minimum Percent Change', \@image_list, $log_y_axis{percent});
	make_2d_gnuplot($filenames{'skip'}, "Maximum Consecutive Skips", 'Maximum Consecutive Skips (Nucleotides)', \@image_list, $log_y_axis{skip});

	my $latex_filename = "$TOP_DIRECTORY/$group.tex";
	my $title = $group;
	$title =~ s/_/ /g;
	open my $tex, '>', $latex_filename or die "Can't write $latex_filename: $!";
	print $tex "\%This file written by $TOP_DIRECTORY/$0\n";#write a comment in the file so I can find the creating perl script later
	print $tex "\\documentclass{article}
\\title{$title Images}
\\usepackage{graphicx,placeins}
\\usepackage{epstopdf}
\\begin{document}
\\maketitle\n";
	print $tex "\\begin{table}[thp]
	\\centering
	\\begin{tabular}{|c|c|}\\hline
		\\textbf{Option} & \\textbf{Setting}\\\\ \\hline\n";
	foreach my $default (sort {$a cmp $b} keys %default) {
		if ($default eq 'CpN') {
			print $tex "\\textbf{$cpn}	&	$default{CpN}\\\\ \n";
		} else {
			print $tex "\\textbf{$default} & $default{$default}\\\\\n";
		}
	}
	print $tex "	\\hline\\end\{tabular}
	
	\\caption{All images show two parameters varied as all others are set to default values, which are shown above.}
	\\label{tab:parameter_cutoffs}
\\end{table}\n";
	foreach my $figure (sort {$a cmp $b} @image_list) {
		print $tex "\\begin{figure}[htp]
	\\centering
	\\includegraphics[width=\\textwidth]{$figure}
	\\caption{}
\\end{figure}
\\FloatBarrier\n";
	}
	print $tex "\\end{document}\n";
	close $tex;
	execute("xelatex $latex_filename");
	if ($latex_filename =~ m/(\S+)\.tex$/) {
		unlink "$1.aux";#remove clutter
	}
}
