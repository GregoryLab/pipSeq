#!/usr/bin/perl -w
# randomly shuffle read alignments between sample/control (used to generate background for CSAR peaks)
# inputs are BAM files

use strict;

my $usage = "Usage: perl shuffle_reads_BAM.pl sample_fn control_fn shuffled_sample_fn shuffled_control_fn\n";
my $sample_file = shift or die $usage;
my $control_file = shift or die $usage;
my $shuffled_sample_file = shift or die $usage;
my $shuffled_control_file = shift or die $usage;

my $shuffled_sample_temp = $shuffled_sample_file;
my $shuffled_control_temp = $shuffled_control_file;
$shuffled_sample_temp =~ s/\.bam/\.sam/;
$shuffled_control_temp =~ s/\.bam/\.sam/;

my $shuffled_sample_temp2 = $shuffled_sample_file;
my $shuffled_control_temp2 = $shuffled_control_file;
$shuffled_sample_temp2 =~ s/\.bam/\.unsorted.bam/;
$shuffled_control_temp2 =~ s/\.bam/\.unsorted.bam/;

system("samtools view -H $sample_file > $shuffled_sample_temp");
system("samtools view -H $control_file > $shuffled_control_temp");

open(SAMPLE_IN, "samtools view $sample_file |") || die "Unable to open: $!";
open(CONTROL_IN, "samtools view $control_file |") || die "Unable to open: $!";
open(SHUFFLED_SAMPLE_OUT, ">>$shuffled_sample_temp") || die "Unable to write to: $!";
open(SHUFFLED_CONTROL_OUT, ">>$shuffled_control_temp") || die "Unable to write to: $!";

# randomly mix alignments from sample and control
while (my $line = <SAMPLE_IN>) {
	chomp $line;
	if (int(rand(2)) == 1) {
		print SHUFFLED_SAMPLE_OUT "$line\n";
	}
	else {
		print SHUFFLED_CONTROL_OUT "$line\n";
	}
}

while (my $line = <CONTROL_IN>) {
	chomp $line;
	if (int(rand(2)) == 1) {
		print SHUFFLED_SAMPLE_OUT "$line\n";
	}
	else {
		print SHUFFLED_CONTROL_OUT "$line\n";
	}
}

close(SAMPLE_IN);
close(CONTROL_IN);
close(SHUFFLED_SAMPLE_OUT);
close(SHUFFLED_CONTROL_OUT);


system("samtools view -bS $shuffled_sample_temp > $shuffled_sample_temp2");
system("samtools view -bS $shuffled_control_temp > $shuffled_control_temp2");
system("rm $shuffled_sample_temp");
system("rm $shuffled_control_temp");

$shuffled_sample_file =~ s/\.bam//;
$shuffled_control_file =~ s/\.bam//;
system("samtools sort $shuffled_sample_temp2 $shuffled_sample_file");
system("samtools sort $shuffled_control_temp2 $shuffled_control_file");
system("rm $shuffled_sample_temp2");
system("rm $shuffled_control_temp2");




