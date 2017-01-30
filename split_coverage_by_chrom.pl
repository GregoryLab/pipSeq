#!/usr/bin/perl -w
use strict;

# split a 3-column input BED file of (coverage) values by chromosome

my $usage = "Usage: perl $0 fn out_prefix";
my $infile = shift or die $usage;
my $outprefix= shift or die $usage;

open(IN, "<$infile") || die "Unable to open $infile: $!";
my $current_chrom = "";

while (my $line = <IN>) {
	chomp $line;
	my ($chrom, $pos, $value) = split(/\t/, $line);
	if ($chrom ne $current_chrom) {
		# open a new file
		if ($current_chrom ne "") {
			close(OUT);
		}
		open(OUT, ">$outprefix.$chrom.coverage.txt") || die "Unable to write to: $!";
		$current_chrom = $chrom;
	}
	printf OUT "$value\n";
}
close(OUT);


