#!/usr/bin/perl

use strict;
use warnings;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_i $opt_o);

# Usage
#perl %VARDB_RUTIL_HOME%\fastq2table.pl -i fastq/PXB0219-0011.wk09.fastq -o 219-11.wk9.txt
my $usage = "

fastq2table.pl - converts fastq files to tab-delimited table format
based on fastq2fasta.pl by Brian J. Knaus

Usage: perl fastq2table.pl options
 required:
  -i	fastq input file.
  -o	txt output file.
";

# command line processing.
getopts('i:o:');
die $usage unless ($opt_i);

my ($infile, $outfile);
$infile	= $opt_i if $opt_i;
$outfile	= $opt_o if $opt_o;

print "$infile\n";
print "$outfile\n";

##### ##### ##### ##### #####
# Globals.

my ($temp, $in, $out);
my @temp;

##### ##### ##### ##### #####
# Main.

# Open input and output files.
open( $in, "<",  $infile)  or die "Can't open $infile: $!";
open( $out, ">",  $outfile)  or die "Can't open $outfile $!";

print $out "sequence\n";
while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  chomp($temp[2] = <$in>);	# Third line is an id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

	print $out "$temp[1]\n";
	#print $out "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# EOF.
