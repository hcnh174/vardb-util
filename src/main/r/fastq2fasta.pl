#!/usr/bin/perl

use strict;
use warnings;
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_i $opt_o $opt_c $opt_v);

# Usage
my $usage = "

fastq2fasta.pl - converts fastq files to fasta format.
                  
		      by
		Brian J. Knaus
		 December 2009

Copyright (c) 2009 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl script.pl options
 required:
  -i	fastq input file.
  -o	txt output file.
 optional:
  -c	inset in nucleotides, used to remove barcodes [default = 0].
  -v	verbose mode [T/F, default is F].

";

# command line processing.
getopts('i:o:c:v:');
die $usage unless ($opt_i);

my ($infile, $inset, $verb);
$infile	= $opt_i if $opt_i;
$outfile	= $opt_o if $opt_o;
$inset	= $opt_c ? $opt_c : 0;
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Globals.

my ($temp, $in, $out, $outfile);
my @temp;

##### ##### ##### ##### #####
# Manage outfile name.
# @temp = split(/\//, $infile);
# $temp = $temp[$#temp];
# @temp = split(/\./, $temp);
# $outfile = $temp[0];

##### ##### ##### ##### #####
# Main.

# Open input and output files.
open( $in, "<",  $infile)  or die "Can't open $infile: $!";
open( $out, ">",  $outfile.".fa")  or die "Can't open $outfile.fa: $!";

while (<$in>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <$in>);	# Second line is a sequence.
  chomp($temp[2] = <$in>);	# Third line is an id.
  chomp($temp[3] = <$in>);	# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

  # Substring to inset value.
  $temp[1] = substr($temp[1], $inset);

  # Print to fasta file.
  print $out ">$temp[0]\n";
  print $out "$temp[1]\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# EOF.
