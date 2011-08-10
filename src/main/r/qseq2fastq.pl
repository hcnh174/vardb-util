#!/usr/bin/perl

use strict;
use warnings;
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_v);

# Usage
my $usage = "
qseq2fastq.pl - converts qseq format files to fastq.
		      by
		Brian J. Knaus
		  May 2010

Usage: perl qseq2fastq.pl options
 required:
  -a	qseq file.
 optional:
  -v	verbose mode [optional T/F, default is F].

";

# command line processing.
getopts('a:v:');
die $usage unless ($opt_a);

my ($inf, $verb);

$inf	= $opt_a if $opt_a;
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Subroutines.

sub outname {
  my $inf = shift;
  my @temp = split("/", $inf);
  $inf = $temp[$#temp];
  @temp = split(/\./, $inf);
  my $outf = $temp[0];
  $outf = join("", $outf, ".fq");

  return($outf);
}

##### ##### ##### ##### #####
# Globals.

my ($in, $name, $outf, $out, $seq, $qual);
my @temp;

##### ##### ##### ##### #####
# Main.

# Manage outfile name.
$outf = outname($inf);

open($in, "<", $inf) or die "Can't open $inf: $!";
open($out, ">", $outf) or die "Can't open $outf: $!";

while(<$in>){
  chomp;
  @temp = split /\t/, $_;
#  print join(":", @temp)."\n";
  $name = $temp[0].":".$temp[2].":".$temp[3].":".$temp[4].":".$temp[5]."\#".$temp[6]."/".$temp[7];
  $seq = $temp[8];
  $qual = $temp[9];
  #
  print $out "\@$name\n";
  print $out "$seq\n";
  print $out "+$name\n";
  print $out "$qual\n";
}

close $in or die "$in: $!";
close $out or die "$out: $!";

##### ##### ##### ##### #####
# EOF.
