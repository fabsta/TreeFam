#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;

#Name of fasta file
my $name = "seq";

#Number of sequences to split fasta by
# we have ~ 2,000,000 sequence
# so we split it into 10 blocks
my $number = 633;

#Set up counters
my $i=0;
my $j=1;

#Name of first file handle
my $first = $name . $j;

#Split fasta file
print STDERR "Splitting $name into files of $number sequences\n";
open(FH, ">$first") or die "Couldn't open fh to $first, $!";
open(PFAMSEQ, $name) or die "Couldn't open $name, $!";  

while(<PFAMSEQ>) {
  if(/^>/){
    $i++;
    if($i > $number) {
      close FH; #Close file handle
      $j++; 
      my $p = $name . $j; #Name of next file handle
      open(FH, ">$p") or die "couldn't open fh to $p, $!";
      $i=0;
    }
  }
  print FH $_;
}
close PFAMSEQ;
close FH;

#Run pfam_scan.pl on farm using job arrays, submit to long queue, select 2Gb memory
print STDERR "Submitting $j pfam_scan.pl jobs to farm\n";
my $num_files = "1-$j";
my $fh = IO::File->new();

$fh->open( "| bsub -q long -o $name.err -Jpfamseq\"[$num_files]\" -R \"select[type==X86_64 && mem>2000] rusage[mem=2000]\" -M 2000000") or die "Couldn't open file handle\n";
$fh->print( "pfam_scan.pl -dir /lustre/scratch101/blastdb/Supported -fasta $name\$\{LSB_JOBINDEX} > $name\$\{LSB_JOBINDEX}.out\n"); 
$fh->close;
