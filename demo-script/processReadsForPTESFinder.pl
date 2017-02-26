#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;
if(scalar(@ARGV) != 3)
{
  print STDERR "Desc : Merge fastq_1 & fastq_2 into one file as the input of PTESFinder.\n";
  print STDERR "       If read identifiers are identical for the mates of a read pair(check it yourself before you run this program), add suffix to distinguish them\n";
  print STDERR "       If read identifiers are already distinguishable for the mates, it's still OK to run this program to merge them, but you can `zcat or cat` them directly, instead.\n";
  print STDERR "Usage: perl $0 fastq_1.gz fastq_2.gz merged.fastq\n";
  exit 1;
}
my $mergedFQ = pop(@ARGV);
open MERGE, ">", $mergedFQ or die $!;
for(my $i = 0; $i < scalar(@ARGV); $i++)
{
  my $fastq = $ARGV[$i];
  open FQ, "<:gzip", $fastq or die $!;
  my $suffix = $i + 1;
  while(<FQ>)
 {
    my $readId = (split(/\s+/))[0];
    $readId = "$readId/$suffix";
    my $seq  = <FQ>;
    my $plus = <FQ>;
    my $qual = <FQ>;
    print MERGE "$readId\n";
    print MERGE $seq;
    print MERGE "+\n";
    print MERGE $qual;

  }
  close FQ;
}
close MERGE;
