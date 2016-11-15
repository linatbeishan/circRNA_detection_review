#!/usr/bin/perl
use strict;
use warnings;

my $input = shift;
my $thresholdRead1 = shift;
my $thresholdRead2 = shift;
if(!defined $input || !defined $thresholdRead1 || !defined $thresholdRead2)
{
  print STDERR "perl $0 KNIFE.reports p_predicted_threshold_for_read1 p_predicted_threshold_for_read2\n";
  print STDERR "threshold used: positive: 0-0, mixed-polyA: 0-0, background-simu(neg): 0.99-0.99, mixed-simu(pos+neg): 0.99-0.99, background-polyA: 0.97-0.97, Hela_plus:0.63-0.73, Hela_minus: 0.78-0.89, Hs68_minus: 0.78-0.75, Hs68_plus: 0.64-0.65, according to the cumalative distribution plot of posterior probability\n";
  exit 1;
}

open INPUT, "<", $input or die $!;
while(<INPUT>)
{
  next unless (/^chr/);
  chomp;
  my ($locus, $p_predictedR1, $p_predictedR2, $totalReads) = (split(/\t/))[0, 2, 4, 5];
  next if($p_predictedR1 < $thresholdRead1 && $p_predictedR2 < $thresholdRead2);
  my @fields = split(/\|/, $locus);
  my $chr = $fields[0];
  my $pos1  = (split(/:/, $fields[1]))[1];
  my $pos2  = (split(/:/, $fields[2]))[1];
  if($pos1 <= $pos2)
  {
    $pos1 += 1;
    print "$chr:$pos1|$pos2\t$totalReads\n";
    print STDERR "$chr:$pos1|$pos2\t$totalReads\n" if($totalReads > 1);
  }
  else
  {
    $pos2 += 1;
    print "$chr:$pos2|$pos1\t$totalReads\n";
    print STDERR "$chr:$pos2|$pos1\t$totalReads\n" if($totalReads > 1);
  }
}
close INPUT;

