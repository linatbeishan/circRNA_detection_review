#!/usr/bin/perl
use strict;
use warnings;

my $input = shift;
my $threshold = shift;
if(!defined $input || !defined $threshold)
{
  print STDERR "Desc : Get supporting backspliced junction reads number for KNIFE's output circRNA candidates\n";
  print STDERR "       Threshold used: polyA: 0.9, RNaseR_plus:0.62, RNaseR_minus: 0.78, according to the cumalative distribution plot of posterior probability\n\n";
  print STDERR "Warn : This program is ONLY for samples: polyA, RNaseR-plus and RNaseR-minus\n";
  print STDERR "Usage: perl $0 KNIFE.glmreports p_predicted_threshold\n";
  print STDERR "E.g. : perl $0 Hela_RNaseR-minus_knife__circJuncProbs.txt 0.78 1>1_circ.num 2>2_circ.num\n";
  exit 1;
}

open INPUT, "<", $input or die $!;
while(<INPUT>)
{
  next unless (/^chr/);
  chomp;
  my ($locus, $numReads, $p_predicted) = (split(/\t/))[0, 1, 2];
  next if($p_predicted < $threshold);
  my @fields = split(/\|/, $locus);
  my $chr = $fields[0];
  my $start = (split(/:/, $fields[1]))[1] + 1;
  my $end   = (split(/:/, $fields[2]))[1];
  print "$chr:$start|$end\t$numReads\n";
  print STDERR "$chr:$start|$end\t$numReads\n" if($numReads > 1);
}
close INPUT;

