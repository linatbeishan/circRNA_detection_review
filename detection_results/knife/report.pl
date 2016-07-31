#!/usr/bin/perl
use strict;
use warnings;

my $input = shift;
if(!defined $input)
{
  print STDERR "Desc : Get supporting backspliced junction reads number for KNIFE's output circRNA candidates\n\n";
  print STDERR "Warn : This program is ONLY for samples: positive and mixed, filtering by naive algorithm instead of GLM,\n";
  print STDERR "       since GLM doesn't seem to work well, when we check the cumulative posterior probability plot\n";
  print STDERR "Usage: perl $0 KNIFE.reports\n";
  print STDERR "E.g. : perl $0 positive_knife_report.txt 1>1_circ.num 2>2_circ.num or perl $0 mixed_knife_report.txt 1>1_circ.num 2>2_circ.num\n";
  exit 1;
}

open INPUT, "<", $input or die $!;
while(<INPUT>)
{
  next unless(/^chr/);
  chomp;
  my($locus, $linear, $circ, $p_value) = (split(/\t/))[0, 1, 5, 7];
  next if( $locus =~ /\|reg\|/);
  next if( $p_value eq "-"  || $p_value < 0.9 || $linear / $circ > 0.1);
  my @fields = split(/\|/, $locus);
  my $chr = $fields[0];
  my $start = (split(/:/, $fields[1]))[1] + 1;
  my $end   = (split(/:/, $fields[2]))[1];
  print "$chr:$start|$end\t$circ\n";
  print STDERR "$chr:$start|$end\t$circ\n" if($circ > 1);
}
close INPUT;

