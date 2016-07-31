#!/usr/bin/perl
use strict;
use warnings;

my $s_filter = shift;
if(!defined $s_filter)
{
  print STDERR "Desc : Get supporting backspliced junction reads number for circRNA_finder's output circRNA candidates\n\n";
  print STDERR "Usage: perl $0 sampleId_s_filteredJunctions.bed\n";
  print STDERR "E.g. : perl $0 positive_s_filteredJunctions_circRNA-finder.bed 1>1_circ.num 2>2_circ.num\n"; 
  exit 1;
}

my %circId_count;
open FILTER, "<", $s_filter or die $!;
while(<FILTER>)
{
  chomp;
  my($chr, $start, $end, $count) = (split(/\t/))[0, 1, 2, 4];
  next if($chr eq "chrM");
  $start = $start + 1; #bed format, adjust to 1-based
  my $circId = "$chr:$start|$end";
  if(!exists($circId_count{$circId}))
  {
    $circId_count{$circId} = $count;
  }
  else
  {
  	$circId_count{$circId} += $count;
  }
}
close FILTER;

foreach my $circId(sort keys%circId_count)
{
  my $count = $circId_count{$circId};
  print "$circId\t$count\n";
  print STDERR "$circId\t$count\n" if($count > 1);
}
