#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use PerlIO::gzip;

my $ofastq = shift;
my $nfastq = shift;
if(!defined $ofastq)
{
	print STDERR "perl $0 fastq\n";
	exit 1;
}

$nfastq = defined($nfastq)? $nfastq : "new.fastq.gz";

if($ofastq =~/\.gz/)
{
	open OFQ, "<:gzip", $ofastq or die $!;
}
else
{
	open OFQ, "<", $ofastq or die $!;
}
open NFQ, ">:gzip", $nfastq or die $!;
my $count = 0;
while(<OFQ>)
{
	chomp;
	my $readType = (split(/\//))[1];
	$readType = (defined $readType)? "/$readType" : "";
	$count++;
	my $readId = "\@background:".$count.$readType;
	my $seq  = <OFQ>;
	my $plus = <OFQ>;
	my $qual = <OFQ>;
	print NFQ "$readId\n";
	print NFQ "$seq";
	print NFQ "$plus";
	print NFQ "$qual";
}
print "Total reads number:\t$count\n";
close OFQ;
close NFQ;
