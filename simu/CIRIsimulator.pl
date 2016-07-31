	use 5.012;
	use Getopt::Long;
	use strict;
	use warnings;
	my ($fq1, $fq2, $out, $gtf, $coverage, $coverage2, $rand_mode, $rand_mode2, $read_length, $seq_err, $insert_length, $ref_dir, $help, $if_chr, $minCircRNA, $circBaseGeneSet, $selectedTransNumPerGene);
	Getopt::Long::GetOptions(
		'1=s'   =>  \$fq1,
		'2=s'   =>  \$fq2,
		'O=s'   =>  \$out,
		'G=s'   =>  \$gtf,
		'C=i'   =>  \$coverage,
		'LC=i'  =>  \$coverage2,
		'R=i'   =>  \$rand_mode,
		'LR=i'  =>  \$rand_mode2,
		'L=i'   =>  \$read_length,
		'E=i'   =>  \$seq_err,
		'I=i'   =>  \$insert_length,
		'D=s'   =>  \$ref_dir,
		'M=i'   =>  \$minCircRNA,
		'DB=s'   =>  \$circBaseGeneSet,
		'CHR1=i'    =>  \$if_chr,
		'H!'    =>  \$help
	);

	$circBaseGeneSet = (defined $circBaseGeneSet)? $circBaseGeneSet : "./hsa_hg19_circRNA.txt";
	$minCircRNA = (defined $minCircRNA)? $minCircRNA : 100;
	my $minNumJuncReads = 2;

	srand(20160316); #reproducible

	my $if_die = 0;
	#my ($fq1, $fq2, $gtf) = (">>./simulate_80_10X_80bp_200_1.fq",">>./simulate_80_10X_80bp_200_2.fq",">>./gtf_80_10X_80bp_200.out");
	#my $coverage = 10;
	#my $coverage2 = 100;
	#my $rand_mode = 1;
	#my $rand_mode2 = 2;
	#my $read_length = 80;
	my $cRNA_size = 0;
	#my $if_PE = 2;
	#my $seq_err = 1;
	#my $insert_length = 200;
	#my $ref_dir = "/panfs/home/zhao/gaoyuan/bwaphage/hg19/";

	if(defined($help)){
		print "This is CIRI-simulator, a simulation tool for circRNAs. Welcome!\n\n";
		print "Written by Yuan Gao. Any questions please mail to gaoyuan06\@mails.ucas.ac.cn.\n\n";
		print "Arguments (all required):\n";
		print "\t-1\t\toutput simulated PE reads file 1 name\n";
		print "\t-2\t\toutput simulated PE reads file 2 name\n";
		print "\t-O\t\toutput simulated reads list name\n";
		print "\t-G\t\tinput gtf formatted annotation file name(hg19 refGene.gtf downloaded from ucsc genome browser)\n";
		print "\t-DB\t\tcircRNA annotation file from circBase\n";
		print "\t-C\t\tset coverage or max coverage (when choosing -R 2) for circRNAs\n";
		print "\t-LC\t\tset coverage or max coverage (when choosing -LR 2) for linear transcripts\n";
		print "\t-R\t\tset random mode for circRNAs: 1 for constant coverage; 2 for random coverage\n";
		print "\t-LR\t\tset random mode for linear transcripts: 1 for constant coverage; 2 for random coverage\n";
		print "\t-L\t\tread length of simulated reads (e.g. 100)\n";
		print "\t-E\t\tpercentage of sequencing error (e.g. 2)\n";
		print "\t-I\t\tinsertion length (should be larger than read length) (e.g. 350)\n";
		print "\t-D\t\tdirectory of reference sequence(s) (please make sure all references referred in gtf file are included in the directory)\n";
		print "\t-CHR1\t\tif only choose chr1 to simulate sequencing reads: 1 for yes; 0 for no\n";
		print "\t-M\t\tminimum size(genomic distance between start and end) of circRNAs to be generated\n";
		print "\t-H\t\tshow help information\n";
		$if_die = 1;
	}elsif(!defined($fq1) or !defined($fq2) or !defined($out) or !defined($gtf) or!defined($coverage) or !defined($coverage2) or !defined($rand_mode) or !defined($rand_mode2) or !defined($read_length) or !defined($seq_err) or !defined($insert_length) or !defined($ref_dir)or !defined($if_chr)){
		$if_die = 1;
		print "Please input complete arguments.\n";
	}elsif($insert_length <= $read_length){
		$if_die = 1;
		print "Insertion length should be larger than read length.\n";
	}
	die if $if_die == 1;

	my $qualChar = chr( 33 - 10 * log10($seq_err/100) );

	my %chr_tranId_circId;
	my %circId_spLen;
	open CIRCBASE, "<", $circBaseGeneSet or die "cannot open circBase file: $circBaseGeneSet\n";
	<CIRCBASE>;
	while(<CIRCBASE>)
	{
		chomp;
		my ($chr, $start, $end, $tranId, $spLen) = (split(/\t/))[0, 1, 2, 10, 6];
		$start = $start + 1; #due to bed format
		my $circId = "$chr\t$start\t$end";
		$circId_spLen{$circId} = $spLen;
		push(@{$chr_tranId_circId{$chr}{$tranId}}, $circId);
	}

	close CIRCBASE;


	$fq1 = ">>".$fq1;
	$fq2 = ">>".$fq2;
	#$out = ">>./".$out;
	#my @allProteinCodingGenes;
	my %chr_gene_trsc_exon;
	my %dupTranId_yes;
	my $pre_gene = '';
	my @gene_anno;
	my @chr;
	my $seqID;
	my $sim_total = 0;
	$ref_dir = $ref_dir."/" unless rindex($ref_dir, "/") == length($ref_dir) - 1;
	open GTF, "<", $gtf or die "cannot open gtf file: $!";

	system("rm $out") if(-e $out);
	open OUT, ">>", $out or die;

	while(<GTF>){
		chomp;
		next if /^#/;
		my @line = split /\t/;
		last if ($if_chr == 1 and $line[0] ne "chr1");
		my @atr = split(/;/, $line[8]);
		my $geneId = (split(/\"/, $atr[0]))[1];
		my $tranId = (split(/\"/, $atr[1]))[1];
		next if(!exists($chr_tranId_circId{$line[0]}{$geneId}));
		if($tranId =~ /(\S+)_dup(\d+)$/) #don't generate from these transcripts
		{
			$dupTranId_yes{$1} = 1;
		}

		if($pre_gene ne $geneId and $pre_gene ne ''){
			&split_transcript(@gene_anno);
			@gene_anno = ();
		}
		push @gene_anno, $_;
		$pre_gene = $geneId;
	}
	close GTF;
	&split_transcript(@gene_anno); #the last flush
	#for my $chr(@chr){
	#   for my $gene(keys %{$chr_$gene_trsc_exon{$chr}}){
	#       for my $trsc(keys %{$chr_$gene_trsc_exon{$chr}{$gene}}){
	#           for $exon(@{$chr_$gene_trsc_exon{$chr}{$gene}{$trsc}}){
	#               print OUT "$chr\t$gene\t$trsc\t@$exon\n";
	#           }
	#       }
	#   }
	#}
	my $ignored = 0;
	for my $chromo(@chr)
	{
		open CHR, "<", $ref_dir."$chromo.fa" or open CHR, "<", $ref_dir."$chromo.fasta" or die "cannot open the chr fasta file $chromo: $!";
		my $uni_seq = 0;
		my $chr_seq;
		while(<CHR>)
		{
			chomp;
			if(/^>/ and $uni_seq == 0)
			{
				$uni_seq = 1;
			}
			elsif(/^>/)
			{
				die "There are more than one sequence in $chromo file. Please check!";
			}
			else
			{
				$chr_seq .= $_;
			}
		}
		close CHR;
		foreach my $trsc(sort keys %{$chr_tranId_circId{$chromo}})
		{
			my $trsc_seq;
			my $gene = $trsc;
			if(!exists($chr_gene_trsc_exon{$chromo}{$gene}{$trsc})) #some are not in annotation files
			{
				next;
			}
			for my $i(0 .. $#{$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}})
			{
			 	my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$i];
				push @$exon, substr ($chr_seq, $$exon[0]-1, $$exon[1]-$$exon[0]+1);
				push @$exon, substr ($chr_seq, $$exon[0]-3, 2);
				push @$exon, substr ($chr_seq, $$exon[1], 2);
				$trsc_seq .= $$exon[-3];
			}
			&simulate_reads2( $rand_mode2, $trsc_seq, $coverage2 ) if (length($trsc_seq) > $insert_length && $coverage2 > 0);
			foreach my $circId(@{$chr_tranId_circId{$chromo}{$trsc}})
			{
				my $cRNA_seq;
				my ($circStart, $circEnd) = (split(/\t/, $circId))[1,2];

				if(exists($dupTranId_yes{$trsc}))
				{
					$ignored++;
					print STDERR "circRNA:\t$circId ignored due to duplicate mRNA id:$trsc\n";
					next;
				}
				
				if($circEnd - $circStart + 1 < $minCircRNA)
				{
					my $length = $circEnd - $circStart + 1;
					$ignored++;
					print STDERR "circRNA:\t$circId ignored due to sequence too short, $length < $minCircRNA bp\n";
					next;
				}
				my ($startInd, $endInd, $startType, $endType);
####################################################################################################################################################				

				my $numExons = scalar(@{$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}});
				for($startInd = 0; $startInd < $numExons; $startInd++)
				{

					if($circStart <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0])
					{
						last;
					}
				}
				for($endInd = $numExons-1; $endInd >= 0; $endInd--)
				{
					if($circEnd >= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1])
					{
						last;
					}
				}

				if($startInd == 0 || ($circStart > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$numExons-1][1])) #intergenic/first exon begin ->
				{
					$startType = "intergenic or first exon begin ";
					if($endInd == -1 || ($circStart > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$numExons-1][1]) ) #intergenic->intergenic/first exonic
					{
						$cRNA_seq = substr($chr_seq, $circStart - 1, $circEnd - $circStart + 1);
						$endType = "intergenic/first exonic";
					}
					elsif($endInd == $numExons - 1) #intergenic or last exon boundary
					{
						$endType = "intergenic or last exon boundary";
						my $length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart;
						$cRNA_seq = substr($chr_seq, $circStart - 1, $length);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					elsif($circEnd > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]) #intergenic -> exonic
					{
						$endType = "exonic";
						$cRNA_seq = substr($chr_seq, $circStart - 1, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						my $length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]+1;
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0], $length);
					}
					elsif($circEnd <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]) #intergenic -> intron or exon boundary
					{
						$endType = "intronic or exon boundary";
						$cRNA_seq = substr($chr_seq, $circStart - 1, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						my $length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					else
					{
						print STDERR "Start with intergenic-Error: circRNA: $circId, $trsc undefined relation\n";
						exit 1;
					}
				}
				elsif( $circStart <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] || $circStart ==  ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0]) #exonic ->
				{
					$startType = "exonic";
					if($circEnd <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1])
					{
						$endType = "same exon";
						$cRNA_seq = substr($chr_seq, $circStart - 1, $circEnd - $circStart + 1);
					}

					elsif($endInd == $numExons - 1) #exonic -> intergenic
					{
						$endType = "intergenic or last exon boundary";
						my $length;
						if(defined ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] && $circStart <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1])
						{
							$length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] - $circStart + 1;
							$cRNA_seq = substr($chr_seq, $circStart - 1, $length);
						}
						foreach my $j($startInd .. $endInd)
						{
							my $exon = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					elsif($circEnd > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0])#exonic->exonic
					{
						$endType = "exonic";
						my $length;
						if(defined ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] && $circStart <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1])
						{
							$length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] - $circStart + 1;
							$cRNA_seq = substr($chr_seq, $circStart - 1, $length);
						}
						foreach my $j($startInd .. $endInd)
						{
							my $exon = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0] + 1;
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0] - 1, $length);
					}
					elsif($circEnd <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]) #exonic->intronic or exon boundary
					{
						$endType = "intronic or exon boundary";
						my $length;
						if(defined ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] && $circStart <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1])
						{
							$length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] - $circStart + 1;
							$cRNA_seq = substr($chr_seq, $circStart - 1, $length);
						}
						foreach my $j($startInd .. $endInd)
						{
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					else
					{
						print STDERR "Start with EXON-Error: circRNA: $circId, $trsc undefined relation\n";
						exit 1;
					}
				}
				##here
				elsif($circStart > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd-1][1] ) #intronic ->
				{
					$startType = "intronic";
					if($endInd == $startInd || $endInd == $startInd - 1) #within same intron/adajacent exon
					{
						$endType = "same intron or adajacent exon to intron";
						$cRNA_seq = substr($chr_seq, $circStart - 1, $circEnd - $circStart + 1);
					}
					elsif($endInd == $numExons - 1 ) #intron->intergenic/last exon end
					{
						$endType = "intergenic or last exon boundary";
						my $length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart;
						$cRNA_seq = substr($chr_seq, $circStart, $length);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					#NM_015078       circRNA:183005268-183015456
					elsif($circEnd > ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]) #intron -> exon
					{
						$endType = "exonic";
						my $length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart;
						$cRNA_seq = substr($chr_seq, $circStart, $length);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0] + 1;
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0] - 1, $length);
					}
					elsif($circEnd <= ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd+1][0]) #intron -> intron or exon boundary
					{
						$endType = "intronic or exon boundary";
						my $length = ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$startInd][0] - $circStart;
						$cRNA_seq = substr($chr_seq, $circStart, $length);
						foreach my $j($startInd .. $endInd)
						{
							my $exon = $chr_gene_trsc_exon{$chromo}{$gene}{$trsc}[$j];
							$cRNA_seq .= $$exon[-3];
						}
						$length = $circEnd - ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1];
						$cRNA_seq .= substr($chr_seq, ${$chr_gene_trsc_exon{$chromo}{$gene}{$trsc}}[$endInd][1], $length);
					}
					else
					{
						print STDERR "Start with Intron-Error: circRNA: $circId, $trsc undefined relation\n";
						exit 1;
					}
				}
				else
				{
					print STDERR "OOOPS, Error, circRNA: $circId, $trsc startType undefined relation\n";
					exit 1;
				}
				my $cLength = length($cRNA_seq);
				if(abs($cLength -$circId_spLen{$circId}) > 10)
				{
					$ignored++;
					print STDERR "circRNA:\t$circId ignored due to sequence inconsistence: dbLength:$circId_spLen{$circId}, currentLength:$cLength, $chromo\t$gene\t$trsc\tcircRNA:$circStart-$circEnd\tExonIndex:$startInd-$endInd\tjunction Type:$startType-$endType\n";
					next;
				}
#				my $inferedLength = $circEnd - $circStart + 1;
#				my $odd = ($cLength > $inferedLength)? 1 : 0;
#				print STDERR "$odd\t$chromo\t$gene\t$trsc\tcircRNA:$circStart-$circEnd\tbingo:$startInd-$endInd\tjunction Type:$startType-$endType\tlength:$cLength-$inferedLength\n";
				#print ">$chromo:$circStart-$circEnd\n";
				#print "$cRNA_seq\n";

				my $numJuncs = &simulate_reads( $rand_mode, $cRNA_seq, $coverage );
				print OUT "summary:\t$gene\t$trsc\tcircRNA:$chromo:$circStart-$circEnd\tsplicedLength:$cLength\tnumofJunctionReadPairs:$numJuncs\tExonIndex:$startInd-$endInd\tjunction Type:$startType-$endType\n";
				$sim_total++;

=head
				#my $outStat = "$chromo\t$gene\t$trsc\tcircRNA:$circStart-$circEnd\tbingo:$startInd-$endInd\tjunction Type:$startType-$endType\n";

				my ($numJuncs, $fq1_reads, $fq2_reads, $reads_stat) = &simulate_reads( $rand_mode, $cRNA_seq, $coverage );
				if($numJuncs > 1)
				{
					$sim_total++;
					$outStat .= $reads_stat;
					print OUT $outStat;
					print FQ1 $fq1_reads;
					print FQ2 $fq2_reads;

				}
				$outStat = "";
=cut
####################################################################################################################################################
			}
		}
	}
	print STDERR "!!total number of circRNA generated from DB: $sim_total\n";
	print STDERR "!!ignored number of circRNA from DB: $ignored\n";
	close FQ1;
	close FQ2;
	close OUT;

	sub simulate_reads2{
			open FQ1, $fq1 or die;
			open FQ2, $fq2 or die;
			my $mode = shift @_;
			my $trsc_coverage;
			my $seq_length = length($_[0]);
			if ($mode == 1){
					$trsc_coverage = $_[1];
			}else{
					$trsc_coverage = rand($_[1]+1);
			}

			my ($read_num, undef) = sort{$b <=> $a}(int( $seq_length * $trsc_coverage / $read_length / 2 ),1);
			my $err_num = int( $seq_length * $trsc_coverage * $seq_err / 100 );
			my %err_read;
			for (1 .. $err_num){
					my $err_loci = int( rand( ($read_num)*2+1 ) );
					$err_read{$err_loci} ++;
			}
			for my $x( 1 .. $read_num ){
					my $strand = int( rand(2) );
					my ($seq1, $seq2);
					my $ins_len_rand = $insert_length;      #insert length can be simulated later
							my $start_loci = int( rand($seq_length - $ins_len_rand) );
					my $start_loci2 = $start_loci + $ins_len_rand - $read_length;
					$seqID ++;
					if ($strand == 0){
							$seq1 = substr( $_[0], $start_loci, $read_length );
							$seq2 = &comp_rev( substr( $_[0], $start_loci2, $read_length ) );
					}else{
							$seq1 = &comp_rev( substr( $_[0], $start_loci, $read_length ) );
							$seq2 = substr( $_[0], $start_loci2, $read_length );
					}
					my @errs1;
					for (1 .. $err_read{$x * 2 - 1}){
							my $err_loci = int( rand($read_length) );
							redo if $err_loci ~~ @errs1;
							push @errs1, $err_loci;
							my $ori_base = substr( $seq1, $err_loci, 1 );
							substr( $seq1, $err_loci, 1 ) = &simulate_seq_error($ori_base);
					}
					my @errs2;
					for (1 .. $err_read{$x * 2}){
							my $err_loci = int( rand($read_length) );
							redo if $err_loci ~~ @errs2;
							push @errs2, $err_loci;
							my $ori_base = substr( $seq2, $err_loci, 1 );
							substr( $seq2, $err_loci, 1 ) = &simulate_seq_error($ori_base);
					}
					print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
					print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
					print FQ1 "$seq1\n";
					print FQ2 "$seq2\n";
					print FQ1 "+\n";
					print FQ2 "+\n";
					print FQ1 ("$qualChar" x $read_length) . "\n";
					print FQ2 ("$qualChar" x $read_length) . "\n";
			}
	}
	sub simulate_reads
	{
		open FQ1, $fq1 or die;
		open FQ2, $fq2 or die;
		my $mode = shift @_;
		my $seq_length = length($_[0]);
		my $cRNA_coverage;
		my $seq4substr;
		if ($mode == 1)
		{
			$cRNA_coverage = $_[1];
		}
		else
		{
			$cRNA_coverage = rand($_[1]+1);
		}
		if($cRNA_size == 1 or $cRNA_size == 0)
		{
			$seq4substr = $_[0] x 10;   #different according to $seq_length my $seq4substr = $_[0] . substr( $_[0], 0, $read_length-1 );
		}
		else
		{
			$seq4substr = $_[0] x 10;
		}
		my $read_num = int( $seq_length * $cRNA_coverage / $read_length / 2 );
		my $err_num = int( $seq_length * $cRNA_coverage * $seq_err / 100 );
		my %err_read;
		for (1 .. $err_num)
		{
			#modified by linwei@20160316: my $err_loci = int( rand( ($read_num)*2+1 ) ); -> my $err_loci = int( rand($read_num*2) +1 );
			my $err_loci = int( rand($read_num*2) +1 );
			$err_read{$err_loci} ++;
		}
		my $fq1_reads;
		my $fq2_reads;
		my $reads_stat;
		my $numJuncs = 0;

		foreach my $n(1..$minNumJuncReads)
		{
			my ($which_to_cross, $start_loci, $start_loci2, $seq1, $seq2);
			$which_to_cross = int(rand(2));
			if($which_to_cross == 0)
			{
				if($seq_length <= $read_length)
				{
					$start_loci = int(rand($seq_length));
				}
				else
				{
					$start_loci = ($seq_length - $read_length) + int(rand($read_length));
				}
				$start_loci2 = ( $start_loci + $insert_length - $read_length ) % $seq_length;
			}
			else
			{
				if($seq_length <= $read_length)
				{
					$start_loci2 = int(rand($seq_length));
				}
				else
				{
					$start_loci2 = ($seq_length - $read_length) + int(rand($read_length));
				}
				$start_loci = ( $start_loci2 + $read_length - $insert_length) % $seq_length;
			}

			my $strand = int( rand(2) );
			$seqID ++;
			if ($strand == 0)
			{
				$seq1 = substr( $seq4substr, $start_loci, $read_length );
				$seq2 = &comp_rev( substr( $seq4substr, $start_loci2, $read_length ) );	
			}else
			{
				$seq1 = &comp_rev( substr( $seq4substr, $start_loci, $read_length ) );
				$seq2 = substr( $seq4substr, $start_loci2, $read_length );
			}
			my $if_error = int(rand(2));
			if($if_error)
			{
				my $err_num = int(($read_length * 2 * $seq_err) / 100);
				my @errs;
				foreach my $n (1..$err_num)
				{
					my $err_loci = int(rand($read_length*2));
					redo if $err_loci ~~ @errs;
					push(@errs, $err_loci);
					if($err_loci >= $read_length)
					{
						$err_loci = $err_loci % $read_length;
						my $ori_base = substr( $seq2, $err_loci, 1 );
						substr( $seq2, $err_loci, 1) = &simulate_seq_error($ori_base);
					}
					else
					{
						my $ori_base = substr( $seq1, $err_loci, 1 );
						substr( $seq1, $err_loci, 1) = &simulate_seq_error($ori_base);
					}
				}
			}
			print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
			print FQ1 "$seq1\n";
			print FQ1 "+\n";
			print FQ1 ("$qualChar" x $read_length) . "\n";

			print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
			print FQ2 "$seq2\n";
			print FQ2 "+\n";
			print FQ2 ("$qualChar" x $read_length) . "\n";

			print OUT ">\t$n\t$seqID\n";
			if($seq_length - $start_loci < $read_length)
#			if($which_to_cross == 0)
			{
				print OUT "**\t1\t$start_loci\t$start_loci2\n";
			}
#			else
			if($seq_length - $start_loci2 < $read_length)
			{
				print OUT "**\t2\t$start_loci\t$start_loci2\n";
			}

			$numJuncs++;
		}
		for my $x( 1 .. $read_num )
		{
			my $isJunc = 0;
			my $start_loci = int( rand($seq_length) );
			my $strand = int( rand(2) );
			my ($seq1, $seq2);
			my $ins_len_rand = $insert_length;      #insert length can be simulated later
			my $start_loci2 = ( $start_loci + $ins_len_rand - $read_length ) % $seq_length;
			$seqID ++;
			if ($strand == 0)
			{
				$seq1 = substr( $seq4substr, $start_loci, $read_length );
				$seq2 = &comp_rev( substr( $seq4substr, $start_loci2, $read_length ) );	
			}else
			{
				$seq1 = &comp_rev( substr( $seq4substr, $start_loci, $read_length ) );
				$seq2 = substr( $seq4substr, $start_loci2, $read_length );
			}
			my @errs1;
			if(exists($err_read{$x * 2 - 1}))
			{
				for (1 .. $err_read{$x * 2 - 1})
				{
					my $err_loci = int( rand($read_length) );
					redo if $err_loci ~~ @errs1;
					push @errs1, $err_loci;
					my $ori_base = substr( $seq1, $err_loci, 1 );
					substr( $seq1, $err_loci, 1) = &simulate_seq_error($ori_base);
				}
			}
			my @errs2;
			if(exists($err_read{$x * 2}))
			{
				for (1 .. $err_read{$x * 2})
				{
					my $err_loci = int( rand($read_length) );
					redo if $err_loci ~~ @errs2;
					push @errs2, $err_loci;
					my $ori_base = substr( $seq2, $err_loci, 1 );
					substr( $seq2, $err_loci, 1) = &simulate_seq_error($ori_base);
				}
			}
			print FQ1 '@simulate:'."$seqID/1 length=$read_length\n";
			print FQ1 "$seq1\n";
			print FQ1 "+\n";
			print FQ1 ("$qualChar" x $read_length) . "\n";

			print FQ2 '@simulate:'."$seqID/2 length=$read_length\n";
			print FQ2 "$seq2\n";
			print FQ2 "+\n";
			print FQ2 ("$qualChar" x $read_length) . "\n";
			
			print OUT ">\t$x\t$seqID\n";
			if($seq_length - $start_loci < $read_length)
			{
				$isJunc = 1;
				print OUT "**\t1\t$start_loci\t$start_loci2\n";
			}
			if($seq_length - $start_loci2 < $read_length)
			{
				$isJunc = 1;
				print OUT "**\t2\t$start_loci\t$start_loci2\n";
			}
			if($isJunc == 1)
			{
				$numJuncs++;
			}
		}
		return $numJuncs;
	}
	sub simulate_seq_error
	{
		my $ori_base = $_[0];
		my @base = ('A', 'T', 'C', 'G');
		my $err_base_index;
		@base = grep { $_ ne $ori_base } @base;
		$err_base_index = int(rand(3));
		$base[$err_base_index];
	}
	sub comp_rev
	{
			my $seq = reverse($_[0]);
			$seq =~ s/[Aa]/X/g;
			$seq =~ s/[Tt]/A/g;
			$seq =~ s/X/T/g;
			$seq =~ s/[Cc]/Y/g;
			$seq =~ s/[Gg]/C/g;
			$seq =~ s/Y/G/g;
			$seq;
	}
	#{chr,geneId,transId} -> [start, end, strand];
	sub split_transcript{
	#my $gene;
		for (@_)
		{
			my @line = split /\t/;
			if($line[2] eq 'exon')
			{
				my @atr = split (/;/, $line[8]);
				my $geneId = (split(/"/, $atr[0]))[1];
				my $tranId = (split(/"/, $atr[1]))[1];
				push @{$chr_gene_trsc_exon{$line[0]}{$geneId}{$geneId}}, [ $line[3], $line[4], $line[6] ];
			}
		}
		my @line2 = split (/\t/, $_[0], 2);
		push @chr, $line2[0] unless $line2[0] ~~ @chr;
	}

	sub log10 {
		my $n = shift;
		return log($n)/log(10);
	}
