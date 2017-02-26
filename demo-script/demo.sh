
#### Demo script for each sample used

### <I>. preprocessing

#fastq-dump version: sratoolkit.2.5.7
#FastQC     version: v0.11.4
#cutadapt   version: 1.9.2.dev0
/path/to/software/fastq-dump --split-files --gzip -O /path/to/fastq/outdir /path/to/SRR445016.sra && \
mv /path/to/fastq/outdir/SRR445016_1.fastq.gz /path/to/fastq/outdir/SRR445016-raw_1.fastq.gz && \
mv /path/to/fastq/outdir/SRR445016_2.fastq.gz /path/to/fastq/outdir/SRR445016-raw_2.fastq.gz && \
/path/to/software/fastqc -t 4 /path/to/fastq/SRR445016-raw_1.fastq.gz /path/to/fastq/SRR445016-raw_2.fastq.gz -o /path/to/fastqc/outdir && \
/path/to/software/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --quality-base 33 --trim-n -q 5 -m 20 -o /path/to/fastq/outdir/SRR445016_1.fastq.gz -p /path/to/fastq/outdir/SRR445016_2.fastq.gz /path/to/fastq/SRR445016-raw_1.fastq.gz /path/to/fastq/SRR445016-raw_2.fastq.gz && echo "clean done" && \
/path/to/software/fastqc -t 4 /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz -o /path/to/fastqc/outdir && echo "preprocess done"


### <II>. circRNA detection

## 1. CIRI version: 2.0.1, URL: https://sourceforge.net/projects/ciri/
#BWA	version: 0.7.12

date && /path/to/software/bwa-0.7.12/bwa mem -t 3 -T 19 /path/to/database/hg19.fa /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz > /path/to/outdir/SRR445016.sam && echo "ciri begin" && date && \
perl /path/to/software/CIRI/CIRI_v2.0.1.pl -I /path/to/SRR445016.sam -F /path/to/database/hg19.fa -A /path/to/database/hg19/gencode.v19.annotation.gtf -G /path/to/outdir/SRR445016.log -T 3 -O /path/to/outdir/SRR445016.circRNA && date && echo "ciri done"

## 2. find_circ version: 1.0, URL: https://github.com/marvin-jens/find_circ
#Bowtie2	version: 2.1.0

date && \
/path/to/software/bowtie2-2.1.0/bowtie2 -p 3 --very-sensitive --phred33 --mm -M20 --score-min=C,-15,0 \
-x /path/to/database/hg19_bowtie/hg19 -q -1 /path/to/fastq/SRR445016_1.fastq.gz -2 /path/to/fastq/SRR445016_2.fastq.gz 2> bowtie2_SRR445016.log \
| samtools view -hbuS - | samtools sort - SRR445016 && \
samtools view -hf 4 SRR445016.bam | samtools view -bS - > unmapped_SRR445016.bam && echo "SRR445016 bowtie2 done" && \
/path/to/software/find_circ/unmapped2anchors.py unmapped_SRR445016.bam | gzip > SRR445016_anchors.qfa.gz && \
if [ ! -d SRR445016 ];then 
	mkdir SRR445016
fi
/path/to/software/bowtie2-2.1.0/bowtie2 -p 3 --reorder --mm -M20 --score-min=C,-15,0 -q -x /path/to/database/hg19_bowtie/hg19 -U SRR445016_anchors.qfa.gz | \
/path/to/software/find_circ/find_circ.py -G /path/to/database/hg19/ucsc/chrs.fa -p SRR445016_ \
-s SRR445016/sites.log > SRR445016/sites.bed 2> SRR445016/sites.reads && echo "SRR445016 find_circ done" && \
grep circ SRR445016/sites.bed | grep -v chrM | /path/to/software/find_circ/sum.py -2,3 | /path/to/software/find_circ/scorethresh.py -16 1 | /path/to/software/find_circ/scorethresh.py -15 2 | /path/to/software/find_circ/scorethresh.py -14 2 | /path/to/software/find_circ/scorethresh.py 7 2 | /path/to/software/find_circ/scorethresh.py 8,9 35 | /path/to/software/find_circ/scorethresh.py -17 100000 > SRR445016/circ_candidates.bed && echo "SRR445016 filter done" && \
date

## 3. mapsplice version: 2.2.0, URL: http://www.netlab.uky.edu/p/bioinfo/MapSplice2
gzip -dc /path/to/fastq/SRR445016_1.fastq.gz > ./SRR445016_1.fastq && \
gzip -dc /path/to/fastq/SRR445016_2.fastq.gz > ./SRR445016_2.fastq && \
echo "mapsplice begin" && date && \
python /path/to/software/MapSplice-v2.2.0/mapsplice.py -1 ./SRR445016_1.fastq -2 ./SRR445016_2.fastq -c /path/to/database/hg19/ucsc/chrs.fa -x /path/to/database/hg19_bowtie/hg19 -p 3 -o /path/to/outdir/mapsplice_SRR445016_output --min-fusion-distance 200 --gene-gtf /path/to/database/hg19/gencode.v19.annotation.gtf --bam --fusion --non-canonical-double-anchor && \
echo "mapsplice end" && date  && \
rm ./SRR445016_2.fastq ./SRR445016_1.fastq

## 4. segemehl version: 0.2.0, URL: http://www.bioinf.uni-leipzig.de/Software/segemehl/
date && \
/path/to/software/segemehl_0_2_0/segemehl/segemehl.x -s -t 3 -d /path/to/database/hg19.fa -i /path/to/database/hg19.idx -q /path/to/fastq/SRR445016_1.fastq.gz -p /path/to/fastq/SRR445016_2.fastq.gz -S | /path/to/software/samtools-0.1.19/samtools view -bS - | /path/to/software/samtools-0.1.19/samtools sort -o - /path/to/outdir/SRR445016 | samtools view -h - > /path/to/outdir/SRR445016.sam && \
/path/to/software/segemehl_0_2_0/segemehl/testrealign.x -t 3 -d /path/to/database/hg19.fa -q /path/to/outdir/SRR445016.sam -n -U /path/to/outdir/SRR445016_splicesites.bed -T /path/to/outdir/SRR445016_transrealigned.bed && \
date && echo "segemehl done" && \
rm /path/to/outdir/SRR445016.sam

## 5. NCLScan version: 1.5, URL: https://github.com/TreesLab/NCLscan
#since the maximum read length in our data is 101, we reset max_read_len = 101, and keep other parameters unchanged in the NCLScan.config file
./NCLscan.sh /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz SRR445016 20 3 50 2>&1 | tee  SRR445016_NCLscan.log

## 6. uroborus version: 0.0.2, URL: https://github.com/WGLab/UROBORUS
#tophat-2.1.1.Linux_x86_64
#bowtie-1.1.2
#bowtie2-2.2.6

date && \
/path/to/software/tophat-2.1.1.Linux_x86_64/tophat -p 3 -o SRR445016_tophat_out -G /path/to/database/hg19/genes.gtf /path/to/database/hg19_bowtie2/hg19 /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz && \
date && \
/path/to/software/samtools-0.1.19/samtools view SRR445016_tophat_out/unmapped.bam > SRR445016_tophat_out/unmapped.sam && \
cd SRR445016_tophat_out && \
perl /path/to/software/UROBORUS-0.0.2/bin/UROBORUS.pl -p 3 -index /path/to/database/hg19_bowtie/hg19 -gtf /path/to/database/hg19/genes.gtf -fasta /path/to/database/hg19/ucsc/chrs.fa unmapped.sam && \
date

## 7. KNIFE version: 1.4, URL: https://github.com/lindaszabo/KNIFE
#bowtie-1.1.2
#bowtie2-2.2.6
#perl
#python2.6
#R 3.1
#samtools 0.1.19

#!/bin/bash
/path/to/software/KNIFE/circularRNApipeline_Standalone/completeRun.sh /path/to/outdir/SRR445016/reads complete /path/to/outdir/SRR445016 SRR445016 13 phred33 circReads 66 2>&1 | tee SRR445016_out.log 

## 8. PTESFinder version: 1.0, URL: https://sourceforge.net/projects/ptesfinder-v1/
#processReadsForPTESFinder.pl: in-house script used to merge fastq_1 & fastq_2 reads from Paired-End sequencing into one file as the input of PTESFinder, provided in this directory.

#!/bin/bash
if [ ! -d SRR445016 ]; then
   mkdir SRR445016;
fi
perl /path/to/script/processReadsForPTESFinder.pl /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz > ./SRR445016.fastq && \
date && \
/path/to/software/PTESFinder/PTESFinder.sh -c /path/to/software/PTESFinder/code -r ./SRR445016.fastq -d SRR445016 -t /path/to/software/PTESFinder/data/ucsc-hg19-refGene.bed -g /path/to/software/PTESFinder/data/hg19.fa -b /path/to/software/PTESFinder/data/hg19 -j 10 -s 90 -u 7  && date && echo "ptesfinder done"

## 9. circRNA_finder version: NA, URL: https://github.com/orzechoj/circRNA_finder
#STAR	version: STAR_2.5.1a

/path/to/software/STAR --genomeDir /path/to/database/hg19_star --readFilesIn /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz --runThreadN 3 --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix /path/to/starout/SRR445016_ --outSAMtype BAM Unsorted --readFilesCommand zcat
cd /path/to/starout/

if [ ! -d circRNA_finder ];then
    mkdir circRNA_finder
fi

date && \
/path/to/software/circRNA_finder/postProcessStarAlignment.pl ./ ./circRNA_finder/ && \
date

## 10. circexplorer version: 1.1.5, URL: https://github.com/YangLab/CIRCexplorer
#We used the STAR output from circRNA_finder as circexplorer's input in our paper
/path/to/software/STAR --genomeDir /path/to/database/hg19_star --readFilesIn /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz --runThreadN 3 --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix /path/to/starout/SRR445016_ --outSAMtype BAM Unsorted --readFilesCommand zcat
cd /path/to/starout/

if [ ! -d CIRCexplorer ];then
        mkdir CIRCexplorer
fi

date && \
/path/to/software/CIRCexplorer-1.1.5/circ/star_parse.py SRR445016_Chimeric.out.junction CIRCexplorer/SRR445016_junction.txt && \
/path/to/software/CIRCexplorer-1.1.5/circ/CIRCexplorer.py -j CIRCexplorer/SRR445016_junction.txt -g /path/to/database/hg19.fa -r /path/to/software/CIRCexplorer-1.1.5/circ/ref.txt -o CIRCexplorer/SRR445016 && \
date

## 11. DCC version: 0.3.2, URL: https://github.com/dieterich-lab/DCC
#STAR   version: STAR_2.5.1a

mkdir -p /path/to/outdir/SRR445016_pairs /path/to/outdir/SRR445016_mate1 /path/to/outdir/SRR445016_mate2

/path/to/software/STAR --runThreadN 3   --genomeDir /path/to/database/hg19_star --outSAMtype BAM Unsorted --readFilesIn /path/to/fastq/SRR445016_1.fastq.gz /path/to/fastq/SRR445016_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix /path/to/outdir/SRR445016_pairs/SRR445016_pairs_ --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --outFilterMultimapNmax 20   --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2  --chimSegmentMin 15 --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15 && echo "SRR445016 pairs done" && \
/path/to/software/STAR --runThreadN 3   --genomeDir /path/to/database/hg19_star  --outSAMtype None --readFilesIn /path/to/fastq/SRR445016_1.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix /path/to/outdir/SRR445016_mate1/SRR445016_mate1_ --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15 && echo "SRR445016 mate1 done" && \
/path/to/software/STAR --runThreadN 3   --genomeDir /path/to/database/hg19_star  --outSAMtype None --readFilesIn /path/to/fastq/SRR445016_2.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix /path/to/outdir/SRR445016_mate2/SRR445016_mate2_ --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15 && echo "SRR445016 mate2 done"

python /path/to/software/DCC-0.3.2/DCC/main.py @/path/to/input/config/SRR445016.samplesheet -mt1 @/path/to/input/config/SRR445016.mate1 -mt2 @/path/to/input/config/SRR445016.mate2 -D -N -R /path/to/database/hg19/hg19_repeats.gtf -F -M -Pi -an /path/to/database/hg19/gencode.v19.annotation.gtf -Nr 1 1 -fg && echo "DCC SRR445016 done"
