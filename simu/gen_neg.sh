#!/bin/bash

#ART version: 2.5.8 downloaded files: artbinmountrainier20160605linux64tgz.tgz
#illumina test examples
art=../art_illumina
#change -f 80 -> -f 200 on 20161002 for deep coverage simulate files
$art -ss HS25 -d simulate -na -i ../refMrna.fa -o ./simulate -l 101 -f 200 -p -m 350 -s 10 -sp -rs 20160830 -qs -13 -qs2 -13
gzip simulate1.fq simulate2.fq && echo "done"
mv simulate1.fq.gz simulate_1.fastq.gz && mv simulate2.fq.gz simulate_2.fastq.gz && echo "done"
#in-house script used to unify readIds
perl changeReadId.pl simulate_1.fastq.gz dsimulate_1.fastq.gz && echo "done" && \
perl changeReadId.pl simulate_2.fastq.gz dsimulate_2.fastq.gz && echo "done"
