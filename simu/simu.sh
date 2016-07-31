date && \
perl ./CIRIsimulator.pl -1 pos_1.fastq -2 pos_2.fastq -O pos_circRNA.list -G ./refGene.gtf -DB ./hela_hg19_circRNA.txt -C 10 -LC 0 -R 1 -LR 1 -L 101 -E 1 -I 350 -D ./database/hg19/ucsc -CHR1 0 -M 50 && \
gzip pos_1.fastq pos_2.fastq && \
date
