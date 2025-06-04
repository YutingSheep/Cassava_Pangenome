ln -s /home/HIC/${sample}_R1.fq.gz r1Reads
ln -s /home/HIC/${sample}_R2.fq.gz r2Reads
cp /home/haplotypes/${sample}.hap1.fa contigsFasta
cat /home/haplotypes/${sample}.hap2.fa >> contigsFasta

cp  ../sample-hap1/yahs.out_scaffolds_final.agp scaffold.joint.agp
cat ../sample-hap2/yahs.out_scaffolds_final.agp >> scaffold.joint.agp

#Align Hi-C data to the assembly, remove PCR duplicates and filter out secondary and supplementary alignments
bwa index contigsFasta
bwa mem -5SP contigsFasta r1Reads r2Reads -t 9 | samblaster | samtools view - -@ 9 -S 
-h -b -F 3340 -o HiC.bam

rm -rf HiC.filtered.bam
#Filter the alignments with MAPQ 0 (mapping quality =1) and NM 3 (edit distance < 3)
/home/software/HapHiC/utils/filter_bam HiC.bam 1 --nm 3 --threads 9 |samtools view - -b -@ 9 -o HiC.filtered.bam

samtools faidx contigsFasta
/home/software/HapHiC/scripts/../utils/juicer pre -a -q 1 -o out_JBAT HiC.filtered.bam scaffold.joint.agp contigsFasta.fai >out_JBAT.log 2>&1