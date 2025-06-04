ln -s /home/${sample}.clean.hap1.CLEAN.fa contigsFasta

ln -s /home/HIC/${sample}_R1.fq.gz r1Reads
ln -s /home/HIC/${sample}_R2.fq.gz r2Reads

#indec
samtools faidx contigsFasta
chromap -i -r contigsFasta -o contigs.index

# alignment
chromap --preset hic -r contigsFasta -x contigs.index --remove-pcr-duplicates -1 r1Reads -2 r2Reads --SAM -o aligned.sam -t 20

#sort  
samtools view -bh aligned.sam | chromapyahs/bin/samtools sort -@ 20 -n > aligned.bam
rm aligned.sam

#step2: scaffolding
yahs contigsFasta aligned.bam

juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp contigsFasta.fai

asm_size=$(awk '{s+=$2} END{print s}' contigsFasta.fai)
java -Xmx36G -jar juicer_tools_1.19.02.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly ${asm_size}")