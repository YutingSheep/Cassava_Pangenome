基因组质量评估——BUSCO
busco -i $GENOME/${i}_haplotype.fa -c $threads -o busco -m geno -l embryophyta_odb10 --offline

基因组质量评估——completeness & QV
meryl k=19 count output ${sample}.sr1.meryl /home/WGS/${sample}_1_clean.fq.gz
meryl k=19 count output ${sample}.sr2.meryl /home/WGS/${sample}_2_clean.fq.gz

meryl union-sum output ${sample}.sr.k19.meryl ${sample}.sr*.meryl

merqury.sh ${sample}.sr.k19.meryl ${sample}_hap1.fa ${sample}_hap2.fa out_prefix

基因组质量评估——switch error
## step1 HiFi reads构建标准数据集
mkdir 01.hifi-mapping
cd 01.hifi-mapping
cp /home/haplotype/${sample}_hap1.fa ${sample}.fasta
samtools faidx ${sample}.fasta
minimap2 -t 5 -ax map-pb --secondary=no ${sample}.fasta ${sample}.ccs.fq > ${sample}.sam
samtools view -bt ${sample}.fasta.fai ${sample}.sam > ${sample}.bam
samtools sort -@ 5 -o ${sample}.pb.sorted.bam ${sample}.bam
samtools index ${sample}.pb.sorted.bam
cd ..

## step2 两个单倍型比较得到实际组装phased 状态
mkdir 02.hap-snp-bcftools
cd 02.hap-snp-bcftools
minimap2 -a -x asm20 --cs -r2k -t 10 /home/haplotype/${sample}_hap1.fa /home/haplotype/${sample}_hap2.fa > ${sample}.hap1hap2.sam
samtools sort -O BAM -o ${sample}.hap1hap2.bam ${sample}.hap1hap2.sam
samtools index ${sample}.hap1hap2.bam
bcftools mpileup -f /home/haplotype/${sample}_hap1.fa --threads 10 -o ${sample}.hap1hap2.vcf -A -O v ${sample}.hap1hap2.bam
bcftools call -v -c ${sample}.hap1hap2.vcf -o ${sample}.hap1hap2.called.vcf
sed -i 's/1\/1/0|1/g' ${sample}.hap1hap2.called.vcf
sed -i 's/${sample}.hap1hap2.bam/${sample}/g' ${sample}.hap1hap2.called.vcf
cd ..

## step3 实际的phase和标准phase进行比较
mkdir 03.compare
cd 03.compare
#chr为染色体list
cp /home/.../chr .
cat chr| parallel -j 10 'whatshap phase --ignore-read-groups -o {}.phased.vcf --reference=../01.hifi-mapping/${sample}.fasta --chromosome ${sample}HA{} ../01.hifi-mapping/output.vcf.gz ../01.hifi-mapping/${sample}.pb.sorted.bam'
cat chr*.phased.vcf|grep -v '#'|grep PS > pb.wh.phase.vcf
cat chr01.phased.vcf | head -n 46 > ${sample}.hifi.phased.vcf
cat pb.wh.phase.vcf >> ${sample}.hifi.phased.vcf
sed -i 's/default/${sample}/g' ${sample}.hifi.phased.vcf
whatshap compare --names truth,whatshap --tsv-pairwise eval.tsv ${sample}.hifi.phased.vcf ../02.hap-snp-bcftools/${sample}.hap1hap2.called.vcf