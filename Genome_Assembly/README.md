## Haplotype-resolved_Assembly

### 1.AM560的组装

`hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2} --ul $ONT $CCS`

### 2.其他木薯基因组组装

```hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2}  $CCS```

### 3.contig过滤
Since the initial contig assembly contained plastid genome sequences, we downloaded previously published cassava mitochondrial (MK176513.1) and chloroplast (EU117376.1) genome sequences from NCBI (GCF_001659605.2) and annotated our assembled contigs based on sequence alignment using Minimap2. Contigs predominantly composed of plastid DNA were excluded. Satellite repeats were mainly identified using Jellyfish (v2.3.1) (https://github.com/gmarcais/Jellyfish) with 41-bp K-mers. GC content analysis was also calculated in identifying satellite contigs, as well as contigs with potential assembly algorithm bias or sequencing bias. After these filtering steps, the remaining contigs were used for downstream Hi-C scaffolding.

`python XXX_filter.py draft.contig.fa > contig.fa`

#两个hap单独挂载排序
Hic_scaffold_script1_hap1.sh
#和hap1方法一样
Hic_scaffold_script1_hap2.sh

#两个hap顺序合并用Hic重新挂载
Hic_scaffold_script1_joint.sh

#用juicebox手动调整

#手动调整之后的新的结果，将染色体根据Hifiasm的结果分到两个Haplotype
