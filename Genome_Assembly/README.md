# Non-AM560 sample
$hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2}  $CCS

# AM560
$hifiasm -t $threads -o TEST1-MS --primary -l3 --h1 ${Hic_1} --h2 ${Hic_2} --ul $ONT  $CCS

# 过滤小于10kb的从提高和质体序列
python XXX_filter.py draft.contig.fa > contig.fa

#两个hap单独挂载排序
Hic_scaffold_script1_hap1.sh
#和hap1方法一样
Hic_scaffold_script1_hap2.sh

#两个hap顺序合并用Hic重新挂载
Hic_scaffold_script1_joint.sh

#用juicebox手动调整

#手动调整之后的新的结果，将染色体根据Hifiasm的结果分到两个Haplotype
