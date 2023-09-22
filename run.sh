# ssh 192.168.1.23
rootpath=/mnt/f/data/chip/liver
samples=(fetal-liver-H3K4me1-1 fetal-liver-H3K4me1-2 fetal-liver-H3K27ac-1 fetal-liver-H3K27ac-2)
# cd $rootpath

# --------------------------->
# QC and pre-processing
# --------------------------->
mkdir -p $rootpath/1.QC
# cd $rootpath/1.QC
for i in ${samples[*]}
do
  fastqc -o $rootpath/1.QC/ $rootpath/fastq/${i}_1.fastq.gz &
  fastqc -o $rootpath/1.QC/ $rootpath/fastq/${i}_2.fastq.gz &
done
wait

# --------------------------->
# mapping with bowtie2
# --------------------------->
# core 线程数
# 文件重命名策略：输出文件名称，SRR8435039:CTCF1，旧名字：新名字
mkdir -p $rootpath/2.map/bowtie2_result
mkdir -p $rootpath/2.map/bowtie2_summary
bt2ref=/home/nnlrl/genome/bovine_UCSC/bowtie2_bovine/bosTau9
cores=4
for i in ${samples[*]}
do

  bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --no-unal --phred33 \
    -I 10 -X 1000 -p $cores -x $bt2ref \
    -1 $rootpath/fastq/${i}_1.fastq.gz -2 $rootpath/fastq/${i}_2.fastq.gz \
    -S $rootpath/2.map/bowtie2_result/${i}_bowtie2.sam \
    &> $rootpath/2.map/bowtie2_summary/${i}_bowtie2.txt &
done
wait
# The paired-end reads are aligned by Bowtie2 using parameters --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 for mapping of inserts 10-700 bp in length.
# Critical step: There is no need to trim reads from out standard 25x25 PE sequencing, as adapter sequences will not be included in reads of inserts >25 bp. However, for users performing longer sequencing, reads will need to be trimmed by Cutadapt and mapped by --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 to ignore any remaining adapter sequence at the 3’ ends of reads during mapping.

# U00096.3
# 可以去ncbi下载E.coli的ref
# spikeinref=/public/workspace/shaojf/Course/NGS/Reference/bowtie2_Index/Ecoli
# for info in SRR8435039:CTCF1 SRR8435040:CTCF2 SRR8435051:IgG1 SRR8435052:IgG2
# do
#   i=`echo $info | cut -d ":" -f 1`
#   name=`echo $info | cut -d ":" -f 2`
#   bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
#     -I 10 -X 700 -p $cores -x $spikeinref \
#     -1 $rootpath/data/${i}_1.fastq.gz -2 $rootpath/data/${i}_2.fastq.gz \
#     -S $rootpath/2.map/bowtie2_result/${name}_bowtie2_spikeIn.sam \
#     &> $rootpath/2.map/bowtie2_summary/${name}_bowtie2_spikeIn.txt &
# done
# wait

# grep "rate" 2.map/bowtie2_summary/*.txt
# 计算有多少对读段比对到E.coli上,再写到一个文件里，方便后面计算scale factor
# for name in ${samples[*]}
# do
#   seqDepthDouble=`samtools view -F 0x04 $rootpath/2.map/bowtie2_result/${name}_bowtie2_spikeIn.sam | wc -l`
#   seqDepth=$((seqDepthDouble/2))
#   echo $seqDepth > $rootpath/2.map/bowtie2_result/${name}_bowtie2_spikeIn.seqDepth
# done

# In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we do not recommend removing the duplicates.
mkdir -p $rootpath/2.map/removeDuplicate/picard_summary
for name in ${samples[*]}
do
	picard SortSam I=$rootpath/2.map/bowtie2_result/${name}_bowtie2.sam \
		O=$rootpath/2.map/removeDuplicate/${name}_bowtie2.sorted.sam SORT_ORDER=coordinate
	picard MarkDuplicates I=$rootpath/2.map/removeDuplicate/${name}_bowtie2.sorted.sam \
		O=$rootpath/2.map/removeDuplicate/${name}_bowtie2.sorted.dupMarked.sam METRICS_FILE=$rootpath/2.map/removeDuplicate/picard_summary/${name}_picard.dupMark.txt
done

# -F unmap
mkdir -p $rootpath/2.map/fragmentLen
for name in ${samples[*]}
do
  samtools view -F 0x04 $rootpath/2.map/bowtie2_result/${name}_bowtie2.sam | \
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
    sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
    > $rootpath/2.map/fragmentLen/${name}_fragmentLen.txt &
done

# --------------------------------------->
# filtering and file format conversion
# --------------------------------------->
# sam转成bam，并且留下Q10以上的
# 转成bedgraph格式
# -F INT   only include reads with none of the FLAGS in INT present [0]

mkdir -p $rootpath/3.filter/
for name in ${samples[*]}
do
  samtools view -bS -F 0x04 -q 10 $rootpath/2.map/bowtie2_result/${name}_bowtie2.sam \
    >$rootpath/3.filter/${name}_bowtie2.q10.bam &
done
wait

# 转成bedpe
for name in ${samples[*]}
do
  bedtools bamtobed -i $rootpath/3.filter/${name}_bowtie2.q10.bam \
  -bedpe > $rootpath/3.filter/${name}_bowtie2.bed &
done
wait

# 上一步已转换成bedpe格式，用于表示双末端bed文件，一行里面有两个坐标，前三列一个坐标（正向测），三到六另一个（反向测）
# awk 第一列等于第四列（在同一染色体上）且终止-起始<1000
for name in ${samples[*]}
do
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $rootpath/3.filter/${name}_bowtie2.bed \
    > $rootpath/3.filter/${name}_bowtie2.clean.bed &
done
wait

# 如果第一列相同，则按第二列从小到大排序，如果第二列相同，按照第三列从小到大排
# -k，1 只看第一列
for name in ${samples[*]}
do
  cut -f 1,2,6 $rootpath/3.filter/${name}_bowtie2.clean.bed | \
  sort -k1,1 -k2,2n -k3,3n > $rootpath/3.filter/${name}_bowtie2.fragments.bed &
done
wait

# scale factor
# if [ expression ]
# then
#    Statement(s) to be executed if expression is true
# fi
# Shell expression求值。如果结果值是true，给定statement被执行。如果 expression 为false ，则没有语句将不会被执行。

# bc-l 相当于一个计算器，用长浮点型来表示
# echo 1+3 | bc 把计算的公式穿给bc去运算
# chromSize=/public/workspace/shaojf/Course/NGS/Reference/bowtie2_Index/GRCh38.chrom.size
# 长度
# for name in CTCF1 CTCF2 IgG1 IgG2
# do
#   seqDepth=`cat $rootpath/2.map/bowtie2_result/${name}_bowtie2_spikeIn.seqDepth`
#   if [[ "$seqDepth" -gt "1" ]]
#   then
#       mkdir -p $rootpath/3.filter/bedgraph
#       scale_factor=`echo "10000 / $seqDepth" | bc -l`
#       echo "Scaling factor for $name is: $scale_factor"
#       bedtools genomecov -bg -scale $scale_factor \
#       	-i $rootpath/3.filter/${name}_bowtie2.fragments.bed \
#       	-g $chromSize > $rootpath/3.filter/bedgraph/${name}_bowtie2.fragments.normalized.bedgraph &
#   fi
# done
# wait
# for name in CTCF1 CTCF2 IgG1 IgG2
# do
#   bedSort $rootpath/3.filter/bedgraph/${name}_bowtie2.fragments.normalized.bedgraph $rootpath/3.filter/bedgraph/${name}_bowtie2.fragments.normalized.srt.bedgraph &
# done
# wait

# 转bigwig方便可视化
# 用IGV打开，调颜色，bed文件第四列在IGV中会标注
mkdir -p $rootpath/3.filter/bigwig
for name in ${samples[*]}
do
  bedGraphToBigWig $rootpath/3.filter/${name}_bowtie2.fragments.bed $chromSize $rootpath/3.filter/bigwig/${name}_bowtie2.fragments.bw &
done
wait

# ------------->
# peak calling
# ------------->
# SEACR
# SEACR is intended to call peaks and enriched regions from sparse CUT&RUN or chromatin profiling data in which background is dominated by "zeroes" (i.e. regions with no read coverage).
# bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] [norm | non] [relaxed | stringent] output prefix

# 加不加0.01差很多，把top1%作为peak
seacr=SEACR_1.3.sh
mkdir -p $rootpath/4.peakCalling/
for pair in ${samples[*]}
do
  name=`echo $pair | cut -d ":" -f 1`
  control=`echo $pair | cut -d ":" -f 2`
#   bash $seacr \
#     $rootpath/3.filter/bedgraph/${name}_bowtie2.fragments.normalized.bedgraph \
#     $rootpath/3.filter/bedgraph/${control}_bowtie2.fragments.normalized.bedgraph \
#     non stringent $rootpath/4.peakCalling/${name}_seacr_control.peaks &
  bash $seacr \
    $rootpath/3.filter/${name}_bowtie2.fragments.normalized.bedgraph \
    0.05 \
    non stringent $rootpath/4.peakCalling/${name}_seacr_top0.05.peaks &
done
wait

### alternative method with MACS2
# 单末端测序：给所有的读段延长到200/300bp的fragment
# 双末端测序：直接基于bam文件，用真正的长度。需要-f BAMPE
# mkdir -p $rootpath/4.peakCalling/MACS2
# for pair in CTCF1:IgG1 CTCF2:IgG2
# do
#   name=`echo $pair | cut -d ":" -f 1`
#   control=`echo $pair | cut -d ":" -f 2`
#   macs2 callpeak -t $rootpath/3.filter/${name}_bowtie2.q10.bam \
#         -c $rootpath/3.filter/${control}_bowtie2.q10.bam \
#         -g hs -f BAMPE -n macs2_peak_q0.1.$name \
#         --outdir $rootpath/4.peakCalling/MACS2 \
#         -q 0.1 --keep-dup all \
#         2>$rootpath/4.peakCalling/MACS2/macs2Peak_summary.$name.txt &
# done
###

# Rscript cuttag.visualization.R /public/workspace/stu18230130/MyCourse/NGS-02/lab1-CUTTag CTCF,IgG 1,2 CTCF

###
# https://github.com/Boyle-Lab/Blacklist/tree/master/lists

# ------------->
# motif
# ------------->
# MEME输入需要fasta格式
# 其他功能：想知道基因组上那些位置有某个motif,去基因组上scan
# mkdir -p $rootpath/5.motif/
# motifdb=/public/workspace/shaojf/Course/NGS/Reference/MEME/JASPAR2018_CORE_vertebrates_non-redundant.meme
# reffasta=/public/workspace/shaojf/Course/NGS/Reference/bowtie2_Index/GRCh38.fa
# blacklist=/public/workspace/shaojf/Course/NGS/Reference/ENCODE/blacklist.v2/hg38-blacklist.v2.bed

# - 表示从标准输入中获取
# sed 's?chr??' 把chr去掉
# bedtools intersect -u \
#   -a $rootpath/4.peakCalling/CTCF2_seacr_top0.05.peaks.stringent.bed \
#   -b $rootpath/4.peakCalling/CTCF1_seacr_top0.05.peaks.stringent.bed | \
#   bedtools intersect -v -a - -b <(sed 's?chr??' $blacklist) > $rootpath/4.peakCalling/CTCF.overlapped.bed

# 给坐标提取fasta，把两个重复的peak去交集
# bedtools getfasta -fi $reffasta -fo $rootpath/5.motif/CTCF.overlapped.fa -bed $rootpath/4.peakCalling/CTCF.overlapped.bed

# -oc 强行用某个文件夹 -o 需要新建的文件夹
# -db denovo找出来的和数据库比对 $motifdb已知的数据库
# meme-chip -o $rootpath/5.motif/meme.CTCF.overlapped -db $motifdb $rootpath/5.motif/CTCF.overlapped.fa
# for name in CTCF1 CTCF2
# do
# 	bedtools getfasta -fi $reffasta \
# 		-fo $rootpath/5.motif/${name}_seacr_top0.01.peaks.stringent.fa \
# 		-bed $rootpath/4.peakCalling/${name}_seacr_top0.01.peaks.stringent.bed &
# done
# wait
# for name in CTCF1 CTCF2
# do
# 	meme-chip -o $rootpath/5.motif/meme.$name -db $motifdb $rootpath/5.motif/${name}_seacr_top0.01.peaks.stringent.fa &
# done
# wait
