export PATH=$PATH:~/anaconda3/bin
#path to the fastq folder
fastq_folder=/home/samba/pihome/SeqData/2020/200814_M00426_0413_000000000-CYRK2analysis/Data/Intensities/BaseCalls
#path to the output folder
all_output_folder=/home/samba/storage0/shintaku/20200814MiSeq011rev
#path to the sample list
sample_ID=$all_output_folder/sample.txt

for sample in $(cat $sample_ID); do
echo $sample
echo $(zcat $fastq_folder/${sample}*1*fastq.gz|wc -l)/8 | bc
echo $(zcat $all_output_folder/*$sample*R1*fastq.gz|wc -l)/4 | bc
echo $(zcat $all_output_folder/UMI_attached/${sample}_umi*R2*fastq.gz|wc -l)/4 | bc
echo $(samtools view -c $all_output_folder/STAR/${sample}.Aligned.sortedByCoord.out.sam)
icnt=0
count=0
for sample in $(zcat $all_output_folder/count/$sample.counts.tsv.gz); do
icnt=$((icnt+1))
#echo $icnt
if test "$icnt" = "3"; then
#echo $sample
count=$((count+sample))
icnt=0
fi; done
echo $count; done
