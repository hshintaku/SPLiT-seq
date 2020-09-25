export PATH=$PATH:~/anaconda3/bin
#path to the fastq folder
fastq_folder=/home/samba/pihome/SeqData/2020/200814_M00426_0413_000000000-CYRK2analysis/Data/Intensities/BaseCalls
#path to the output folder
all_output_folder=/home/samba/storage0/shintaku/20200814MiSeq011
#path to the sample list
sample_ID=$all_output_folder/sample.txt
#path to the gtf file
gtf_file=/home/samba/storage0/genome/pig11.1/hmp_genes.gtf
#number of threads for the computation
core=16
index=/home/samba/storage0/shintaku/SPLiT-seq/star

#input_dir list_of_files output_dir barcode_pattern num_of_cores i/r(index/read)
python /home/samba/storage0/shintaku/SPLiT-seq/fastq_index_edit.py $fastq_folder fastq_list.txt $all_output_folder XXXXCCCCCCCCCCCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCC 2 r
gunzip *.fastq.gz
cp $fastq_folder/*R2_001.fastq.gz ./
gunzip *fastq.gz
# create list of files
#ls *.fastq > ./rawdatalist.txt
#demulti=./rawdatalist.txt
# fastq_pair all the pairs of fastq files
for sample in $(cat $sample_ID); do
fastq_pair ${sample}_R1_001.fastq ${sample}_R2_001.fastq
# rename fastq files 
mv ${sample}_R1_001.fastq.paired.fq ${sample}_R1_paired.fastq
mv ${sample}_R2_001.fastq.paired.fq ${sample}_R2_paired.fastq
# compress fastq files
gzip ${sample}_R2_paired.fastq
gzip ${sample}_R1_paired.fastq; done
# remove raw demultiplexed data
rm *.fastq
# create file list for further analysis
#ls *_R1_paired.fastq.gz > filelist.txt

for sample in $(cat $sample_ID); do umi_tools whitelist --stdin ./$sample*R1*gz --bc-pattern=NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCC --set-cell-number=5000 --log2stderr > $all_output_folder/$sample.whitelist.txt; done

UMI_attached_read=$all_output_folder/UMI_attached
mkdir $UMI_attached_read
for sample in $(cat $sample_ID); do umi_tools extract --bc-pattern=NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCC --stdin ./$sample*R1*paired*gz --stdout $UMI_attached_read/${sample}_umi.R1.fastq.gz --read2-in ./$sample*R2*paired*gz --read2-out $UMI_attached_read/${sample}_umi.R2.fastq.gz --filter-cell-barcode --whitelist=$all_output_folder/$sample.whitelist.txt; done

star_output=$all_output_folder/STAR
mkdir $star_output
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do
 STAR --runThreadN $core --genomeDir $index --readFilesIn $UMI_attached_read/${sample}_umi.R2.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterMismatchNmax 30 --genomeLoad LoadAndKeep --outFileNamePrefix $star_output/${sample}. --outReadsUnmapped Fastx ;done
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do samtools sort -T $star_output/${sample} -o $star_output/${sample}.Aligned.sortedByCoord.out.bam -O bam $star_output/${sample}.Aligned.out.sam;done

fC_output=$all_output_folder/featureCountDir
mkdir $fC_output

##for sample in $(cat $sample_ID); do featureCounts -a $gtf_file -o $fC_output/${sample}.gene_assigned -R BAM $star_output/${sample}.Aligned.sortedByCoord.out.bam -T $core
##done

##for sample in $(cat $sample_ID)
##do
##samtools sort $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam -o $fC_output/${sample}.assigned_sorted.bam
##samtools index $fC_output/${sample}.assigned_sorted.bam
##done

################# Generate the sparse gene count matrix
# count reads mapping to genes
input_folder=$star_output
#sample_ID=$all_output_folder/barcode_samples.txt

script=/home/samba/storage0/shintaku/ELPseq_pipeline/primary_pipeline_scripts/sciRNAseq_count3.py
echo "Start the gene count...."
python_path="/home/shintaku/anaconda3/bin"
$python_path/python $script $gtf_file $input_folder $sample_ID $core


for sample in $(cat $sample_ID)
do
samtools view -bS $fC_output/${sample}_annotated.sam > $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam
samtools sort $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam -o $fC_output/${sample}.assigned_sorted.bam
samtools index $fC_output/${sample}.assigned_sorted.bam
done

count_output=$all_output_folder/count
mkdir $count_output
for sample in $(cat $sample_ID)
do
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $fC_output/${sample}.assigned_sorted.bam -S $count_output/${sample}.counts.tsv.gz
done
