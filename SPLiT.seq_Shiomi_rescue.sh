export PATH=$PATH:/home/shintaku/anaconda3/bin
#path to the fastq folder
fastq_folder=/home/samba/pihome/SeqData/2020/200801_M00426_0407_000000000-CVF6Panalysis/Data/Intensities/BaseCalls
#path to the output folder
all_output_folder=/home/samba/storage0/shintaku/20200802MiSeq010
#path to the sample list
sample_ID=$all_output_folder/sample.txt
#path to the gtf file
gtf_file=/home/samba/storage0/genome/pig11.1/hmp_genes2.gtf
#gtf_file=/home/samba/storage0/shintaku/refdata-cellranger-GRCh38-and-mm10-3.1.0/genes/genes.gtf
#number of threads for the computation
core=16
index=/home/samba/storage0/shintaku/SPLiT-seq/star

#input_dir list_of_files output_dir barcode_pattern num_of_cores i/r(index/read)

#/home/shintaku/anaconda3/bin/python /home/samba/storage0/shintaku/SPLiT-seq/fastq_index_edit.py $fastq_folder $all_output#_folder/fastq_list.txt $all_output_folder XXXXCCCCCCCCCCCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCC 2 r

#cp $fastq_folder/*R2_001.fastq.gz $all_output_folder
#gunzip $all_output_folder/*fastq.gz
# fastq_pair all the pairs of fastq files
#for sample in $(cat $sample_ID); do
#fastq_pair $all_output_folder/${sample}_R1_001.fastq $all_output_folder/${sample}_R2_001.fastq
# rename fastq files 
#mv $all_output_folder/${sample}_R1_001.fastq.paired.fq $all_output_folder/${sample}_R1_paired.fastq
#mv $all_output_folder/${sample}_R2_001.fastq.paired.fq $all_output_folder/${sample}_R2_paired.fastq
# compress fastq files
#gzip $all_output_folder/${sample}_R2_paired.fastq
#gzip $all_output_folder/${sample}_R1_paired.fastq; done


for sample in $(cat $sample_ID); do umi_tools whitelist --stdin $all_output_folder/$sample*R1*gz --set-cell-number 5000 --bc-pattern=NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCC --log2stderr > $all_output_folder/$sample.whitelist.txt; done

UMI_attached_read=$all_output_folder/UMI_attached
mkdir $UMI_attached_read
for sample in $(cat $sample_ID); do umi_tools extract --bc-pattern=NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCC --stdin $all_output_folder/$sample*R1*paired*gz --stdout $UMI_attached_read/${sample}_umi.R1.fastq.gz --read2-in $all_output_folder/$sample*R2*paired*gz --read2-out $UMI_attached_read/${sample}_umi.R2.fastq.gz --filter-cell-barcode --whitelist=$all_output_folder/$sample.whitelist.txt; done

star_output=$all_output_folder/STAR
mkdir $star_output
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do
 STAR --runThreadN $core --genomeDir $index --readFilesIn $UMI_attached_read/${sample}_umi.R2.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 1 --outFilterMismatchNmax 30 --genomeLoad LoadAndKeep --outFileNamePrefix $star_output/${sample}. --outReadsUnmapped Fastx ;done
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do samtools sort -T $star_output/${sample} -o $star_output/${sample}.Aligned.sortedByCoord.out.sam -O sam $star_output/${sample}.Aligned.out.sam;done

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
samtools view -bS $star_output/${sample}_annotated.sam > $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam
samtools sort $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam -o $fC_output/${sample}.assigned_sorted.bam
samtools index $fC_output/${sample}.assigned_sorted.bam
done

count_output=$all_output_folder/count
mkdir $count_output
for sample in $(cat $sample_ID)
do
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $fC_output/${sample}.assigned_sorted.bam -S $count_output/${sample}.counts.tsv.gz
done
