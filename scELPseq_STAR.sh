export PATH=$PATH:~/anaconda3/bin
#path to the fastq folder
#fastq_folder=/home/samba/pihome/SeqData/20190914_MiSeq007/Analysis2
fastq_folder=/home/samba/storage0/Shiomi/data_0609
#path to the output folder
all_output_folder=/home/samba/storage0/shintaku/20200609MiSeq009
#path to the sample list
sample_ID=$all_output_folder/sampleidall.txt
#path to the gtf file
gtf_file=/home/samba/storage0/shintaku/refdata-cellranger-GRCh38-and-mm10-3.1.0/genes/genes.gtf
#number of threads for the computation
core=16
index=/home/samba/storage0/shintaku/homo_sapiens/star_db

#for sample in $(cat $sample_ID); do umi_tools whitelist --stdin $fastq_folder/$sample*R1*gz --bc-pattern=CCCCCCCCCCNNNNNNNN --set-cell-number=16 --log2stderr > $all_output_folder/$sample.whitelist.txt; done

#UMI_attached_read=$all_output_folder/UMI_attached
#mkdir $UMI_attached_read
#for sample in $(cat $sample_ID); do umi_tools extract --bc-pattern=CCCCCCCCCCNNNNNNNN --stdin $fastq_folder/$sample*R1*gz --stdout $UMI_attached_read/${sample}_umi.R1.fastq.gz --read2-in $fastq_folder/$sample*R2*gz --read2-out $UMI_attached_read/${sample}_umi.R2.fastq.gz --filter-cell-barcode --whitelist=$all_output_folder/$sample.whitelist.txt; done

star_output=$all_output_folder/STAR
#mkdir $star_output
STAR --genomeDir $index --genomeLoad Remove
for sample in $(cat $sample_ID); do
 STAR --runThreadN $core --genomeDir $index --readFilesIn ./unmatched_paired_R2.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 20 --outFilterMismatchNmax 30 --genomeLoad LoadAndKeep --outFileNamePrefix $star_output/${sample}. --outReadsUnmapped Fastx ;done
STAR --genomeDir $index --genomeLoad Remove
#for sample in $(cat $sample_ID); do samtools sort -T $star_output/${sample} -o $star_output/${sample}.Aligned.sortedByCoord.out.bam -O bam $star_output/${sample}.Aligned.out.sam;done

#fC_output=$all_output_folder/featureCountDir
#mkdir $fC_output

#for sample in $(cat $sample_ID); do featureCounts -a $gtf_file -o $fC_output/${sample}.gene_assigned -R BAM $star_output/${sample}.Aligned.sortedByCoord.out.bam -T $core
#done

#for sample in $(cat $sample_ID)
#do
#samtools sort $fC_output/${sample}.Aligned.sortedByCoord.out.bam.featureCounts.bam -o $fC_output/${sample}.assigned_sorted.bam
#samtools index $fC_output/${sample}.assigned_sorted.bam
#done

#count_output=$all_output_folder/count
#mkdir $count_output
#for sample in $(cat $sample_ID)
#do
#umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I $fC_output/${sample}.assigned_sorted.bam -S $count_output/${sample}.counts.tsv.gz
#done
