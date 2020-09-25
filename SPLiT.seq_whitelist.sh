export PATH=$PATH:~/anaconda3/bin
#path to the fastq folder
fastq_folder=/home/samba/storage0/shintaku/20200814MiSeq011/cutadapt
#path to the output folder
all_output_folder=/home/samba/storage0/shintaku/20200814MiSeq011
#path to the sample list
sample_ID=$all_output_folder/sample.txt
#path to the gtf file
gtf_file=/home/samba/storage0/genome/pig/hmp_genes.gtf
#number of threads for the computation
core=16
index=/home/samba/storage0/shintaku/SPLiT-seq/star

#input_dir list_of_files output_dir barcode_pattern num_of_cores i/r(index/read)
#python /home/samba/storage0/shintaku/SPLiT-seq/fastq_index_edit.py $fastq_folder fastq_list.txt $all_output_folder CCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCC 2 r
#gunzip $fastq_folder/*.fastq.gz
#cp $fastq_folder/*R2_001.fastq.gz ./
#gunzip *R2_001.fastq.gz
# create list of files
#ls *.fastq > ./rawdatalist.txt
#demulti=./rawdatalist.txt
# fastq_pair all the pairs of fastq files
#for sample in $(cat $sample_ID); do
fastq_pair $fastq_folder/${sample}_R1_001_cutadapt.fastq $fastq_folder/${sample}_R2_001.fastq
# rename fastq files 
#mv $fastq_folder/${sample}_R1_001_cutadapt.fastq.paired.fq ${sample}_R1_paired.fastq
#mv $fastq_folder/${sample}_R2_001.fastq.paired.fq ${sample}_R2_paired.fastq
# compress fastq files
#gzip ${sample}_R1_paired.fastq
#gzip ${sample}_R2_paired.fastq; done
# remove raw demultiplexed data
#rm *.fastq
# create file list for further analysis
#ls *_R1_paired.fastq.gz > filelist.txt

for sample in $(cat $sample_ID); do umi_tools whitelist --stdin ./cutadapt/$sample*R1*gz --bc-pattern=NNNNNNNNNNCCCCCCCC --set-cell-number=5000 --log2stderr > $all_output_folder/$sample.whitelist.txt; done


