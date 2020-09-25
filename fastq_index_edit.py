"""
Created on Tue Apr  5 22:16:54 2016

@author: Junyue
"""

import subprocess
import sys
from Levenshtein import distance
import gzip
from multiprocessing import Pool
from functools import partial

'''
    this script accept a read1 file, a read2 file, a output_file, a oligodT, barcode list, 
    then it open the read1 and read2, output file,
    then extract the barcode and UMI sequence in the read 1 file, and convert the
    barcode to the real barcode in the barcode list if it is within the mismatch rate of a real barcode,
    then it attach the barcode and UMI sequence to the read name of the read2 file
'''
def index_edit(line,barcode_index):
    decodeline=line.decode('utf-8')
    index1=decodeline.rindex(':')
    index2=decodeline.rindex('+')
    target=decodeline[:index1]
    for i in range(0,len(barcode_index)):
        target=target+decodeline[index1+1+barcode_index[i]]                                        
    return_line=target+decodeline[index2:]
    return return_line.encode()
def read_edit(line,barcode_index):
    decodeline=line.decode('utf-8')
    target=""
    for i in range(0,len(barcode_index)):
        target=target+decodeline[barcode_index[i]]
    return target.encode()+'\n'.encode()
def N_padding (NumN):
    addN=''
    for i in range(0,NumN):
        addN=addN + 'N'
    return addN
def expect_barcode(line2, barcodelist):
    import Levenshtein
    num_barcode=len(barcodelist)
    distance_array=[0]*num_barcode
    part_line=line2[:len(line2)-1]
    i=0
    for code in barcodelist:
        part_code=code[:len(line2)-1]
        distance_array[i]=Levenshtein.distance(part_line,part_code)
        i+=1
    # print(min(distance_array))
    matched_code=barcodelist[distance_array.index(min(distance_array))]
    add_code=matched_code[len(line2)-1:]
    # print(part_line)
    # print(add_code)
    return add_code
def Q_padding (NumN):
    addN=''
    for i in range(0,NumN):
        addN=addN + '!'
    return addN
def fastq_edit_list(sample, input_folder, output_folder, barcode_index, barcode, linetype):
    #import time
    #open the read1, read2, and output file
    #import re
    Read1 = input_folder + "/" + sample# + "_R1_001.fastq.gz"
    with open(barcode, "r") as f:
        barcodelist = [v.rstrip() for v in f.readlines()]
    # print(barcodelist)
    barcodelength=len(barcodelist[0])
    output_file1 = output_folder + "/" + sample# + "_R1.fastq.gz"
    output_file2= output_folder + "/short." + sample
    #mismatch_rate = int(mismatch_rate)
    f1 = gzip.open(Read1)
    f3 = gzip.open(output_file1, 'wb')
    #f4= gzip.open(output_file2,'wb')
    
    line1 = f1.readline()
    #print(line1)
    True_num=0
    False_num=0
    read_length=max(barcode_index)+2
    while (line1):

        if linetype=="i":
            first_line=index_edit(line1,barcode_index)
            read_length_flag=True
            second_line=f1.readline()
            third_line=f1.readline()
            fourth_line = f1.readline()
        elif linetype=="r":
            first_line=line1#index_edit(line1,barcode_index)
            line2=f1.readline()
            if len(line2) >= read_length:
                second_line=read_edit(line2,barcode_index)
                third_line=f1.readline()
                line4=f1.readline()
                fourth_line = read_edit(line4,barcode_index)
                read_length_flag=True
            elif len(line2)+barcodelength>=read_length:
                NumN=max(barcode_index)+2-len(line2);
                # addN=N_padding(NumN)
                #print(line2[len(line2)-1-barcodelength+NumN:])
                addN=expect_barcode(line2[len(line2)-1-barcodelength+NumN:].decode(),barcodelist)
                second_line=line2[:len(line2)-1]+addN.encode()+'\n'.encode()
                second_line=read_edit(second_line,barcode_index)
                #print(addN)
                #print(second_line)
                
                third_line=f1.readline()
                line4=f1.readline()
                addQ=Q_padding(NumN)
                fourth_line=line4[:len(line4)-1]+addQ.encode()+'\n'.encode()
                fourth_line=read_edit(fourth_line,barcode_index)
                read_length_flag=False
            else:
                NumN=max(barcode_index)+2-len(line2);
                addN=N_padding(NumN)
                #print(line2[len(line2)-1-barcodelength+NumN:])
                #addN=expect_barcode(line2[len(line2)-1-barcodelength+NumN:],barcodelist)
                second_line=line2[:len(line2)-1]+addN.encode()+'\n'.encode()
                second_line=read_edit(second_line,barcode_index)
                #print(addN)
                #print(second_line)
                
                third_line=f1.readline()
                line4=f1.readline()
                addQ=Q_padding(NumN)
                fourth_line=line4[:len(line4)-1]+addQ.encode()+'\n'.encode()
                fourth_line=read_edit(fourth_line,barcode_index)
                read_length_flag=False
                
                
#        print(first_line)
#        print(second_line)
#        print(third_line)
#        print(fourth_line)
        if read_length_flag == True:
            f3.write(first_line)
            f3.write(second_line)
            f3.write(third_line)
            f3.write(fourth_line)
            True_num += 1
        else:
            f3.write(first_line)
            f3.write(second_line)
            f3.write(third_line)
            f3.write(fourth_line)
            False_num += 1
        
#                break
#        if find == False:
#            line1 = f1.readline()
#            line1 = f1.readline()
#            line1 = f1.readline()
#            False_num += 1

        line1 = f1.readline()

    print(sample + ", id detected/edited: " + str(True_num) + "/" + str(False_num))
    f1.close()
    f3.close()
    #f4.close()

# this function accept an input folder and a output folder and then generate the output file with the index
def fastq_edit(input_folder, sampleID, output_folder, barcode_pattern, barcode_list, core_number, linetype):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    barcode file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, barcode_pattern)
    
    print(init_message)
    
    # generate the barcode list:
    #barcode_list = []
#    barcodes = open(barcode_file)
    cindex=[]
    for i in range(0,len(barcode_pattern)):
        if barcode_pattern[i]=='C':
            cindex.append(i)
    #    barcode_list.append(barcode.strip())
    #barcodes.close()
    
    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core_number))
    #print("Processing core number: ", core_number)
    func = partial(fastq_edit_list, input_folder = input_folder, output_folder=output_folder, barcode_index=cindex,barcode=barcode_list,linetype=linetype)
    #sciRNAseq_count(sample, input_folder, exons, genes, gene_end)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~bridge trimming is done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    barcode_pattern = sys.argv[4]
    core=sys.argv[5]
    linetype=sys.argv[6]
#    input_folder = "/home/samba/pihome/SeqData/2020/200918_M00426_0421_000000000-JB356analysis/Data/Intensities/BaseCalls"
#    sampleID = "/home/samba/storage0/shintaku/20200919MiSeq012/fastq_list.txt"
#    output_folder="/home/samba/storage0/shintaku/20200919MiSeq012"
#    barcode_pattern="XXXXCCCCCCCCCCCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCC"
    barcode_list_file="/home/samba/storage0/shintaku/SPLiT-seq/SPLiTbarcode.txt"
#    core=1
#    linetype="r"
    fastq_edit(input_folder, sampleID, output_folder, barcode_pattern, barcode_list_file, core, linetype)
