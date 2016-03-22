#Process raw ribo-seq data: remove barcode; remove polyA tail by allowing one mismatch; and retain reads between 25nt to 35nt
#Usage: python /home/jwan/workspace/code/preprocess_bacoded_GTI.py --file=5927_5504_26675_C7PMNANXX_10346359_R1.fastq.gz --trim3=7 --pa=True --prefix=10346359

import string
from optparse import OptionParser
import sys
import gzip 

def process_command_line(argv):
   parser = OptionParser()
   parser.add_option("-f", "--file", dest="input_filename", help="Bacoded GTI file") #input file
   parser.add_option("--trim3", dest="trim3_length", help="Trim a fixed length of sequence at the 3' end", default=0) # how many reads are trimmed before searching for A
   parser.add_option("--pa",dest="trim_polyA", help="Trim polyA tail at the 3' end", default=False) #whether trim polyA or not; most of the time set it to True
   parser.add_option("--prefix", dest="output_prefix", help = "prefix of output files") #output prefix file; 5 files will be generated corresponding to 5 barcodes

   (options, args) = parser.parse_args()

   return options, args

if __name__ == "__main__":
   

   (options, args) = process_command_line(sys.argv)
   valid_barcode = ["TG", "AC", "GA","CT","AG"]
   outfile_prefix = options.output_prefix
   file_handle1 = open(outfile_prefix+'_'+'TG','w')
   file_handle2 = open(outfile_prefix+'_'+'AC','w')
   file_handle3 = open(outfile_prefix+'_'+'GA','w')
   file_handle4 = open(outfile_prefix+'_'+'CT','w')
   file_handle5 = open(outfile_prefix+'_'+'AG','w')

   file_output = open(outfile_prefix+'_'+'RL','w')

   outfile_handle_list = [file_handle1, file_handle2, file_handle3, file_handle4,file_handle5]
   flag = 0
   read_length_dic = {}

   for line in gzip.open(options.input_filename,'rb'):
      if flag == 0 and line[0] == "@":
         read_name = line
         flag  = 1
      elif flag == 1:
         read_seq = string.strip(line)
         flag = 2
      elif flag == 2 and line[0] == '+':
         flag = 3
      elif flag  == 3:
         read_quality_seq = line
         flag = 4   
   
      if flag == 4:
         flag = 0
         barcode = read_seq[:2]
         
         ##trim read by trim3_length nucleotide before removing polyA tail from the 3' end

         if int(options.trim3_length) > 0 and options.trim_polyA == "True":
            trimmed_length = int(options.trim3_length)
            read_length = len(read_seq) - 1 #\n is also counted
            num_mismatch = 0
            
            ##recursively find the last position which is not A by allowing one mismatch

            index_3end_valid_read = 0
            for i in range(read_length-trimmed_length-1,-1,-1):
               if read_seq[i] == "A":
                  continue
               elif read_seq[i] != "A":
                  num_mismatch += 1
                  if num_mismatch > 1:
                     index_3end_valid_read = i #inclusive
                     break
            
            if 25 < index_3end_valid_read < 37: # if after trimming, the read length is between 25 and 35nt, it is retained for further process

               try:
                  barcode_index = valid_barcode.index(barcode)
               except ValueError:
                  continue

               read_length_dic.setdefault(barcode,{})
               read_length_dic[barcode].setdefault(index_3end_valid_read-1,0)
               read_length_dic[barcode][index_3end_valid_read-1] += 1 

               outfile_handle_list[barcode_index].write(read_name)
               outfile_handle_list[barcode_index].write(read_seq[2:(index_3end_valid_read+1)]+'\n')
               outfile_handle_list[barcode_index].write('+'+'\n')
               outfile_handle_list[barcode_index].write(read_quality_seq[2:(index_3end_valid_read+1)]+'\n')
         
         elif int(options.trim3_length) == 0 and options.trim_polyA == "False":
            read_length = len(read_seq) - 1 #\n is also counted
            try:
               barcode_index = valid_barcode.index(barcode)
            except ValueError:
               continue

            read_length_dic.setdefault(barcode,{})
            read_length_dic[barcode].setdefault(read_length-2,0)
            read_length_dic[barcode][read_length-2] += 1 

            outfile_handle_list[barcode_index].write(read_name)
            outfile_handle_list[barcode_index].write(read_seq[2:read_length]+'\n')
            outfile_handle_list[barcode_index].write('+'+'\n')
            outfile_handle_list[barcode_index].write(read_quality_seq[2:read_length]+'\n')

   
   for barcode in read_length_dic:         
      for read_length in read_length_dic[barcode]:
         file_output.write(str(barcode)+'\t'+str(read_length)+'\t'+str(read_length_dic[barcode][read_length])+'\n')    