#1. Obtain read density information for each gene
#2. Consider all transcripts instead of longest one
#python process_vitro_GTI_all_transcript_sam.py sam_file annotation_file psite_mode normalization_mode directory_to_write_output_files cutoff_total_psite_occupancy
#python process_vitro_GTI_all_transcript_sam.py example.sam refseq_annotation_hg19 middle original output_example

import string
import sys
import os
import re

#Convert mapped reads to psite coor
def read_tophat_sam(line,psite_mode):
   pattern = re.compile(r"^(\d+)M(\d+)N(\d+)M$")
   line_ele = string.split(string.strip(line))
   flag_strand = int(line_ele[1])
   if (flag_strand/16)%2==0:
      strand = "+"
   else:
      strand = "-"

   if string.find(line_ele[2], 'chr') == -1:
      chr = 'chr'+line_ele[2]
   else:
      chr = line_ele[2]

   coor_start = int(line_ele[3])
   match_info = line_ele[5]
   if string.find(match_info,'I') != -1:
      return (chr,strand,0)
   if string.find(match_info,'D') != -1:
      return (chr,strand,0)

   if string.find(match_info,'N') != -1:
      try:
         match_regions = [int(ele) for ele in pattern.search(match_info).groups()]
      except AttributeError:
         return (chr,strand,0)

      if psite_mode == "middle":
         if strand == "+":
            if match_regions[0] >= 13:
               psite_coor = coor_start + 12
            else:
               psite_coor = coor_start+match_regions[1]+12
         elif strand == "-":
            if match_regions[2] >= 13:
               psite_coor = coor_start+match_regions[0]+match_regions[1]+match_regions[2]-13
            else:
               psite_coor = coor_start+match_regions[0]+match_regions[2]-13
      elif psite_mode == "5prime":
         if strand == "+":
               psite_coor = coor_start
         elif strand == "-":
               psite_coor = coor_start+match_regions[0]+match_regions[1]+match_regions[2]-1
   else:
      if match_info[-1] == "M":
         match_region = int(match_info[:-1])
         if psite_mode == "middle":
            if strand == "+":
               psite_coor = coor_start + 12
            else:
               psite_coor = coor_start+match_region-13
         elif psite_mode == "5prime":
            if strand == "+":
               psite_coor = coor_start
            else:
               psite_coor = coor_start+match_region-1
         
   return (chr,strand,psite_coor)

#read ensembl GTF file
def read_ensembl_by_nt(ensembl_file):
   tx_dic = {}
   gene_dic = {}
   start_codon_dic = {}
   stop_codon_dic = {}

   for line in open(ensembl_file):
      line_ele = string.split(string.strip(line),'\t')

      if string.find(line_ele[0], 'chr') == -1:
         chr = 'chr'+line_ele[0]
      else:
         chr = line_ele[0]   

      #if chr != "chrMT":
      #   continue
      
      tx_type = line_ele[1]
      #if tx_type != "rRNA":
      #   continue

      if tx_type != "lincRNA":
         continue

      #if string.find(tx_type,'RNA') == -1:
      #   continue
      
      exon_type = line_ele[2]
      exon_start_coor = int(line_ele[3])
      exon_end_coor = int(line_ele[4])
      strand = line_ele[6]
      exon_info = [string.split(ele) for ele in string.split(line_ele[8][:-1],';')]
      exon_info_dic = {}
      
      try:
         for ele in exon_info:
            exon_info_dic.setdefault(ele[0],ele[1].replace("\"",""))
      except IndexError:
         continue

      tx_id = exon_info_dic["transcript_id"]
      gene_id = exon_info_dic["gene_id"]
      gene_name = exon_info_dic["gene_name"]
      
      print gene_id

      if exon_type == "start_codon":
         if strand == "+":
            start_codon_coor = exon_start_coor
         elif strand == "-":
            start_codon_coor = exon_end_coor

         start_codon_dic.setdefault(gene_name,{})
         start_codon_dic[gene_name].setdefault(tx_id,[])
         if (chr,strand,start_codon_coor) not in start_codon_dic[gene_name][tx_id]:
               start_codon_dic[gene_name][tx_id].append((chr,strand,start_codon_coor))

      elif exon_type == "stop_codon":
         if strand == "+":
            stop_codon_coor = exon_end_coor
         elif strand == "-":
            stop_codon_coor = exon_start_coor

         stop_codon_dic.setdefault(gene_name,{})
         stop_codon_dic[gene_name].setdefault(tx_id,[])
         if (chr,strand,stop_codon_coor) not in stop_codon_dic[gene_name][tx_id]:
            stop_codon_dic[gene_name][tx_id].append((chr,strand,stop_codon_coor))         
      
      elif exon_type == "exon":
         bin_id = exon_start_coor/100000
         tx_dic.setdefault(chr,{})
         tx_dic[chr].setdefault(strand,{})
         tx_dic[chr][strand].setdefault(bin_id,{})
         if bin_id < 10:
            bin_range = range(0,bin_id+10)
         else:
            bin_range = range(bin_id-10,bin_id+10)
         flag = 0
         for i in bin_range:
            if tx_dic[chr][strand].has_key(i):
               if tx_dic[chr][strand][i].has_key(gene_name): 
                  tx_dic[chr][strand][i].setdefault(gene_name,{})
                  tx_dic[chr][strand][i][gene_name].setdefault(tx_id,{})
                  for coor in range(exon_start_coor,exon_end_coor+1):
                     tx_dic[chr][strand][i][gene_name][tx_id].setdefault(coor,1)
                  flag = 1
                  break
         if flag == 0:
            tx_dic[chr][strand][bin_id].setdefault(gene_name,{})
            tx_dic[chr][strand][bin_id][gene_name].setdefault(tx_id,{})               
            for coor in range(exon_start_coor,exon_end_coor+1):
               tx_dic[chr][strand][bin_id][gene_name][tx_id].setdefault(coor,1)   
            tx_dic[chr][strand][bin_id][gene_name][tx_id].setdefault(-1,1)

 
         gene_dic.setdefault(chr,{})
         gene_dic[chr].setdefault(strand,{})
         gene_dic[chr][strand].setdefault(bin_id,{})
         
         flag = 0
         for i in bin_range:
            if gene_dic[chr][strand].has_key(i):
               if gene_dic[chr][strand][i].has_key(gene_name):
                  gene_dic[chr][strand][i].setdefault(gene_name,{})
                  for coor in range(exon_start_coor,exon_end_coor+1):
                      gene_dic[chr][strand][i][gene_name].setdefault(coor,1)
                  flag = 1
                  break
         if flag == 0:    
            gene_dic[chr][strand][bin_id].setdefault(gene_name,{})
            for coor in range(exon_start_coor,exon_end_coor+1):
               gene_dic[chr][strand][bin_id][gene_name].setdefault(coor,1) 
               
   return (tx_dic,gene_dic,start_codon_dic,stop_codon_dic)           

#read UCSC GenePred table format
             
def read_refseq_by_nt(refseq_file):
   tx_dic = {}
   gene_dic = {}
   start_codon_dic = {}
   stop_codon_dic = {}
   
   for line in open(refseq_file):
      if line[0] == "#":
         continue
      else:
         line_ele = string.split(string.strip(line),'\t')

         if not len(line_ele) == 16 and not len(line_ele) == 14:
            continue
 
         if line_ele[6] == line_ele[7]: #filterout non-coding RNA (without CDS)
            continue

         #if line_ele[6] != line_ele[7]:    #non-coding RNA without coding capacity
         #   continue

         tx_name = line_ele[1]
         gene_name = line_ele[12]

         if string.find(line_ele[2],'chr') != -1:
            chr = line_ele[2]
         else:
            chr = "chr"+line_ele[2]

         #if len(string.split(chr,'_')) > 1:
         #   continue

         strand = line_ele[3]
         exon_start_list = [int(ele) for ele in string.split(line_ele[9][:-1],',')]
         exon_end_list = [int(ele) for ele in string.split(line_ele[10][:-1],',')]
         exon_coordinate_list = zip(exon_start_list, exon_end_list)

         CDS_length = 0
         for exon_info in exon_coordinate_list:
            exon_length = exon_info[1] - exon_info[0]
            CDS_length += exon_length

         bin_id = int(line_ele[4])/100000
         tx_dic.setdefault(chr,{})
         tx_dic[chr].setdefault(strand,{})
         tx_dic[chr][strand].setdefault(bin_id,{}) 
         tx_dic[chr][strand][bin_id].setdefault(gene_name,{})
         tx_dic[chr][strand][bin_id][gene_name].setdefault(tx_name,{})
         
         gene_dic.setdefault(chr,{})
         gene_dic[chr].setdefault(strand,{})
         gene_dic[chr][strand].setdefault(bin_id,{}) 
         gene_dic[chr][strand][bin_id].setdefault(gene_name,{})
          
         num_exon = len(exon_coordinate_list)
         if num_exon == 1:
            exon_info = exon_coordinate_list[0]
            exon_start = int(exon_info[0]) +1    #- 99
            exon_end = int(exon_info[1])    #+ 100
            #exon_start = int(exon_info[0]) + 1 - 500   #- 99
            #exon_end = int(exon_info[1]) + 500   #+ 100
            for i in range(exon_start, exon_end+1):
               tx_dic[chr][strand][bin_id][gene_name][tx_name].setdefault(i,0)
               tx_dic[chr][strand][bin_id][gene_name][tx_name][i] += 1
               gene_dic[chr][strand][bin_id][gene_name].setdefault(i,0)
               gene_dic[chr][strand][bin_id][gene_name][i] += 1 
         else:
            for i in range(num_exon):
               exon_info = exon_coordinate_list[i]
               exon_start = int(exon_info[0]) +1 #- 99
               exon_end = int(exon_info[1])
               #if i == 0:
               #   exon_start = exon_start - 500 #- 99
               #if i == (num_exon-1):
               #   exon_end = exon_end + 500
                  
               for i in range(exon_start, exon_end+1):
                  tx_dic[chr][strand][bin_id][gene_name][tx_name].setdefault(i,0)
                  tx_dic[chr][strand][bin_id][gene_name][tx_name][i] += 1

                  gene_dic[chr][strand][bin_id][gene_name].setdefault(i,0)
                  gene_dic[chr][strand][bin_id][gene_name][i] += 1 
                 
         tx_dic[chr][strand][bin_id][gene_name][tx_name][-1] = CDS_length

         try:
            if strand == "+":
               start_codon_coor = int(line_ele[6]) + 1 #1-base
               stop_codon_coor = int(line_ele[7])
               #start_codon_seq = string.upper(str(hg19[chr][(start_codon_coor-1):(start_codon_coor + 2)])) #pygr uses 0-base
            elif strand == "-":
               start_codon_coor = int(line_ele[7])
               stop_codon_coor = int(line_ele[6]) + 1
               #start_codon_seq = string.upper(str(-hg19[chr][(start_codon_coor-3):(start_codon_coor)]))

            start_codon_dic.setdefault(gene_name,{})
            start_codon_dic[gene_name].setdefault(tx_name,[])
            if (chr,strand,start_codon_coor) not in start_codon_dic[gene_name][tx_name]:
               start_codon_dic[gene_name][tx_name].append((chr,strand,start_codon_coor))

            stop_codon_dic.setdefault(gene_name,{})
            stop_codon_dic[gene_name].setdefault(tx_name,[])
            if (chr,strand,stop_codon_coor) not in stop_codon_dic[gene_name][tx_name]:
               stop_codon_dic[gene_name][tx_name].append((chr,strand,stop_codon_coor))
         except IndexError:
            continue
   return (tx_dic,gene_dic,start_codon_dic,stop_codon_dic)           

#def calculate_distance_psite_start_codon(annotation_dic,start_codon_dic,psite_coor,gene,tx_name,chr,strand): #for single psite
#   
#   start_codon_coor = start_codon_dic[gene][tx_name][0][2]
#   for bin_id in annotation_dic[chr][strand]:
#      if annotation_dic[chr][strand][bin_id].has_key(gene):
#       if annotation_dic[chr][strand][bin_id][gene].has_key(tx_name):
#         sorted_CDS_region = sorted(annotation_dic[chr][strand][bin_id][gene][tx_name].keys())
#         distance = 0
#         try:
#            start_codon_index = sorted_CDS_region.index(start_codon_coor)
#            psite_codon_index = sorted_CDS_region.index(psite_coor)
#         except ValueError:
#            return "None"
#         if strand == "+":
#            distance = psite_codon_index - start_codon_index
#         else:
#            distance = start_codon_index - psite_codon_index
#         #print chr,strand,psite_coor,start_codon_coor,distance
#         return distance

def calculate_distance_psite_start_codon(start_codon_dic,gene,tx_name,chr,strand,sorted_CDS_coor_list): #for single psite   start_codon_dic,gene,tx_name,chr,strand,sorted_CDS_coor_list
   dic = {}
   start_codon_coor = start_codon_dic[gene][tx_name][0][2]
   for psite_coor in sorted_CDS_coor_list:
         try:
            start_codon_index = sorted_CDS_coor_list.index(start_codon_coor)
            psite_codon_index = sorted_CDS_coor_list.index(psite_coor)
         except ValueError:
            return "None"
         if strand == "+":
            distance = psite_codon_index - start_codon_index
         else:
            distance = start_codon_index - psite_codon_index
         dic.setdefault(psite_coor,distance)
   
   return dic

if __name__ == "__main__":
   (tx_annotation_dic,gene_annotation_dic,start_codon_dic,stop_codon_dic) = read_refseq_by_nt(sys.argv[2]) #read gene annotationa as well as start/stop codon information from gene annotation file
   #(tx_annotation_dic,gene_annotation_dic,start_codon_dic,stop_codon_dic) = read_ensembl_by_nt(sys.argv[2])
   
#Initiation as well as parameter setting
   GTI_total_by_gene = {}
   number_gene_dic_by_relative_position_to_TIS = {}
   normalized_count_by_dis = {}
   normalization_mode = sys.argv[4] #original_read or norm_read_total_read for original or normalized value
   psite_mode = sys.argv[3] #5prime or middle, 5prime for the position at the 5' end of RPF and middle is the 13th position of the RPF
   total_mapped_read = 0

#Read sam file and quantify each merged gene as well as mRNA isoforms
   for line in open(sys.argv[1]):
      if line[0] == "@":
         continue
      (chr, strand, psite_coor) = read_tophat_sam(line,psite_mode) #psite_coor is the converted psite coordinate from RPF
      #if chr != "chrMT":
      #   continue

      if not gene_annotation_dic.has_key(chr):
         continue
      if not gene_annotation_dic[chr].has_key(strand):
         continue
      
      bin_id = psite_coor/100000 #Assign psite to a chromosome bin for speeding up the program 
      bin_list = range(bin_id-5,bin_id+6)

#Search gene containing current psite from gene annotation and calculate the abundace of psite over gene body

      for tmp_bin_id in bin_list:
         if gene_annotation_dic[chr][strand].has_key(tmp_bin_id):
            for gene in gene_annotation_dic[chr][strand][tmp_bin_id]:
              if gene_annotation_dic[chr][strand][tmp_bin_id][gene].has_key(psite_coor):
                 GTI_total_by_gene.setdefault(gene,{})
                 GTI_total_by_gene[gene].setdefault(chr,{})
                 GTI_total_by_gene[gene][chr].setdefault(strand,{})
                 GTI_total_by_gene[gene][chr][strand].setdefault(psite_coor,0)
                 GTI_total_by_gene[gene][chr][strand][psite_coor] += 1
                 GTI_total_by_gene[gene][chr][strand].setdefault(-1, 0)
                 GTI_total_by_gene[gene][chr][strand][-1] += 1    
                 total_mapped_read += 1
   #print GTI_total_by_gene     
             
   try:
      os.mkdir(sys.argv[5]) #Specify a directory for output files
   except OSError:
      pass
   
   cutoff = int(sys.argv[6]) #Specify a cutoff value to filter low abundant gene

#Read psite abundance to different mRNA disoforms   
   for gene in GTI_total_by_gene:
         for chr in GTI_total_by_gene[gene]:
            for strand in GTI_total_by_gene[gene][chr]:
               psite_density_dic = {}
               sum_GTI_by_gene = GTI_total_by_gene[gene][chr][strand][-1]

               if sum_GTI_by_gene < cutoff:
                  continue

               for psite_coor in GTI_total_by_gene[gene][chr][strand]:
                  if psite_coor == -1:
                     continue

                  #Normalization

                  if normalization_mode == "original_read":
                     GTI_total_by_gene[gene][chr][strand][psite_coor] = GTI_total_by_gene[gene][chr][strand][psite_coor]
                  elif normalization_mode == "norm_read_total_read":
                     GTI_total_by_gene[gene][chr][strand][psite_coor] = float(GTI_total_by_gene[gene][chr][strand][psite_coor]*1000000)/total_mapped_read

                  psite_density_dic.setdefault(psite_coor, GTI_total_by_gene[gene][chr][strand][psite_coor])
                  
               for bin_id in tx_annotation_dic[chr][strand]:
                  if tx_annotation_dic[chr][strand][bin_id].has_key(gene):
                    for tx_name in tx_annotation_dic[chr][strand][bin_id][gene]:
                     sorted_CDS_coor_list = sorted(tx_annotation_dic[chr][strand][bin_id][gene][tx_name])
                     try:
                        start_codon_list = start_codon_dic[gene][tx_name]
                        stop_codon_list = stop_codon_dic[gene][tx_name]
                     except KeyError:
                        if strand == "+":
                           start_codon_list = [(chr,strand,sorted_CDS_coor_list[0])]   #start codon of the transcript
                           stop_codon_list = [(chr,strand,sorted_CDS_coor_list[-1])]   #stop codon of the transcript
                        else:
                           start_codon_list = [(chr,strand,sorted_CDS_coor_list[-1])]
                           stop_codon_list = [(chr,strand,sorted_CDS_coor_list[0])]
                        start_codon_dic.setdefault(gene,{})
                        start_codon_dic[gene].setdefault(tx_name,start_codon_list)
                           
                     distance_dic = calculate_distance_psite_start_codon(start_codon_dic,gene,tx_name,chr,strand,sorted_CDS_coor_list) #relative distance to start codon

                     if distance_dic == "None":
                        continue

                     try:
                        outfile = open(sys.argv[5]+'/'+gene+'_'+tx_name+'_density_info','w')
                     except IOError:
                        continue

                     num_position = len(sorted_CDS_coor_list)

#output psite occupancy as well as gene structure information
#The column information is 1. Position in transcript; 2-4. chromosome coordinate; 5. psite occupancy; 6. Flag value for start and stop codon; 7. relative position to start codon
                     if strand == "+":
                        for i in range(1,num_position):
                           current_psite_coor = sorted_CDS_coor_list[i]
                           if psite_density_dic.has_key(current_psite_coor):
                              outfile.write(str(i)+'\t'+chr+'\t'+strand+'\t'+str(current_psite_coor)+'\t'+str(psite_density_dic[current_psite_coor])+'\t')                          
                           else:
                              outfile.write(str(i)+'\t'+chr+'\t'+strand+'\t'+str(current_psite_coor)+'\t'+'0'+'\t')
                           try:
                              if (chr,strand,current_psite_coor) in start_codon_list:
                                 outfile.write('1'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                              elif (chr,strand,current_psite_coor) in stop_codon_list:
                                 outfile.write('-1'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                              else:
                                 outfile.write('0'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                           except KeyError:
                              outfile.write('0'+'\t'+str(distance_dic[current_psite_coor])+'\n') 
                     elif strand == "-":
                        for i in range(num_position-1,0,-1):
                           current_psite_coor = sorted_CDS_coor_list[i]
                           if psite_density_dic.has_key(current_psite_coor):
                              outfile.write(str(num_position - i)+'\t'+chr+'\t'+strand+'\t'+str(current_psite_coor)+'\t'+str(psite_density_dic[current_psite_coor])+'\t')
                           else:
                              outfile.write(str(num_position - i)+'\t'+chr+'\t'+strand+'\t'+str(current_psite_coor)+'\t'+'0'+'\t')
                            
                           try:
                              if (chr,strand,current_psite_coor) in start_codon_list:
                                 outfile.write('1'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                              elif (chr,strand,current_psite_coor) in stop_codon_list:
                                 outfile.write('-1'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                              else:
                                 outfile.write('0'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                           except KeyError:
                              outfile.write('0'+'\t'+str(distance_dic[current_psite_coor])+'\n')
                     outfile.close() 

               