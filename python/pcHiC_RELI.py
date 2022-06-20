#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 22:08:01 2020

@author: hongchanyoung
"""

import pandas as pd
import numpy as np
import random
import subprocess as sub
import glob
import os
import argparse
import sys
#import scipy
import scipy.stats


#getting per chromosome TSS filename
def chr (chr):
    filename = chr+"_TSS.bed"
    return filename
    
#Function for reading file to array
def read_file_array (file):
    arr = []
    with open (file,'r') as f:
        for l in f:
            l = l.strip('\n')
            l = l.split("\t")
            arr.append(l)          
    return arr

#Function for randomizing the input with fixed looping distance and anchor width from both sides
def randomize (chr_num, loop_dis,anchor_width_left,anchor_width_right):
    rand_arr = []
    #read TSS location files
    chr_arr =read_file_array(chr(chr_num))
    #Randomly choose any TSS location within same chromosome
    rand_line = random.choice(chr_arr)
    rand_chr = rand_line[0]
    #In case gene is + strands
    if rand_line[2] == "+":
        rand_e_pos_left = int(rand_line[1])  #ending position of left anchor (which is TSS site)   
        rand_s_pos_left = rand_e_pos_left - anchor_width_left # staring position of left anchor 
        rand_s_pos_right = rand_e_pos_left + loop_dis # starting position of right anchor (TSS + looping distance)
        rand_e_pos_right = rand_s_pos_right + anchor_width_right # ending position of right anchor
    #In case gene is - strands
    else:
        rand_s_pos_right = int(rand_line[1]) #starting position of right anchor (which is TSS site)
        rand_e_pos_right = rand_s_pos_right + anchor_width_right # ending position of right anchor
        rand_e_pos_left = rand_s_pos_right - loop_dis # ending position of left anchor (which is TSS - looping distance)
        rand_s_pos_left = rand_e_pos_left - anchor_width_left # starting posiiton of left anchor
    #appending the random loop
    rand_arr.append([rand_chr,
                    rand_s_pos_left,
                    rand_e_pos_left,
                    rand_chr,
                    rand_s_pos_right,
                    rand_e_pos_right])
    return rand_arr[0]

#making random loop ***any loop anchor position should be greater than 0 and right position should be less than the maximum of chromosome maximum base (values in chr_max_dic)
def make_random_loop_array (input_file):
    rand_loop = []
    input_arr = read_file_array(input_file)
    iter_num = 0
    chr_max_dic = {"chr10":135440299,
    "chr11":134856693,
    "chr12":133812422,
    "chr13":115099423,
    "chr14":107259214,
    "chr15":102519301,
    "chr16":90244014,
    "chr17":81188573,
    "chr18":78005397,
    "chr19":59095762,
    "chr1":249213345,
    "chr20":62934707,
    "chr21":48085036,
    "chr22":51238065,
    "chr2":243102476,
    "chr3":197949384,
    "chr4":190989019,
    "chr5":180795226,
    "chr6":170893780,
    "chr7":158937649,
    "chr8":146281416,
    "chr9":141093903,
    "chrM":3230,
    "chrX":155257848,
    "chrY":59360854}
    for line in input_arr:
        new_line_num = 0
        chr_num = line[0]
        chr_max = int(chr_max_dic[chr_num])
        s_pos_left = int(line[1])
        e_pos_left = int(line[2])
        s_pos_right = int(line[4])
        e_pos_right = int(line[5])
        loop_dis = s_pos_right - e_pos_left
        anchor_width_left = e_pos_left - s_pos_left
        anchor_width_right = e_pos_right - s_pos_right
        rand_arr = randomize(chr_num, loop_dis,anchor_width_left,anchor_width_right)
        while True:
            if ((rand_arr[1]<=0)  or (rand_arr[4]<=0)) or ((rand_arr[2]>chr_max) or (rand_arr[5]>chr_max)):

                rand_arr = randomize(chr_num, loop_dis,anchor_width_left,anchor_width_right)
                if rand_arr[4] <=0:
                    print(rand_arr)
            else:
                break
        rand_loop.append(rand_arr)
    return rand_loop        

#Slice looping (6 columns) into two bed format (2 of 3 columns)
def seperate_loop_into_bed (arr,side,kind):
    bed_arr = []
    i = 1
    for array in arr:
        loop_name = kind+"_loop"+str(i)
        if side == "left":
            bed_arr.append(['\t'.join(str(v) for v in array[0:3]),loop_name])
        else: 
            bed_arr.append(['\t'.join(str(v) for v in array[3:6]),loop_name])
        i += 1
    return bed_arr

#Save looping bed format array into bedfile in order to perform bedtools intersect
def saving_loop_arr_into_bed(arr,side,kind,iter_num=""): # In future, add path as argument
    if iter_num != "":
        iter_num = str(iter_num) + "th_"
    if kind =="ref":
        iter_num == ""
        #print(iter_num)
    elif kind == "random":
        kind = kind +".tmp"
    elif kind == "input":
        kind = kind +".tmp"
        iter_num == ""
    file_name =  side +"_"+ iter_num +kind+".bed" # future, change the path
    f = open(file_name,'w')
    for item in arr:
        f.write('\t'.join(str(v) for v in item) + '\n')
    f.close()



#Perform bedtools intersect and use program "overlap", get overlap

def get_overlap(left_file,right_file,iter_num=""):
    left_ref = "left_ref.bed"
    right_ref = "right_ref.bed"
    if iter_num != "" :
        iter_num = str(iter_num)+"th_"
        left_file_name = "left_"+iter_num +"intersected.tmp.bed" # future, change the path
        right_file_name = "right_"+iter_num +"intersected.tmp.bed" # future, change the path
        
    else:
        left_file_name = "left_input_intersected.tmp.bed" # future, change the path
        right_file_name = "right_input_intersected.tmp.bed" # future, change the path
   
    left_cmd = f"bedtools intersect -a {left_file} -b {left_ref} -wo "
    right_cmd = f"bedtools intersect -a {right_file} -b {right_ref} -wo "
    left_bt = sub.run(left_cmd, shell = True, stdout = sub.PIPE, stderr = sub.PIPE)
    # Collect bedtools output
    left_res = left_bt.stdout.decode('utf-8')
    left_err = left_bt.stderr.decode('utf-8')
 
    left_f = open(left_file_name,'w')
    left_f.write(left_res)
    left_f.close()
  
    right_bt = sub.run(right_cmd, shell = True, stdout = sub.PIPE, stderr = sub.PIPE)
    
    # Collect bedtools output
    right_res = right_bt.stdout.decode('utf-8')
    right_err = right_bt.stderr.decode('utf-8')
   
    right_f = open(right_file_name,'w')
    right_f.write(right_res)
    right_f.close()
    
    #getting overlap using "overlap" program ***Prerequisite
    overlap_cmd = f"overlap -a {left_file_name} -b {right_file_name} -fa 4,8 |wc -l"
    overlap = sub.run(overlap_cmd, shell=True,stdout=sub.PIPE, stderr=sub.PIPE)
    overlap_res = overlap.stdout.decode('utf-8')
    overlap_err = overlap.stderr.decode('utf-8')
    overlap_res = overlap_res.strip("\n")
    return int(overlap_res)

#Get the mean of randome overlap after iteration
def get_mean_of_random_overlap_array_from_iterating(iter_num,input_file,actual_overlap):
    random_overlap_array = []
    for i in range(iter_num):
        print("Iteration # %s" % (str(i+1)))
        print("Calculating Random Overlap")
        new_rand_arr = make_random_loop_array(input_file)
        left_rand_arr = seperate_loop_into_bed(new_rand_arr,"left","random")
        right_rand_arr = seperate_loop_into_bed(new_rand_arr,"right","random")
        saving_loop_arr_into_bed(right_rand_arr,"right","random",i)
        saving_loop_arr_into_bed(left_rand_arr,"left","random",i)
        left_file = "left_"+ str(i) + "th_random.tmp.bed"
        right_file = "right_"+ str(i) + "th_random.tmp.bed"
        overlap = get_overlap(left_file,right_file,i)
        random_overlap_array.append(overlap)
        print("Overlap is %s" % (overlap))
    print("Random Overlap Calculation is done successfully!")
    random_overlap_array.append(actual_overlap)
    #return np.mean(random_overlap_array),np.std(random_overlap_array)
    return np.mean(random_overlap_array),np.std(random_overlap_array,ddof=1)

#Get actual input overlap
def get_actual_overlap(input_arr):
    
    left_input_arr = seperate_loop_into_bed (input_arr,"left","input")
    right_input_arr = seperate_loop_into_bed (input_arr,"right","input")
    saving_loop_arr_into_bed(right_input_arr,"right","input")
    saving_loop_arr_into_bed(left_input_arr,"left","input")
    left_input_file = "left_input.tmp.bed"
    right_input_file =  "right_input.tmp.bed"
    actual_overlap = get_overlap(left_input_file,right_input_file)
    return int(actual_overlap)

#Make reference loop to bed file format
def reference_loop_to_bed(reference_file):
    ref_loop_arr = read_file_array(reference_file)
    right_ref_arr = seperate_loop_into_bed(ref_loop_arr,"right","ref")
    saving_loop_arr_into_bed(right_ref_arr,"right","ref")
    left_ref_arr = seperate_loop_into_bed(ref_loop_arr,"left","ref")
    saving_loop_arr_into_bed(left_ref_arr,"left","ref")
    return True

#cleaning up the temporary files
def removing_temp_files():
    files = sorted(glob.glob("*.tmp.bed"))
    for file in files:
        if os.path.exists("tmp/") == False:
            tmp_dir = "mkdir tmp"
            sub.call(tmp_dir,shell=True)
        cmd = "mv %s tmp/" % file
        sub.call(cmd,shell=True)
    return True


#Getting Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file", help="input filename")
parser.add_argument("-s","--iteration_number", help="simulation iteration number")
parser.add_argument("-l", "--left_reference",  help="left reference loop bed")
parser.add_argument("-r", "--right_reference",  help="right reference loop bed")
parser.add_argument("-c", "--control_reference",  help="reference loop file")


args = parser.parse_args()


if not args.input_file:
        sys.exit("you must give input bed file")
#we need left_ref.bed and right_ref.bed for prerequisite, or it can generate based on the reference loopfile
elif os.path.exists("left_ref.bed") ==False or os.path.exists("right_ref.bed") == False:
    if not args.control_reference:
        sys.exit("you must give control reference loop")
    else:
        reference_file = args.control_reference
        reference_loop_to_bed(reference_file)
        left_ref = "left_ref.bed"
        right_ref = "right_ref.bed"
#Setting the input arguments
left_ref = "left_ref.bed"
right_ref = "right_ref.bed"
input_file = args.input_file
iter_num = int(args.iteration_number)

#Getting actual overlap
actual_overlap = get_actual_overlap(read_file_array(input_file))

#Get mean of random overlap
print("------------Getting Mean of Random overlap------------")
rand_overlap_mean,rand_overlap_stdev = get_mean_of_random_overlap_array_from_iterating(iter_num,input_file,actual_overlap)
print("Mean of random overlap: % s, Standard Deviation of random overlap: %s " % (rand_overlap_mean,rand_overlap_stdev))

#Printing actual overlap
print("------------Getting actual overlap------------")
print("Actual overlap: %s" % (actual_overlap))

#Calculate Z-Score
print("Calculating Z-Score")
if rand_overlap_stdev == 0:
    print("you need more iteration time")
    quit()
elif actual_overlap == 0:
    z_score = 0
else:
    z_score = (actual_overlap - rand_overlap_mean) / rand_overlap_stdev
fold_enrich = (actual_overlap / rand_overlap_mean)
#Convert Z-score to P-value
if z_score == 0:
    p_value = 1
else:
    p_value = scipy.stats.norm.sf(abs(z_score))


#Reporting the z-score and p-value
print("Z-score: %s, p-value: %s , fold_enrichment: %s "%(z_score,p_value,fold_enrich))
print("Removing temporary files")
removing_temp_files()



