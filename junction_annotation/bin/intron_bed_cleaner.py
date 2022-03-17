#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:55:34 2021

@author: root
"""

import pandas as pd
import time
import sys
import numpy as np

file = sys.argv[1]
output = sys.argv[2]

#import the intron bed file
#file = '/Users/lkluegel/Downloads/leafviz_all_introns.bed'
#output = '/Users/lkluegel/Downloads/leafviz_all_introns_cleaned.bed'

df = pd.read_csv(file, sep='\t', header=None)
#get list of genes
genes = list(set(list(df.iloc[:,3])))

df_list = []
#iterate over the genes
#for each gene find all tags

count = 0
start_time = time.time()
for gene in genes:
    curr_df = df[df[3] == gene]
    tags = list(set(list(curr_df[9])))
    tags = [x for x in tags if str(x) != 'nan']
    #iterate over the relevant tags
    #count how often exon 1 appears for each tag in the gene, that is how often
    #the tag is duplicated
    for tag in tags:
        curr_tag_df = curr_df[curr_df[9] == tag]
        num_exon_1 = curr_tag_df[7].value_counts()[1]
        #if there is a duplicate, find the duplicates
        if num_exon_1 > 1:
            starts = []
            ends = []
            starts.append(curr_tag_df.index[0])
            for i in range(len(curr_tag_df.index)-1):
                curr_index = curr_tag_df.index[i]
                next_index = curr_tag_df.index[i+1]

                intron_i = curr_tag_df.loc[curr_index,7]
                intron_j = curr_tag_df.loc[next_index,7]
                if intron_j <= intron_i:
                    ends.append(curr_index)
                    starts.append(next_index)
            ends.append(next_index)
            
            #get the transcript with the largest number of introns
            #if there is a tie, keep the one at the top of the file
            max_num_introns = 0
            best_i = 0

            for i in range(len(starts)):
                num_introns = ends[i] - starts[i]
                if num_introns > max_num_introns:
                    max_num_introns = num_introns
                    best_i = i

            output_curr_tag_df = curr_tag_df.loc[starts[best_i]:ends[best_i]]
            df_list.append(output_curr_tag_df)
        else:
            df_list.append(curr_tag_df)
    if count%100 == 0:
        print('Genes processed: %s' %count)
        print('Time elapsed: %ss' %round((time.time() - start_time), 0))
    count += 1
    
df2 = pd.concat(df_list)

#df2.to_csv('/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Annotator/Data/leafviz_all_introns_cleaned.bed', header=None, index=False, sep='\t')
df2.to_csv(output, header=None, index=False, sep='\t')
