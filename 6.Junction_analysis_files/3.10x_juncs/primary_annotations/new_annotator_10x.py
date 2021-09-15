#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 09:27:11 2020

@author: root
"""

import re
import os
import time
import warnings
import pandas as pd
import sys
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)

def gene_counter(input_string):
    """get the count of genes per cluster"""
    num_genes = len(input_string.split(','))
    return num_genes

def transcript_tag_finder(rel_ref_df):
    """takes the relevant reference df and returns a list of the exons with the
    highest confidence"""
    patterns = ['(.*)appris_principal_%s' %x for x in range(1,6)]
    patterns.extend(['(.*)appris_alternative', '(.*)appris_candidate'])
    tags = rel_ref_df[9]
    match_list = 0
    for pattern in patterns:
        r = re.compile(pattern)
        match_list = list(filter(r.search, tags))
        if len(match_list) > 0:
            match_list = list(set(match_list))
            break
    if isinstance(match_list, int):
        match_list = list()
    return(match_list)

def tag_sorter(tag_list):
    '''sorts a list of tags by length'''
    len_list = [len(x) for x in tag_list]
    tag_list = [x for _,x in sorted(zip(len_list, tag_list))]    
    return(tag_list)

def exon_lister(rel_ref_df, tag_list):
    '''gets the relevant exons for a tag'''
    curr_rel_ref_df = rel_ref_df[rel_ref_df[9].isin(tag_list)]
    starts = list(curr_rel_ref_df[1])
    ends = list(curr_rel_ref_df[2])
    nums = list(curr_rel_ref_df[7])
    tags = list(curr_rel_ref_df[9])
    exons = [[starts[x], ends[x]] for x in range(curr_rel_ref_df.shape[0])]
    return(starts, ends, exons, nums, tags)
        
def relevant_exon_finder(row, start_end, exons, tag_list):
    """finds the nearest reference exon for each junction"""
    start = row['start']
    end = row['end']
    if ((start_end == 'end')):
        exon_ends = [x[1] for x in exons]
        distances = [abs(end-x) for x in exon_ends]
        min_distances = min(distances)
        occurrences = distances.count(min_distances)
        if occurrences == 1:
            relevant_exon = distances.index(min(distances)) + 1
            output = relevant_exon
        else:
            occurrence_indexes = [i for i, x in enumerate(distances) if x == min_distances]
            rel_tags = [len(tag_list[x]) for x in occurrence_indexes]
            shortest_tag = min(rel_tags)
            rel_index = rel_tags.index(shortest_tag)
            output = occurrence_indexes[rel_index] + 1
    elif ((start_end == 'start')):
        exon_starts = [x[0] for x in exons]
        distances = [abs(x-start) for x in exon_starts]
        min_distances = min(distances)
        occurrences = distances.count(min_distances)
        if occurrences == 1:
            relevant_exon = distances.index(min(distances)) + 1
            output = relevant_exon
        else:
            occurrence_indexes = [i for i, x in enumerate(distances) if x == min_distances]
            rel_tags = [len(tag_list[x]) for x in occurrence_indexes]
            shortest_tag = min(rel_tags)
            rel_index = rel_tags.index(shortest_tag)
            output = occurrence_indexes[rel_index] + 1
    else:
        output = 'Error: Invaldid site_start value'
    return output
    
def site_annotator(site, ref_sites, start_end, strand):
    '''classifies site based on strand'''
    output = 'no conditions met'
    if (strand in ['+', '-']) & (start_end in ['start', 'end']):
        if site in ref_sites:
            output = 'main'
        else:
            if start_end == 'start':
                if strand == '+':
                    output = 'not_main_5_prime'
                elif strand == '-':
                    output = 'not_main_3_prime'
            elif start_end == 'end':
                if strand == '+':
                    output = 'not_main_3_prime'
                elif strand == '-':
                    output = 'not_main_5_prime'
    else:
        output = 'wrong_strand_or_site_type'
    return(output)

def distance_finder(site, rel_exon, ref_sites, strand):
    '''gets distance to relevant exon site'''
    output = 'no conditions met'
    if strand == '+':
        output  = site - ref_sites[rel_exon-1]
    elif strand == '-':
        output = site - ref_sites[rel_exon-1]
    else:
        output = 'wrong_strand_type'
    return(output)

def alt_tag_finder(tag):
    '''finds an alternative_3/5_UTR tag'''
    tag_list = list()
    tag_list.append(tag)
    alt_tag_patterns = ['(.*)alternative_3_UTR', '(.*)alternative_5_UTR']
    alt_tags = ['alternative_3_UTR', 'alternative_5_UTR']
    output = ''
    for i in range(len(alt_tags)):
        r = re.compile(alt_tag_patterns[i])
        match_list = list(filter(r.search, tag_list))
        if len(match_list) > 0:
            output = output + alt_tags[i]
    return(output)


def principal_tag_finder(tag):
    '''finds the principal part of each tag'''
    tag_list = list()
    tag_list.append(tag)
    patterns = ['(.*)appris_principal_%s' %x for x in range(1,6)]
    patterns.extend(['(.*)appris_alternative', '(.*)appris_candidate'])  
    tags = ['appris_principal_%s' %x for x in range(1,6)]
    tags.extend(['appris_alternative', 'appris_candidate'])
    output = ''
    for i in range(len(patterns)):
        r = re.compile(patterns[i])
        match_list = list(filter(r.search, tag_list))
        if len(match_list) > 0:
            output = output + tags[i]
    return(output)      
    
def splice_annotator(input_df, ref_df):
    '''combines all functions to annotate'''
    start_time = time.time()
    genes = list(set(list(input_df['gene'])))
    print('num genes = %s' %len(genes))
    output_df_list = list()
    count = 0
    for gene in genes:
        curr_input_df = input_df[input_df['gene'] == gene]
        curr_ref_df = ref_df[ref_df[3] == gene]
        tags = transcript_tag_finder(curr_ref_df)
        tags = tag_sorter(tags)
        num_tags = len(tags)
        curr_input_df['numTags'] = len(tags)
        curr_input_df['tags'] = ', '.join(tags)
        nans = [np.nan]*curr_input_df.shape[0]
        if num_tags > 0:
            rel_starts, rel_ends, rel_exons, rel_nums, tags = exon_lister(curr_ref_df, tags)
            curr_input_df['relStart'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, 'start', rel_exons, tags), axis=1)
            curr_input_df['relEnd'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, 'end', rel_exons, tags), axis=1)
            curr_input_df['relStartTag'] = curr_input_df.apply(lambda x: tags[x['relStart']-1], axis=1)
            curr_input_df['relEndTag'] = curr_input_df.apply(lambda x: tags[x['relEnd']-1], axis=1)
            curr_input_df['relStartTagClass'] = curr_input_df.apply(lambda x: principal_tag_finder(x['relStartTag']), axis=1)
            curr_input_df['relEndTagClass'] = curr_input_df.apply(lambda x: principal_tag_finder(x['relEndTag']), axis=1)
            curr_input_df['relStartTagAlt'] = curr_input_df.apply(lambda x: alt_tag_finder(x['relStartTag']), axis=1)
            curr_input_df['relEndTagAlt'] = curr_input_df.apply(lambda x: alt_tag_finder(x['relEndTag']), axis=1)
            curr_input_df['startClass'] = curr_input_df.apply(lambda x: site_annotator(x['start'], rel_starts, 'start', x['strand']), axis=1)
            curr_input_df['endClass'] = curr_input_df.apply(lambda x: site_annotator(x['end'], rel_ends, 'end', x['strand']), axis=1)
            curr_input_df['startDistance'] = curr_input_df.apply(lambda x: distance_finder(x['start'], x['relStart'], rel_starts, x['strand']), axis=1)
            curr_input_df['endDistance'] = curr_input_df.apply(lambda x: distance_finder(x['end'], x['relEnd'], rel_ends, x['strand']), axis=1)   
        elif num_tags ==0:
            curr_input_df['relStart'] = nans
            curr_input_df['relEnd'] = nans
            curr_input_df['relStartTag'] = nans
            curr_input_df['relEndTag'] = nans
            curr_input_df['relStartTagClass'] = nans
            curr_input_df['relEndTagClass'] = nans
            curr_input_df['relStartTagAlt'] = nans
            curr_input_df['relEndTagAlt'] = nans
            curr_input_df['startClass'] = nans
            curr_input_df['endClass'] = nans
            curr_input_df['startDistance'] = nans
            curr_input_df['endDistance'] = nans
        if count%100 == 0:
            print('Genes processed: %s' %count)
            print('Time elapsed: %ss' %round((time.time() - start_time), 0))
        count += 1
        output_df_list.append(curr_input_df)
    output_df = pd.concat(output_df_list)
    return(output_df, output_df_list)

path = sys.argv[1]
effectFile = sys.argv[2]
juncFile = sys.argv[3]
mainOutput = sys.argv[4]
    
#path = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Annotator/Data'
#effectFile = 'MDS_P1_intron_metadata.txt'
#juncFile = 'leafviz_all_introns.bed'
#mainOutput = 'annotations11.csv'
#exonOutput = 'exons8.csv'

#qcOutput = 'test_quality_control4.csv'

#change path
os.chdir(path)

#import the effect file
effectDf = pd.read_csv(effectFile,sep='\t')

#drop the missing values
effectDf = effectDf.dropna(axis=0,how='any')

effectDf['end'] = effectDf['end']
#import the reference dataframe
#refDf = pd.read_csv(refFile)
#
##keep only the exon information
#refDf = refDf[refDf['type'] == 'exon']

#remove the exons that have 0 length
#refDf['length'] = refDf['end'] - refDf['start']
#refDf = refDf[refDf['length'] > 0]

#import the junction annotation file
juncDf = pd.read_csv(juncFile, header=None, sep='\t')

#run the outer function to generate all data
output_df, output_list = splice_annotator(effectDf, juncDf)
output_df['end'] = output_df['end']
output_df.to_csv(mainOutput)

#qc_df.to_csv(qcOutput, index = False)
#exon_df.to_csv(exonOutput, index = False)

#test
#rel_ref_df = juncDf[juncDf[3] == 'ADAT1']
#curr_input_df = effectDf[effectDf['gene'] == 'ADAT1']
#tags = transcript_tag_finder(rel_ref_df)
#tags = tag_sorter(tags)
#rel_tag = tags[0]
#rel_starts, rel_ends, rel_exons, rel_nums= exon_lister(rel_ref_df, rel_tag)
#print(curr_input_df.apply(lambda x: relevant_exon_finder(x, 'start', rel_exons), axis=1))