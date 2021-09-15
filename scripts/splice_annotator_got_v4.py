"""text"""
import fnmatch
import re
import os
import copy
import bisect
import time
import warnings
import pandas as pd
import sys
import numpy as np
from collections import Counter

warnings.simplefilter(action='ignore', category=FutureWarning)

def gene_counter(input_string):
    """get the count of genes per cluster"""
    num_genes = len(input_string.split(','))
    return num_genes

def canonical_splice_site_lister(ref_dataframe, gene):
    """takes the reference GTF df and gets a list of canonical splice sites for a
    given gene
    """
    canonical_starts = list(set(list(ref_dataframe['start'][ref_dataframe['gene'].isin(gene)])))
    canonical_ends = list(set(ref_dataframe['end'][ref_dataframe['gene'].isin(gene)]))
    canonical_sites = canonical_starts + canonical_ends
    canonical_sites = list(set(canonical_sites))
    canonical_sites.sort()
    return(canonical_starts, canonical_ends, canonical_sites)

def canonical_splice_finder(ref_dataframe, gene):
    """takes the reference GTF df and gets a list of tuples of canonical splices for
    a given gene
    """
    starts = list(ref_dataframe['start'][ref_dataframe['gene'].isin(gene)])
    ends = list(ref_dataframe['end'][ref_dataframe['gene'].isin(gene)])
    canonical_splices = [(starts[i], ends[i]) for i in range(len(starts))]
    return canonical_splices

def minimal_exon_finder(ref_dataframe, gene, canonical_sites):
    """takes the reference GTF df and finds the smallest version of the exons for a
    given gene
    """
    starts = list(ref_dataframe['start'][ref_dataframe['gene'].isin(gene)])
    ends = list(ref_dataframe['end'][ref_dataframe['gene'].isin(gene)])
    exons = list()
    for i in range(len(canonical_sites) - 1):
        left = canonical_sites[i]
        right = canonical_sites[i+1]
        if (left in starts) & (right in ends):
            exons.append((left, right))
    return exons

def splice_region_lister(input_list):
    """#turns an ordered list of splice sites into a list of regions
    between each splice site
    """
    splice_regions = [(input_list[i], input_list[i+1]) for i in range(len(input_list)-1)]
    return splice_regions


def region_code_string_generator(canonical_sites, canonical_splices_list, exons, splice_regions):
    """converts the list of canonical splice sites, canonical splices and exons into
    a code string that represents the region
    1 means region that is part of an alternative donor site
    3 means region that is part of an alternative acceptor site
    2 is the minimal canonical exon
    0 is the minimal canonical intron
    """
    #set up a list of zeros, where zero denotes a minimal intron
    #this will later be joined to form a string
    code_string = ['0']*len(splice_regions)
    #label the exons as 2
    for j, item in enumerate(exons):
        exon_left_index = canonical_sites.index(item[0])
        code_string[exon_left_index] = '2'
    #find the alternative donor and acceptor sites
    for i, item in enumerate(canonical_splices_list):
        splice_index_left = canonical_sites.index(item[0])
        splice_index_right = canonical_sites.index(item[1])
        region = code_string[splice_index_left:splice_index_right]
        num_exons = 0
        #the GTF contains splices that have skipped exons. The following
        #separates skipped exon cases and normal cases
        for k in region:
            if k == '2':
                num_exons += 1
        if num_exons > 1:
            left_exon = region.index('2')
            right_exon = len(region)-1-region[::-1].index('2')
            if region[0] != '2':
                code_string[splice_index_left:splice_index_left+left_exon] = \
                ['1'] * left_exon
            if region[-1] != '2':
                code_string[splice_index_left+right_exon+1:splice_index_right] = \
                ['3'] * (splice_index_right-splice_index_left-right_exon-1)
        elif len(region) == 1:
            pass
        else:
            exon = region.index('2')
            if region[0] != '2':
                code_string[splice_index_left:splice_index_left+exon] = ['1'] * exon
            if region[-1] != '2':
                code_string[splice_index_left+exon+1:splice_index_right] = \
                ['3'] * (splice_index_right-splice_index_left-exon-1)
    code_string = ''.join(code_string)
    return code_string

def start_site_annotator_GOT_old(start, strand, canonical_starts, canonical_ends):
    '''classifies the end as cryptic, annotated'''
    if strand == '+':
        if start not in canonical_ends:
            output = 'cryptic'
        else:
            output = 'annotated'
    elif strand == '-':
        if start not in canonical_ends:
            output = 'cryptic'
        else:
            output = 'annotated'
    else:
        output = 'Error: wrong strand indicator'
    return output

def start_site_annotator_GOT(start, canonical_ends):
    '''classifies the end as cryptic, annotated'''
    if start not in canonical_ends:
        output = 'cryptic'
    else:
        output = 'annotated'
    return output

def relevant_exon_finder_old(start, end, site_start, strand, exons):
    if ((site_start == True) & (strand == '+')) or ((site_start == False) & (strand == '-')):
        exon_ends = [x[1] for x in exons]
        distances = [start-x for x in exon_ends]
        relevant_distances = [x for x in distances if x>=0]
        if len(relevant_distances) > 0:
            relevant_exon = distances.index(min(relevant_distances)) + 1
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'
    elif ((site_start == True) & (strand == '-')) or ((site_start == False) & (strand == '+')):
        exon_starts = [x[0] for x in exons]
        distances = [x-end for x in exon_starts]
        relevant_distances = [x for x in distances if x>=0]
        if len(relevant_distances) > 0: 
            relevant_exon = distances.index(min(relevant_distances)) + 1       
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'          
    else:
        return 'Error: Bad strand or site_start'
    
def relevant_exon_finder_old2(start, end, site_start, strand, exons):
    if ((site_start == True) & (strand == '+')):
        exon_ends = [x[1] for x in exons]
        distances = [start-x for x in exon_ends]
        relevant_distances = [x for x in distances if x>=0]
        if len(relevant_distances) > 0:
            relevant_exon = distances.index(min(relevant_distances)) + 1
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'
    elif ((site_start == True) & (strand == '-')):
        exon_starts = [x[0] for x in exons]
        distances = [x-end for x in exon_starts]
        relevant_distances = [x for x in distances if x>=0]
        if len(relevant_distances) > 0: 
            relevant_exon = distances.index(min(relevant_distances)) + 1       
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'  
    elif  ((site_start == False) & (strand == '-')):
        exon_ends = [x[1] for x in exons]
        distances = [start-x for x in exon_ends]
        relevant_distances = [x for x in distances if x<=0]
        if len(relevant_distances) > 0:
            relevant_exon = distances.index(max(relevant_distances)) + 1
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'
    elif ((site_start == False) & (strand == '+')):
        exon_starts = [x[0] for x in exons]
        distances = [x-end for x in exon_starts]
        relevant_distances = [x for x in distances if x <= 0]
        if len(relevant_distances) > 0: 
            relevant_exon = distances.index(max(relevant_distances)) + 1       
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'          
    else:
        return 'Error: Bad strand or site_start'
  
def relevant_exon_finder(row, site_start, exons, canonical_sites):
    start = row['start']
    end = row['end']
    if ((site_start == True)):
        exon_ends = [x[1] for x in exons]
        distances = [abs(start-x) for x in exon_ends]
        relevant_exon = distances.index(min(distances)) + 1
        output = relevant_exon
        if start < min(canonical_sites):
            output =  'cryptic site outside of gene locus'
    elif ((site_start == False)):
        exon_starts = [x[0] for x in exons]
        distances = [abs(x-end) for x in exon_starts]
        relevant_exon = distances.index(min(distances)) + 1       
        output = relevant_exon
        if end > max(canonical_sites):
            output =  'cryptic site outside of gene locus'   
    else:
        output = 'Error: Invaldid site_start value'
    return output

def relevant_exon_test(row):
    if row['relStartExon'] == row['relEndExon']:
        return 'Fail'
    else:
        return 'Pass'

def exon_test(row, site_start, exons, canonical_splices):
    if site_start == True:
        item = row['relStartExon']
    else:
        item = row['relEndExon']
    if item == 'cryptic site outside of gene locus':
        return 'Not applicable'
    elif type(item) == int:
        rel_exon = exons[item-1]
        if rel_exon in canonical_splices:
            return 'Pass'
        else:
            return 'Fail'
        
'''
elif ((site_start == False)):
        exon_starts = [x[0] for x in exons]
        distances = [x-end for x in exon_starts]
        relevant_distances = [x for x in distances if x>=0]
        if len(relevant_distances) > 0: 
            relevant_exon = distances.index(min(relevant_distances)) + 1       
            return relevant_exon
        else:
            return 'cryptic site outside of gene locus'          
    else:
        return 'Error: Invaldid site_start value'
'''    
def end_site_annotator_GOT(end, canonical_starts):
    '''classifies the end as cryptic, annotated'''
    if end not in canonical_starts:
        classification = 'cryptic'
    else:
        classification = 'annotated'
    return classification

def alt_annotated_finder(curr_df, site_start):
    if site_start == True:
        exon = 'relStartExon'
        siteClass = 'startClass'
        index = 'start'
    elif site_start == False:
        exon = 'relEndExon'
        siteClass = 'endClass'
        index = 'end'
    curr_df = curr_df.sort_values(by=[exon,index])
    annotated_exons_counts = curr_df[exon][curr_df[siteClass] == 'annotated'].value_counts()
    pot_relevant_exons = list(annotated_exons_counts[annotated_exons_counts > 1].keys())
    for key in pot_relevant_exons:
        pot_alt_annotateds = list(set(list(curr_df[index][(curr_df[exon] == key) & (curr_df[siteClass] == 'annotated')])))
        pot_alt_annotateds.sort()
        if len(pot_alt_annotateds) > 1:
            count = 1
            for i in pot_alt_annotateds:
                curr_df[siteClass][(curr_df[exon] == key) & (curr_df[index] == i)] = 'annotated_alt_exon%s_num%s' %(key, count)
                count += 1
    return curr_df    

def most_common_alt_annotated_finder(curr_df, site_start):
    if site_start == True:
        exon = 'relStartExon'
        siteClass = 'startClass'
        index = 'start'
        output_col = 'distanceToMostCommonAltAnnotatedStart'
        qcCol = 'mostCommonAltAnnotationStartQualityTest'
    elif site_start == False:
        exon = 'relEndExon'
        siteClass = 'endClass'
        index = 'end'
        output_col = 'distanceToMostCommonAltAnnotatedEnd'
        qcCol = 'mostCommonAltAnnotationEndQualityTest'
    curr_df[output_col] = [np.nan] * curr_df.shape[0]
    exons = list(set(list(curr_df[exon]))) 
    r = re.compile('annotated_alt_exon*.')
    for e in exons:
        alt_annotateds = list(filter(r.match, list(curr_df[siteClass][curr_df[exon] == e])))
        if len(alt_annotateds) > 0:
            c = Counter(alt_annotateds)
            most_common_alt_annotated = c.most_common(1)[0][0]
            if c.most_common(2)[0][1] == c.most_common(2)[1][1]:
                curr_df[qcCol][curr_df[exon] == e] = 'Fail'
            ref_index = curr_df[curr_df[siteClass] == most_common_alt_annotated].iloc[0,list(curr_df.columns.values).index(index)]
            indexes = list()
            for row in curr_df.iterrows():
                if row[1][siteClass] in alt_annotateds:
                    indexes.append(row[0])
            for i in indexes:
                row = curr_df.loc[i]
                if row['strand'] == '+':
                    output = row[index] - ref_index
                elif row['strand'] == '-' :
                    output = ref_index - row[index]
                output = int(output)
                curr_df.loc[i,output_col] = output
    return curr_df

            
def alt_annotated_finder_right(curr_df):
    curr_df = curr_df.sort_values(by=['rightExon','end'])
    annotated_exons_counts = curr_df['rightExon'][curr_df['endClass'] == 'annotated'].value_counts()
    pot_relevant_exons = list(annotated_exons_counts[annotated_exons_counts > 1].keys())
    for key in pot_relevant_exons:
        pot_alt_annotateds = list(set(list(curr_df['end'][(curr_df['rightExon'] == key) & (curr_df['endClass'] == 'annotated')])))
        if len(pot_alt_annotateds) > 1:
            count = 1
            for i in pot_alt_annotateds:
                curr_df['endClass'][(curr_df['rightExon'] == key) & (curr_df['end'] == i)] = 'annotated_alt_exon%s_num%s' %(key, count)
                count += 1
    return curr_df

def verdict_judger(row, site_start, exon_starts, exon_ends):
    output = 'Error: no condition triggered'
    strand = row['strand']
    if site_start == True:
        site = row['startClass']
        if strand == '+':
            if site == 'cryptic':
                output = 'cryptic_fiveprime'
            else:
                if row['start'] in exon_ends:
                    output = 'minimal_exon_fiveprime'
                elif site == 'annotated':
                    output = 'annotated'
                else:
                    output = 'alternative_fiveprime'
        elif strand == '-':
            if site == 'cryptic':
                output = 'cryptic_threeprime'
            else:
                if row['start'] in exon_ends:
                    output = 'minimal_exon_threeprime'
                elif site == 'annotated':
                    output = 'annotated'
                else:
                    output = 'alternative_threeprime'
    elif site_start == False:
        site = row['endClass']
        if strand == '+':
            if site == 'cryptic':
                output = 'cryptic_threeprime'
            else:
                if row['end'] in exon_starts:
                    output = 'minimal_exon_threeprime'
                elif site == 'annotated':
                    output = 'annotated'
                else:
                    output = 'alternative_threeprime'
        elif strand == '-':
            if site == 'cryptic':
                output = 'cryptic_fiveprime'
            else:
                if row['end'] in exon_starts:
                    output = 'minimal_exon_fiveprime'
                elif site == 'annotated':
                    output = 'annotated'
                else:
                    output = 'alternative_fiveprime'
    return output

def minimal_exon_distance_finder_old(row, site_start, exons):
    output = 'Error: no condition triggered'
    if site_start == True:
        verdict = row['startVerdict']
        relevant_exon = row['relStartExon']
        if (verdict == 'minimal_exon_fiveprime') or (verdict == 'minimal_exon_threeprime'):
            output = 0 
        elif relevant_exon == 'cryptic site outside of gene locus':
            output = np.nan
        else:
            strand = row['strand']
            if strand == '+':
                output = row['start'] - exons[relevant_exon-1][1]
            elif strand == '-':
                output = exons[relevant_exon-1][1] - row['start']       
    else:
        verdict = row['endVerdict']
        relevant_exon = row['relEndExon']
        if (verdict == 'minimal_exon_fiveprime') or (verdict == 'minimal_exon_threeprime'):
            output = 0 
        elif relevant_exon == 'cryptic site outside of gene locus':
            output = np.nan
        else:
            strand = row['strand']
            if strand == '+':
                output = row['end'] - exons[relevant_exon-1][1]
            elif strand == '-':
                output = exons[relevant_exon-1][0] - row['end']    
    return output
 
def minimal_exon_distance_finder(row, site_start, exons):
    output = 'Error: no condition triggered'
    if site_start == True:
        verdict = row['startVerdict']
        relevant_exon = row['relStartExon']
        if (verdict == 'minimal_exon_fiveprime') or (verdict == 'minimal_exon_threeprime'):
            output = 0 
        elif relevant_exon == 'cryptic site outside of gene locus':
            output = np.nan
        else:
            strand = row['strand']
            if strand == '+':
                output = row['start'] - exons[relevant_exon-1][1]
            elif strand == '-':
                output = exons[relevant_exon-1][1] - row['start']       
    else:
        verdict = row['endVerdict']
        relevant_exon = row['relEndExon']
        if (verdict == 'minimal_exon_fiveprime') or (verdict == 'minimal_exon_threeprime'):
            output = 0 
        elif relevant_exon == 'cryptic site outside of gene locus':
            output = np.nan
        else:
            strand = row['strand']
            if strand == '+':
                output = row['end'] - exons[relevant_exon-1][0]
            elif strand == '-':
                output = exons[relevant_exon-1][0] - row['end']    
    return output
                       
def alt_annotated_finder_left(curr_df):
    curr_df = curr_df.sort_values(by=['leftExon','start'])
    annotated_exons_counts = curr_df['leftExon'][curr_df['startClass'] == 'annotated'].value_counts()
    pot_relevant_exons = list(annotated_exons_counts[annotated_exons_counts > 1].keys())
    for key in pot_relevant_exons:
        pot_alt_annotateds = list(set(list(curr_df['start'][(curr_df['leftExon'] == key) & (curr_df['startClass'] == 'annotated')])))
        if len(pot_alt_annotateds) > 1:
            count = 1
            for i in pot_alt_annotateds:
                curr_df['startClass'][(curr_df['leftExon'] == key) & (curr_df['start'] == i)] = 'annotated_alt_exon%s_num%s' %(key, count)
                count += 1
    return curr_df
    
def start_site_annotator(start, code_string, canonical_starts, canonical_ends,
                         canonical_sites):
    """determines the type of LeafCutter splice start"""
    num_matches = 0
    #this chunk deals with splices outside of the gene locus
    if start < min(canonical_sites):
        output = 'cryp. acceptor before gene locus start'
        num_matches += 1
    elif start > max(canonical_sites):
        output = 'cryp. acceptor after gene locus end'
        num_matches += 1
    #this chunk deals with splices at the boundary of the gene locus
    elif start == min(canonical_sites):
        first_region = code_string[0]
        if first_region == '1':
            output = 'acceptor at can. alt. donor'
            num_matches += 1
        elif first_region == '2':
            output = 'acceptor at min. exon start'
            num_matches += 1
        else:
            output = 'ERROR 1'
            num_matches += 1
    elif start == max(canonical_sites):
        last_region = code_string[-1]
        if last_region == '3':
            output = 'can. alt. acceptor'
            num_matches += 1
        elif last_region == '2':
            output = 'min. exon end'
            num_matches += 1
        else:
            output = 'ERROR 2'
            num_matches += 1
    #this chunk deals with sites that coincide with canonical sites
    else:
        left_site = list(filter(lambda i: i <= start, canonical_sites))[-1]
        left_index = canonical_sites.index(left_site)
        left_region = code_string[left_index-1]
        right_region = code_string[left_index]
        if start in canonical_ends:
            if start in canonical_starts:
                if (left_region == '2') & (right_region == '2'):
                    output = 'can. acceptor in min. exon'
                    num_matches += 1
                elif (left_region == '2') & (right_region != '2'):
                    output = 'min. exon end'
                    num_matches += 1
                elif (left_region != '2') & (right_region == '2'):
                    output = 'can. acceptor at min exon start'
                    num_matches += 1
                elif (left_region in ['0', '1']) & (right_region != '2'):
                    output = 'can. acceptor at can. alt. donor'
                    num_matches += 1
                elif left_region == '3':
                    output = 'can. alt. acceptor'
                    num_matches += 1
            else:
                if left_region == '2':
                    output = 'min. exon end'
                    num_matches += 1
                elif left_region == '3':
                    output = 'can. alt. acceptor'
                    num_matches += 1
                else:
                    output = 'ERROR 3'
                    num_matches += 1
        elif start in canonical_starts:
            if right_region == '2':
                output = 'cryp. acceptor at min exon start'
                num_matches += 1
            elif right_region == '1':
                output = 'cryp. acceptor at can. alt. donor'
                num_matches += 1
            else:
                output = 'Error 4'
                num_matches += 1
        #this chunk deals with purely cryptic sites
        else:
            curr_site = list(filter(lambda i: i <= start, canonical_sites))[-1]
            curr_site = canonical_sites.index(curr_site)
            curr_site = code_string[left_index]
            if curr_site == '2':
                output = 'cryp. acceptor in min. exon'
                num_matches += 1
            elif curr_site == '3':
                output = 'cryp. alt. acceptor'
                num_matches += 1
            elif curr_site == '1':
                output = 'funky cryp. alt acceptor'
                num_matches += 1
            elif curr_site == '0':
                output = 'cryp. alt. acceptor in can. intron'
                num_matches += 1
    if num_matches == 0:
        output = 'no matches'
    elif num_matches > 2:
        output = 'multiple conditions met'
    return output

def end_site_annotator(end, code_string, canonical_starts, canonical_ends,
                       canonical_sites):
    """determines the type of LeafCutter splice end"""
    num_matches = 0
    #this chunk deals with splices outside of the gene locus
    if end < min(canonical_sites):
        output = 'cryp. donor before gene locus start'
        num_matches += 1
    elif end > max(canonical_sites):
        output = 'cryp. donor after gene locus end'
        num_matches += 1
    #this chunk deals with splices at the boundary of the gene locus
    elif end == min(canonical_sites):
        first_region = code_string[0]
        if first_region == '1':
            output = 'donor at can. alt. donor'
            num_matches += 1
        elif first_region == '2':
            output = 'donor at min. exon start'
            num_matches += 1
        else:
            output = 'ERROR 1'
            num_matches += 1
    elif end == max(canonical_sites):
        last_region = code_string[-1]
        if last_region == '3':
            output = 'donor at can. alt. acceptor'
            num_matches += 1
        elif last_region == '2':
            output = 'donor at min. exon end'
            num_matches += 1
        else:
            output = 'ERROR 2'
            num_matches += 1
    #this chunk deals with sites that coincide with canonical sites
    else:
        left_site = list(filter(lambda i: i <= end, canonical_sites))[-1]
        left_index = canonical_sites.index(left_site)
        left_region = code_string[left_index-1]
        right_region = code_string[left_index]
        if end in canonical_ends:
            if end in canonical_starts:
                if (left_region == '2') & (right_region == '2'):
                    output = 'can. donor in min. exon'
                    num_matches += 1
                elif (left_region == '2') & (right_region != '2'):
                    output = 'can. donor at min. exon end'
                    num_matches += 1
                elif (right_region == '2') & (left_region != '2'):
                    output = 'min. exon start'
                    num_matches += 1
                elif (left_region in ['0', '1']) & (right_region != '2'):
                    output = 'can. alt. donor'
                    num_matches += 1
                elif left_region == '3':
                    output = 'can. donor at can. alt. acceptor'
                    num_matches += 1
            else:
                if left_region == '2':
                    output = 'cryp. donor at min. exon end'
                    num_matches += 1
                elif left_region == '3':
                    output = 'cryp. donor at can. alt. acceptor'
                    num_matches += 1
                else:
                    output = 'ERROR 3'
                    num_matches += 1
        elif end in canonical_starts:
            if right_region == '2':
                output = 'min. exon start'
                num_matches += 1
            elif right_region == '1':
                output = 'can. alt. donor'
                num_matches += 1
            else:
                output = 'Error 4'
                num_matches += 1
        #this chunk deals with purely cryptic sites
        else:
            curr_site = list(filter(lambda i: i <= end, canonical_sites))[-1]
            curr_site = canonical_sites.index(curr_site)
            curr_site = code_string[left_index]
            if curr_site == '2':
                output = 'cryp. donor in min. exon'
                num_matches += 1
            elif curr_site == '3':
                output = 'funky cryp. alt. donor'
                num_matches += 1
            elif curr_site == '1':
                output = 'cryp. alt acceptor'
                num_matches += 1
            elif curr_site == '0':
                output = 'cryp. alt. donor in can. intron'
                num_matches += 1
    if num_matches == 0:
        output = 'no matches'
    elif num_matches > 2:
        output = 'multiple conditions met'
    return output

def quality_control_string(code_string, df, canonical_sites):
    min_start = min(df['start'])
    max_end = max(df['end'])
    if min_start <= min(canonical_sites):
        start = 0
    else:
        temp_sites = copy.copy(canonical_sites)
        temp_sites.append(min_start)
        temp_sites.sort()
        start = temp_sites.index(min_start) - 1
    if max_end >= max(canonical_sites):
        end = len(canonical_sites) - 1
    else:
        temp_sites = copy.copy(canonical_sites)
        temp_sites.append(max_end)
        temp_sites.sort()
        end = temp_sites.index(max_end) + 1        
    relevant_code_string = code_string[start:end+1]
    """checks the region code string for impossible strings"""
    output = 'Fail'
    if fnmatch.fnmatch(relevant_code_string, '*00*'):
        output = output + '1'
    if fnmatch.fnmatch(relevant_code_string, '*03*'):
        output = output + '2'
    if fnmatch.fnmatch(relevant_code_string, '*13*'):
        output = output + '3'
    if fnmatch.fnmatch(relevant_code_string, '*22*'):
        output = output + '4'
    if fnmatch.fnmatch(relevant_code_string, '*31*'):
        output = output + '5'
    if len(output) == 4:
        output = 'Pass'
    return output

def inserter(input_list, i):
    """#insert an item in a sorted list in the correct location"""
    bisect.insort(input_list, i)
    return input_list

def cryptic_site_finder(input_effect_df, canonical_sites):
    """finds the cryptic sites"""
    cryptic_starts = set(list(input_effect_df['start']))
    cryptic_ends = set(list(input_effect_df['end']))

    #find the cryptic sites
    cryptic_starts = list(cryptic_starts.difference(canonical_sites))
    cryptic_ends = list(cryptic_ends.difference(canonical_sites))
    return(cryptic_starts, cryptic_ends)

def quality_control_cryp_sites(canonical_sites, cryptic_starts, cryptic_ends):
    """if all LeafCutter sites are cryptic, fail QC, probably an alignment
    issue
    """
    cryptic_sites = copy.copy(cryptic_starts)
    cryptic_sites = cryptic_sites.extend(cryptic_ends)
    try:
        if set(canonical_sites).difference(set(cryptic_sites)) == set(canonical_sites):
            output = 'Fail'
        else:
            output = 'Pass'
    except:
        output = 'Pass'
    return output

def quality_control_start_end_sites(ref_df, gene, input_df):
    """check that no start sites are end sites"""
    sample_starts = list(input_df['start'])
    sample_ends = list(input_df['end'])
    canonical_starts = list(ref_df['start'][(ref_df['gene'].isin(gene)) & (ref_df['start'] <= max(sample_ends))])
    canonical_ends = list(ref_df['end'][(ref_df['gene'].isin(gene) & (ref_df['end'] >= min(sample_starts)))])
    starts = canonical_starts + sample_ends
    ends = canonical_ends + sample_starts
    if len(set(starts).intersection(set(ends))) == 0:
        return 'Pass'
    else:
        return 'Fail'

def exon_skip_counter(code_string, start, end, canonical_starts, canonical_ends, canonical_sites):
    """count the number of skipped exons"""
    if start < min(canonical_starts):
        start = min(canonical_starts)
    if end > max(canonical_ends):
        end = max(canonical_ends)

    left_site = list(filter(lambda i: i <= start, canonical_sites))[-1]
    right_site = list(filter(lambda i: i >= end, canonical_sites))[0]
    left_index = canonical_sites.index(left_site)
    right_index = canonical_sites.index(right_site)
    temp_string = code_string[left_index:right_index]
    pattern = '*2*'
    if len(temp_string) == 0:
        return 0
    elif fnmatch.fnmatch(temp_string, pattern):
        count = temp_string.count('2')
        repeats = [x for x in re.findall('2+', temp_string) if len(x) > 1]
        count2 = ''.join(repeats)
        count = count - len(count2) + len(repeats)
        return count
    else:
        return 0

def quality_control_slice_len(series):
    """checks if the minimum slice is nonzero"""
    if min([len(x) for x in series]) < 1:
        return 'Fail'
    else:
        return 'Pass'

def quality_control_cryp_site_location(cryptic_starts, cryptic_ends, canonical_sites):
    """checks if the cryptic sites fall within the known gene locus"""
    status = 'Pass'
    if len(cryptic_starts) > 0:
        if (min(cryptic_starts) < min(canonical_sites)) or (max(cryptic_starts) > max(canonical_sites)):
            status = 'Fail'
    elif len(cryptic_ends) > 0:
        if (max(cryptic_ends) > max(canonical_sites)) or (min(cryptic_ends) < min(canonical_sites)):
            status = 'Fail'
    return status

def convert_exon_list(curr_genes_list, curr_starts_list, curr_ends_list,
                      curr_names_list, gene, exons_list):
    """takes a list of genes and a list of corresponding minimal exons and
    converts it into a dataframe"""
    starts = [x[0] for x in exons_list]
    ends = [x[1] for x in exons_list]
    genes = [gene]*len(exons_list)
    names = [x for x in range(1, len(genes)+1)]
    curr_genes_list.extend(genes)
    curr_starts_list.extend(starts)
    curr_ends_list.extend(ends)
    curr_names_list.extend(names)
    return(curr_genes_list, curr_starts_list, curr_ends_list, curr_names_list)

def quality_control_overlapping_min_exon(exons, canonical_splices, min_starts, max_ends):
    """checks whether there are potential issues with the minimal exons"""
    output = 'Pass'
    for i in exons:
        if i not in canonical_splices:
            if (i[0] >= min_starts) & (i[1] <= max_ends):
                output = 'Fail'
                return output
    return output

def splice_annotator(input_df, ref_df):
    """combines all functions to produce the splice annotation"""
    start_time = time.time()
    genes = list(set(list(input_df['gene'])))
    print('Number of genes: %s' %len(genes))
    output_df_list = list()
    quality_control_string_list = list()
    quality_control_start_end_sites_list = list()
    quality_control_cryp_sites_list = list()
    quality_control_cryp_site_location_list = list()
#    quality_control_min_exon_list = list()
    interim_genes_list = list()
    interim_starts_list = list()
    interim_ends_list = list()
    interim_exon_name_list = list()
    count = 0
    for gene in genes:
        gene_list = gene.split(',')
        curr_input_df = input_df[input_df['gene'] == gene]
        can_starts, can_ends, can_sites = canonical_splice_site_lister(ref_df,  gene_list)
        can_splices = canonical_splice_finder(ref_df,  gene_list)
        exons = minimal_exon_finder(ref_df,  gene_list, can_sites)
        interim_genes_list, interim_starts_list, interim_ends_list, interim_exon_name_list = \
        convert_exon_list(interim_genes_list, interim_starts_list,
                          interim_ends_list, interim_exon_name_list,
                           gene_list, exons)
        regions = splice_region_lister(can_sites)
        code_string = region_code_string_generator(can_sites,
                                                   can_splices,
                                                   exons,
                                                   regions)
        cryp_starts, cryp_ends = cryptic_site_finder(curr_input_df, can_sites)
        curr_input_df['startClass'] = curr_input_df.apply(lambda x: start_site_annotator_GOT(x['start'], can_ends), axis=1)
        curr_input_df['endClass'] = curr_input_df.apply(lambda x: end_site_annotator_GOT(x['end'], can_starts), axis=1)
        curr_input_df['relStartExon'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, True, exons, can_sites), axis=1)
        curr_input_df['relEndExon'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, False, exons, can_sites), axis=1)
        curr_input_df = alt_annotated_finder(curr_input_df, True)
        curr_input_df = alt_annotated_finder(curr_input_df, False)
        exon_starts = [x[0] for x in exons]
        exon_ends = [x[1] for x in exons]
        curr_input_df['startVerdict'] = curr_input_df.apply(lambda x: verdict_judger(x, True, exon_starts, exon_ends), axis=1)
        curr_input_df['endVerdict'] = curr_input_df.apply(lambda x: verdict_judger(x, False, exon_starts, exon_ends), axis=1)
        curr_input_df['start_min_exon_distance'] = curr_input_df.apply(lambda x: minimal_exon_distance_finder(x, True, exons), axis=1)
        curr_input_df['end_min_exon_distance'] = curr_input_df.apply(lambda x: minimal_exon_distance_finder(x, False, exons), axis=1)
        curr_input_df['mostCommonAltAnnotationStartQualityTest'] = ['Pass'] * curr_input_df.shape[0]
        curr_input_df['mostCommonAltAnnotationEndQualityTest'] = ['Pass'] * curr_input_df.shape[0]
        curr_input_df = most_common_alt_annotated_finder(curr_input_df, True)
        curr_input_df = most_common_alt_annotated_finder(curr_input_df, False)
        curr_input_df['numSkippedExons'] = curr_input_df.apply(lambda x: \
            exon_skip_counter(code_string, x['start'], x['end'], can_starts,
                              can_ends, can_sites), axis=1)
        curr_input_df['relExonTest'] = curr_input_df.apply(lambda x: relevant_exon_test(x), axis=1)
        curr_input_df['relStartExonQualityTest'] = curr_input_df.apply(lambda x: exon_test(x, True, exons, can_splices), axis=1)
        curr_input_df['relEndExonQualityTest'] = curr_input_df.apply(lambda x: exon_test(x, False, exons, can_splices), axis=1)
        output_df_list.append(curr_input_df)
        quality_control_string_list.append(quality_control_string(code_string, curr_input_df, can_sites))
        quality_control_start_end_sites_list.append(\
            quality_control_start_end_sites(ref_df,  gene_list, curr_input_df))
        quality_control_cryp_sites_list.append(quality_control_cryp_sites(\
            can_sites, cryp_starts, cryp_ends))
        quality_control_cryp_site_location_list.append(\
            quality_control_cryp_site_location(cryp_starts, cryp_ends, can_sites))
        count += 1
        if count%100 == 0:
            print('Genes processed: %s' %count)
            print('Time elapsed: %ss' %round((time.time() - start_time), 0))
    output_df = pd.concat(output_df_list)
    qc_df = pd.DataFrame(list(zip(genes, quality_control_string_list,
                                  quality_control_start_end_sites_list,
                                  quality_control_cryp_sites_list,
                                  quality_control_cryp_site_location_list)),
                         columns = ['gene',
                                    'ref GFT test',
                                    'start/end site overlap test',
                                    'some annotated sites test',
                                    'sites outside of locus test'])
    exon_df = pd.DataFrame(list(zip(interim_genes_list,
                                    interim_starts_list,
                                    interim_ends_list,
                                    interim_exon_name_list)),
                           columns = ['gene', 'start', 'end', 'name'])
    exon_df['gene'] = exon_df.apply(lambda x: x['gene'][0], axis=1)
    print("--- Execution time: %s seconds ---" %round((time.time() - start_time), 0))
    return(output_df, qc_df, exon_df)

path = sys.argv[1]
effectFile = sys.argv[2]
refFile = sys.argv[3]
mainOutput = sys.argv[4]
exonOutput = sys.argv[5]
#qcOutput = sys.argv[6]

#path = '/Users/lkluegel/Documents/Splicing/Fede_Ally_Paulina/Data'
#effectFile = 'MDS_P1_intron_metadata.txt'
#refFile = 'gencode.v31.basic.annotation.csv'
#mainOutput = 'annotations4.csv'
#exonOutput = 'exons4.csv'

#qcOutput = 'test_quality_control4.csv'

#change path
os.chdir(path)

#import the effect file
effectDf = pd.read_csv(effectFile,sep='\t')

#drop the missing values
effectDf = effectDf.dropna(axis=0,how='any')

effectDf['end'] = effectDf['end'] + 1
#import the reference dataframe
refDf = pd.read_csv(refFile)

#keep only the exon information
refDf = refDf[refDf['type'] == 'exon']

#remove the exons that have 0 length
refDf['length'] = refDf['end'] - refDf['start']
refDf = refDf[refDf['length'] > 0]

#run the outer function to generate all data
output_df, qc_df, exon_df = splice_annotator(effectDf, refDf)
output_df['end'] = output_df['end'] - 1
output_df.to_csv(mainOutput)
#qc_df.to_csv(qcOutput, index = False)
exon_df.to_csv(exonOutput, index = False)

'''
test_df = output_df[output_df['strand'] == '-']
starts_minus1 = 0
starts_neutral = 0
starts_plus1 = 0

ends_minus1 = 0
ends_neutral = 0
ends_plus1 = 0
for index, row in test_df.iterrows():
    gene = row['gene']
    gene_list = gene.split(',')
    curr_input_df = effectDf[effectDf['gene'] == gene]
    can_starts, can_ends, can_sites = canonical_splice_site_lister(refDf, gene_list)
    if row['start'] in can_starts:
        starts_neutral += 1
    if row['start'] - 1 in can_starts:
        starts_minus1 += 1
    if row['start'] + 1 in can_starts:
        starts_plus1 += 1
    if row['start'] in can_ends:
        ends_neutral += 1
    if row['start'] - 1 in can_ends:
        ends_minus1 += 1
    if row['start'] + 1 in can_ends:
        ends_plus1 += 1     

print(starts_minus1)
print(starts_neutral)
print(starts_plus1)

print(ends_minus1)
print(ends_neutral)
print(ends_plus1)

print(131391108 in can_ends)
print(left_exon_finder(159789888,exons))
'''
'''
gene = 'WAC'
gene_list = gene.split(',')
curr_input_df = effectDf[effectDf['gene'] == gene]
can_starts, can_ends, can_sites = canonical_splice_site_lister(refDf, gene_list)
can_splices = canonical_splice_finder(refDf, gene_list)
exons = minimal_exon_finder(refDf, gene_list, can_sites)
regions = splice_region_lister(can_sites)
code_string = region_code_string_generator(can_sites,
                                           can_splices,
                                           exons,
                                           regions)
cryp_starts, cryp_ends = cryptic_site_finder(curr_input_df, can_sites)
curr_input_df['startClass'] = curr_input_df.apply(lambda x: start_site_annotator_GOT(x['start'], can_ends), axis=1)
curr_input_df['endClass'] = curr_input_df.apply(lambda x: end_site_annotator_GOT(x['end'], can_starts), axis=1)
curr_input_df['relStartExon'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, True, exons, can_sites), axis=1)
curr_input_df['relEndExon'] = curr_input_df.apply(lambda x: relevant_exon_finder(x, False, exons, can_sites), axis=1)
curr_input_df = alt_annotated_finder(curr_input_df, True)
curr_input_df = alt_annotated_finder(curr_input_df, False)
exon_starts = [x[0] for x in exons]
exon_ends = [x[1] for x in exons]
curr_input_df['startVerdict'] = curr_input_df.apply(lambda x: verdict_judger(x, True, exon_starts, exon_ends), axis=1)
curr_input_df['endVerdict'] = curr_input_df.apply(lambda x: verdict_judger(x, False, exon_starts, exon_ends), axis=1)
curr_input_df['start_min_exon_distance'] = curr_input_df.apply(lambda x: minimal_exon_distance_finder(x, True, exons), axis=1)
curr_input_df['end_min_exon_distance'] = curr_input_df.apply(lambda x: minimal_exon_distance_finder(x, False, exons), axis=1)
            
print([str(x)[2:] for x in can_sites])
print(code_string)
starts = curr_input_df['start']
ends = curr_input_df['end']
print([(str(x[0])[4:],str(x[1])[4:]) for x in can_splices])
print(curr_input_df[['start','end','start_min_exon_distance','end_min_exon_distance']])
print([str(x)[4:] for x in can_sites])
for i, row in curr_input_df.iterrows():
    print('index = %s' %i)
    print(row['relEndExon'])
    print(minimal_exon_distance_finder(row, False, exons))
   
print(len(code_string))
print(len(can_sites))
print(can_sites)
'''