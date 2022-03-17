#Python script to grab all junctions from a bam file of aligned reads

import pysam
import sys
import collections

def robust_get_tag(read, tag="ts"):
    try: 
        return(read.get_tag(tag))
    except KeyError: 
        return(".")

BAM_CREF_SKIP = 3
match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position

def find_introns_single_cell(read_iterator, strandtag = "GS", cell_barcode_tag = "BC"):

    res = {}
    cell_barcodes = set()
    umis = collections.Counter()
    for r in read_iterator:
        base_position = r.pos
        #chrom_strand = (r.reference_name, robust_get_tag(r, strandtag))
	#chrom_strand = robust_get_tag(r, strandtag)
        chrom_strand = (r.reference_name, "-" if r.is_reverse else "+")
        cell_barcode = robust_get_tag(r, cell_barcode_tag)
        umi = robust_get_tag(r, "U8")
        umis[umi] += 1
        if cell_barcode == ".": continue
        cell_barcodes.add(cell_barcode)
        for op, nt in r.cigartuples:
            if op in match_or_deletion: 
                base_position += nt
            elif op == BAM_CREF_SKIP: 
                junc_start = base_position
                base_position += nt
                if not chrom_strand in res: 
                    res[chrom_strand] = {}
                junc = (junc_start, base_position)
                if not junc in res[chrom_strand]: 
                    res[chrom_strand][junc] = collections.Counter()
                res[chrom_strand][junc][cell_barcode] += 1
    return(res,cell_barcodes,umis)

##MDS samples:
input_file_path = sys.argv[1]    
samfile = pysam.AlignmentFile(input_file_path, "rb")


res_sc,cell_barcodes,umis = find_introns_single_cell(samfile.fetch(), strandtag = "GS")

for chrom_strand,r in res_sc.items(): 
    print(chrom_strand,len(r))

import gzip

# for viewing in IGV
#with gzip.open("counts_pysam.bed.gz","w") as f: 
#    for chrom_strand,r in res.items(): 
#        for j,c in r.items():
#            if c >= 10: 
#                f.write( ("%s\t%i\t%i\t.\t%i\n" % (chrom_strand[0], j[0], j[1], c)).encode() )
                
cell_barcodes = list(cell_barcodes)

output_file_path = sys.argv[2]
with gzip.open(output_file_path, "wt", encoding='utf-8') as f: 
    f.write( "chrom\tstart\tend\tstrand" )
    for cb in cell_barcodes: 
        f.write( "\t%s" % cb )
    f.write("\n")
    for chrom_strand,r in res_sc.items(): 
        for junc,counts in r.items():
            f.write( "%s\t%i\t%i\t%s" % (chrom_strand[0], junc[0], junc[1], chrom_strand[1]) )
            for cb in cell_barcodes: 
                f.write( "\t%i" % counts[cb] )
            f.write("\n")
            
import numpy as np
umi_counts = list(umis.values())
print("Average value for the number of reads per umi:")
print(np.mean(umi_counts)) # knowles: 1.21 , Paulina: 1.611843008688349 					## Average value for the number of reads per umi's 
print("Median value for the number of reads per umi:")
print(np.median(umi_counts)) # knowles: 1.0, Paulina : 1.0   							## Median value for the number of reads per umi
print("Percent of umi's that have 2 supporting reads:")
print(np.sum(np.array(umi_counts)==2) / len(umi_counts)) # knowles: 0.139, Paulina: 0.16633508489002077  	## Percent of umi's that have 2 supporting reads
