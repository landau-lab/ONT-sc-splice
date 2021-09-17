
# **Single Cell Splicing Analysis of ONT Reads**

These tools allow for identification of differentially spliced transcripts in single cells from long read data. Full length cDNA produced from the 10X 3' single cell RNA sequencing kit is sequenced using oxford nanopore technologies (ONT) to get reads corresponding to full length transcripts. Prior to use of this pipeline, reads are aligned and cell barcodes and UMI's are identified in each read using the previously published pipeline (SiCeLoRe - https://www.nature.com/articles/s41467-020-17800-6). A tagged bam file containing tags for cell barcodes and UMI are then processed using a similar method to leafcutter (https://davidaknowles.github.io/leafcutter/) to call intron junctions found in each cell. Intron junctions are annotated as either alternative 3' or alternative 5' prior to differential transcript usage. 

Additionally, genotype status of mutation of interest (in this case SF3B1 in patients with myeloid displastic syndrome) was determined using the previously published method GoT (https://www.nature.com/articles/s41586-019-1367-0). Differential transcript usage can then be determined within each sample, comparing mutant and wild type cells, and can be further broken down by cell type. We also allow here an option to integrate across multiple single cell samples.

**Dependencies**
* [SiCeLoRe](https://github.com/ucagenomix/sicelore)
* R
* Samtools
* Python

**Installation** 

Clone the files. 

To run as one combined pipeline on slurm cluster, use splice_pipeline.sh. 

```
sbatch splice_pipeline.sh <output directory> <full path to fastq, gzipped> <directory where fastq is stored> <sample_name> <path to all scripts> <genotype_table> <pattern/sample ID on cell barcode> <number of permutations> 
```

# **Intron Junction Calling in Single Cells**

Intron junction calling is modelled after methods described in leafcutter, originally developed for identifying junctions in short read RNA-sequencing data. 

Step 1: Before calling introns on your own sample, you must generate your own annotation reference files. This can be done on your own by using the function gtf2leafcutter.pl found on leafcutter's page. Or you may use our annotation references found in the annotation_reference folder generated using Hg38. 

Step 2: Count intron junctions. Input bam must be bam file with cell barcode tags (BC) and UMI tags (U8) 
```
python count_introns_ONT.py <path/to/input/bam> <path_to_output_file>
```
Step 3: Creating metadata table and adding junction information to counts matrix. This step also includes the assigning of cluster IDs. A Cluster ID is an ID that is unique to all junctions that have the same 5p or the same 3p end. Junctions within the same cluster represent a group of junctions that have been alternatively spliced. 
```
Rscript junc_calling_script.R <path/to/output_folder> <sample_ID> <path/to/annotation_reference>
```
Final ouptput: Counts matrix with each row as a junction and each column as a cell, metadata containing junction information, including gene, transcript ID, and chromosomal location. 

# **Annotation of Intron Junctions**

Annotation of intron junctions identified in sample. Here uses the reference gtf to classify each 5' and 3' end of the intron. Output will add on additional columns to the metadata table including startClass and endClass, classifying the junction ends as the canonical (main) or alternative (not_main_3_prime/not_main_5_prime) end. 
```
python new.annotator_v2.py <path/to/junction/calling/outputs> <input datafile (output from previous step)> <bed file output from intron junction calling> <output file name> 
```

# **Differential Transcript Usage** 

Strand adjustment of metadata file. Adds in columns to assign start/end as five prime or three prime ends of the gene for annotations. Output is metadata with new columns fivep_class and threep_class. 
```
Rscript strand_adjustment.R <metadata> <path/to/output>
```
## **Option A: Individual Patient**

Run Differential transcript usage analysis. This requires the strand adjusted metadata matrix, full counts matrix from intron junction calling, and a table with cell barcodes, genotype information, and cell type assignment (produced from short read data analysis using GoT). From these inputs, junctions with no alternative sites and low coverage (n<5 - this threshold will be able to be altered by the user in a future version) will be filtered out prior to calculating an odds ratio for each junction in comparison to all other junctions with the same three prime or five prime end. Genotype assignments are then permuted x number of times (we recommend doing a test with 100-1000 but using at least 100,000 for making any final conclusions) and odds ratios are recalculated before determining the likelihood of the observed odds ratio being statistically significant. A final table for both Alt_5P and Alt_3P is output along with the log(odds ratio) for each junction and total observed reads across mutant and wild type cells. 

Format: Rscript Junc_Permute_logOR_alt3p_alt5p_within_cell_type.R <counts matrix> <genotype information> <strand adjusted metadata> <pattern/sample ID added to cell barcodes> <number of permutations> <path/to/output>
  
With a high number of junctions, it is likely that this may take a long time. If that is the case, we recommend splitting by cluster ID and then running for each group of smaller clusters a modified version that we have provided. See example below along with an example of how to run and submit on a slurm cluster in the examples folder. 

```
Rscript split_clusters_v2.R <counts> <genotype> <strand adjusted metadata> <pattern/ sample ID added to cell barcodes> <path/to/output>

Rscript split_JuncPermute_LogOR_perm_within_celltype_5p_3p.r <path/to/split/files> <genotype> <number of permtuations> <pattern> <output directory> <output file name> 
  
Rscript merge_final_output_ind_patient.R <path/to/split/output> <path/to/metadata> <path/to/output>

```
## **Option B: Combine Samples**

If you have multiple samples that you would like to compare and look for differential transcript usage across all samples, we have also implemented a way to integrate the same analysis described above across our samples. Here, we ony keep clusters that are found in all samples and again only junctions that have n > 5 reads (this threshold will be able to be altered by the user in a future version). The output includes the individual odds ratio for each sample for each junction as well as a weighted odds ratio based on total cells found in each sample. There will be one p-value reported for each junction corresponding to the significance that this junction is differentially used across all samples. There is no minimum or maximum number of samples needed. Here you will need to move all files for all patients into one folder - i.e. all counts matrices should be in one folder and all genotype information should be in another folder. Again, to increase speed, we suggest that we split across clusters so have provided the scripts to do so. 

```
Rscript create_combined_metadata.R <path/to/metadata> <path/to/output> <patient.names> 

Rscript split_clusters_combined_patient.R <path/to/folder/with/counts> <path/to/folder/with/genotype> <path/to/combined/metadata> <path/to/output> <patient.names> <pattern/unique identifiers on cell barcodes>

Rscript split_JuncPermute_LogOR_combined_patient_within_celltype.R <path/to/split/files> <path/to/genotype> <number of permutations> <path/to/output> <output file name>
  
Rscript merge_final_output_comb_patient.R <path/to/outputs> <path/to/metadata> <patient.names> <path/to/output>
```

## **Option C: Within Cell Types/ Clusters** 

Finally, using information defined from short read illumina sequencing, you can identify differentially used transcripts across cell types. Here, we integrate across samples and across cell types in order to identify key transcripts that are differentially spliced between mutant and wild type cells in various cell types. Permutations happen within each sample and within each cell type. The output table contains a combined table with weighted odds ratio for each junction found in each cell type with one p-value being reported for the likelihood that junction is significantly differentially used in each cell type. This can also be applied to clusters within a sample and does not need to be within cell types. 

To run individual patient within each cell type, run the following: 
```
Rscript split_clusters_combined_patient.R <path/to/folder/with/counts> <path/to/folder/with/genotype> <path/to/combined/metadata> <path/to/output> <patient.names> <pattern/unique identifiers on cell barcodes>

Rscript split_JuncPermute_ind_patient_mut_wt_all_celltypes.R <path/to/split/files> <path/to/genotype> <number of permutations> <path/to/output> <output file name>
```

To run combined patient within each cell type, run the following: 
```
Rscript JuncPermute_LogOR_combined_patient_celltype.R <path/to/counts> <path/to/genotype> <path/to/metadata> <number of permutations> <path/to/output>
```
