
# Single Cell Splicing Analysis of ONT Reads

These tools allow for identification of differentially spliced transcripts in single cells from long read data. Full length cDNA produced from the 10X 3' single cell RNA sequencing kit is sequenced using oxford nanopore technologies (ONT) to get reads corresponding to full length transcripts. Prior to use of this pipeline, reads are aligned and cell barcodes and UMI's are identified in each read using the previously published pipeline (SiCeLoRe - https://www.nature.com/articles/s41467-020-17800-6). A tagged bam file containing tags for cell barcodes and UMI are then processed using a similar method to leafcutter (https://davidaknowles.github.io/leafcutter/) to call intron junctions found in each cell. Intron junctions are annotated as either alternative 3' or alternative 5' prior to differential transcript usage. 

Additionally, genotype status of mutation of interest (in this case SF3B1 in patients with myeloid displastic syndrome) was determined using the previously published method GoT (https://www.nature.com/articles/s41586-019-1367-0). Differential transcript usage can then be determined within each sample, comparing mutant and wild type cells, and can be further broken down by cell type. We also allow here an option to integrate across multiple single cell samples.

## Requirements 
- [SiCeLoRe](https://github.com/ucagenomix/sicelore)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://www.htslib.org/)
- [racon]
- Java 1.9 or higher

## Reference files and general setup

You will need to provide the following reference files: 

1. Minimap2 reference index in the `mmi` format.
In the below example `ref.fa` is the input fasta sequence and `ref.mmi` is the output index that will be used for mapping.

```
minimap2 -d ref.mmi ref.fa # indexing
```

2. Reference junctions bed file - To account for the high error rate in ONT, we utilize a reference of splice junctions identified from short read sequencing from a library without a splicing mutation. 
This is highly recommended as to identify the most accurate splice junction calls. 
Using the [splice aware alignment with STAR](https://github.com/alexdobin/STAR), we created a reference `.bed` file. 

3. Complementary short read data must be processed through [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) prior to using this pipeline to obtain the `posssorted_bam.bam` and `outs/filtered_feature_bc_matrix/barcodes.tsv` file. 

4. Short read data must also be processed through [IronThrone](https://github.com/dan-landau/IronThrone-GoT) to obtain genotyping results. 
A genotyping table containing the following columns is required: 

## GoT-Splice Pipeline overview 

The GoT-Splice pipeline is set up has three large steps, CB/UMI calling in long reads using SiCeLoRe, junction calling, and differential transcript usage. 
Currently, junction calling and differential transcript usage are two separate workflows that can be run independently, or they can be run as one larger pipeline within the GoT-splice workflow where SiCeLoRe is performed first prior to intron junction and differential transcript usage. 

We recommend first running SiCeLoRe following the user instructions. After which intron junction and then differential transcript usage can be performed. 

### Junction Calling in Single Cells

Junction calling is modelled after methods described in leafcutter, originally developed for identifying junctions in short read RNA-sequencing data. 

Step 1: Before calling junctions on your own sample, you must generate your own annotation reference files. This can be done on your own by using the function `gtf2leafcutter.pl()` from [leafcutter](). Or you may use our annotation references found in the [`junction_annotation/annotation_reference`](./junction_annotation/annotation_reference/) folder generated using Hg38. 

Step 2: Count intron junctions. Input bam must be bam file with cell barcode tags (BC) and UMI tags (U8) 
```
python count_introns_ONT.py \
  <path/to/input/bam> \
  <path_to_output_file>
```
Step 3: Creating metadata table and adding junction information to counts matrix. This step also includes the assigning of cluster IDs. A Cluster ID is an ID that is unique to all junctions that have the same 5p or the same 3p end. Junctions within the same cluster represent a group of junctions that have been alternatively spliced. 
```
Rscript junc_calling_script.R \
  <path/to/output_folder> \
  <sample_ID> \
  <path/to/annotation_reference>
```
Final ouptput: Counts matrix with each row as a junction and each column as a cell, metadata containing junction information, including gene, transcript ID, and chromosomal location. 

### Annotation of Intron Junctions**

Annotation of intron junctions identified in sample. Here uses the reference gtf to classify each 5' and 3' end of the intron. Output will add on additional columns to the metadata table including startClass and endClass, classifying the junction ends as the canonical (main) or alternative (not_main_3_prime/not_main_5_prime) end. 
```
python new.annotator_v2.py \
  <path/to/junction/calling/outputs> \
  <input datafile (output from previous step)> \
  <bed file output from intron junction calling> \
  <output file name> 
```

To run the junction annotation pipeline run the following: 

```
```

### Differential Transcript Usage

Strand adjustment of metadata file. Adds in columns to assign start/end as five prime or three prime ends of the gene for annotations. Output is metadata with new columns fivep_class and threep_class. 
```
Rscript strand_adjustment.R <metadata> <path/to/output>
```
#### Option A: Individual Patient

Run Differential transcript usage analysis. This requires the strand adjusted metadata matrix, full counts matrix from intron junction calling, and a table with cell barcodes, genotype information, and cell type assignment (produced from short read data analysis using GoT). From these inputs, junctions with no alternative sites and low coverage (n<5 - this threshold will be able to be altered by the user in a future version) will be filtered out prior to calculating an odds ratio for each junction in comparison to all other junctions with the same three prime or five prime end. Genotype assignments are then permuted x number of times (we recommend doing a test with 100-1000 but using at least 100,000 for making any final conclusions) and odds ratios are recalculated before determining the likelihood of the observed odds ratio being statistically significant. A final table for both Alt_5P and Alt_3P is output along with the log(odds ratio) for each junction and total observed reads across mutant and wild type cells. 

Format: Rscript Junc_Permute_logOR_alt3p_alt5p_within_cell_type.R <counts matrix> <genotype information> <strand adjusted metadata> <pattern/sample ID added to cell barcodes> <number of permutations> <path/to/output>
  
With a high number of junctions, it is likely that this may take a long time. If that is the case, we recommend splitting by cluster ID and then running for each group of smaller clusters a modified version that we have provided. See example below along with an example of how to run and submit on a slurm cluster in the examples folder. 

```
Rscript split_clusters_v2.R <counts> <genotype> <strand adjusted metadata> <pattern/ sample ID added to cell barcodes> <path/to/output>

Rscript split_JuncPermute_LogOR_perm_within_celltype_5p_3p.r <path/to/split/files> <genotype> <number of permtuations> <pattern> <output directory> <output file name> 
  
Rscript merge_final_output_ind_patient.R <path/to/split/output> <path/to/metadata> <path/to/output>

```
#### Option B: Combine Samples

If you have multiple samples that you would like to compare and look for differential transcript usage across all samples, we have also implemented a way to integrate the same analysis described above across our samples. Here, we ony keep clusters that are found in all samples and again only junctions that have n > 5 reads (this threshold will be able to be altered by the user in a future version). The output includes the individual odds ratio for each sample for each junction as well as a weighted odds ratio based on total cells found in each sample. There will be one p-value reported for each junction corresponding to the significance that this junction is differentially used across all samples. There is no minimum or maximum number of samples needed. Here you will need to move all files for all patients into one folder - i.e. all counts matrices should be in one folder and all genotype information should be in another folder. Again, to increase speed, we suggest that we split across clusters so have provided the scripts to do so. 

```
Rscript create_combined_metadata.R <path/to/metadata> <path/to/output> <patient.names> 

Rscript split_clusters_combined_patient.R <path/to/folder/with/counts> <path/to/folder/with/genotype> <path/to/combined/metadata> <path/to/output> <patient.names> <pattern/unique identifiers on cell barcodes>

Rscript split_JuncPermute_LogOR_combined_patient_within_celltype.R <path/to/split/files> <path/to/genotype> <number of permutations> <path/to/output> <output file name>
  
Rscript merge_final_output_comb_patient.R <path/to/outputs> <path/to/metadata> <patient.names> <path/to/output>
```

#### Option C: Within Cell Types/ Clusters

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

## Running the full GoT-Splice pipeline

To run as one combined pipeline on slurm cluster, use splice_pipeline.sh. 

```
sbatch splice_pipeline.sh \
  --fastq <full path to fastq, gzipped> \
  --short_read_files <directory to outs folder from Cell Ranger> \
  --sample_name <unique sample identifier> \
  --genotype_info <genotype_table> \
  --output_dir <directory to write all output files> \
  --pattern <pattern/sample ID on cell barcodes from integrated objects> \
  --ref_genome <path to reference minimap2 index> \
  --ref_junc_bed <path to reference bed file for splice junction correction> \
  --sicelore_dir <path to sicelore directory> \
  --minimap_dir <path to minimap2 directory>
```

:warning: This is only set up to run on a slurm HPC. 
It also assumes queue names of `pe2` and `bigmem`. If other queue names are used by your HPC those should be changed before running the workflow. 

This also relies on loading modules with the following commands: 

```
source /etc/profile.d/modules.sh
module load java/1.9
module load samtools
module load racon
```
These lines should either be changed or commented out before running and modules should be load in the correct format for your HPC. 