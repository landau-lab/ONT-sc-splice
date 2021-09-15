#!/bin/bash

#SBATCH --job-name=generate.gene.cards
#SBATCH --partition=bigmem,pe2
#SBATCH --mem=100g

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RforPlots

sampleID='MDS.combined.2WT'
#crypticListPrefix='MDS.cryptic.junction.list'
#crypticListPrefix='cryptic.junction.list'
#crypticListPrefix='cryptic.junction.list.individual'
#crypticListPrefix='MDS.cryptic.junction.list.for.figure'
crypticListPrefix='MDS.cryptic.junction.list.BAX'


##------------------------------------------------- Adjust the crypitc junctions for the sashimiplots -------------------------------------------------##
Rscript bin/adjust.cryptic.junctions.R cryptic.junction.lists/"$crypticListPrefix".txt cryptic.junction.lists/"$crypticListPrefix".for.sashimi.txt

##-------------------------------------------------- For loop through file with cryptic site info -------------------------------------------------##
IFS=$'\n'
for junction_info in $(cat cryptic.junction.lists/"$crypticListPrefix".for.sashimi.txt);
do

gene=$(echo "$junction_info" | cut -d' ' -f1)
cryptic_site=$(echo "$junction_info" | cut -d' ' -f2)
cryptic_site_sashimi=$(echo "$junction_info" | cut -d' ' -f3)

echo "$gene"
echo "$cryptic_site"
echo "$cryptic_site_sashimi"

##-------------------------------------------------- Generate the dPSI vs Expression Plot and save as png --------------------------------------------------##
Rscript bin/generate.dPSI.expression.plot.R "$gene" "$cryptic_site" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/6.Comb_patients_2WT/5_read_threshold/MDS_P5_P6/ /gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/7.Comb_patients_ind_celltype_merge_counts_2WT/MDS_P5_P6/ "HSPC,IMP,MkP,MEP,EP,NP" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/4.Integrated_seurat_objs/ONT.seruat.objects/mds.ont.samples.integrated.rna.normscale.rds "Cell.Assignment" "Genotype_1WT" "$sampleID" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/MDSP5.MDSP6.all.intron.metadata.with.bulkPSI.ONT.tsv "bulk.psi" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/12.Bulk.PSI/output/CH502_DNMT3A.WT.only.all.intron.metadata.with.bulkPSI.ONT.tsv 

##-------------------------------------------------- Generate Sashimi plot and save as png --------------------------------------------------##
bash bin/generate.sashimi.plot.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/10.Bam_splitting/make.sashimi.plots/input.bams/MDS.untreated.combined.input.bams.bulk.and.ct.2UMI.tsv "$cryptic_site_sashimi" /gpfs/commons/groups/landau_lab/SF3B1_splice_project/15.Gene.Cards/sashimi.plots/"$sampleID"."$gene"."$cryptic_site".sashimi.png

##-------------------------------------------------- Merge the 2 putput files and save as pdf --------------------------------------------------##
Rscript bin/merge.plots.R dPSI.Expression.plots/"$sampleID"."$gene"."$cryptic_site".dPSI.Exp.plot.png sashimi.plots/"$sampleID"."$gene"."$cryptic_site".sashimi.png final.pdfs/"$sampleID"."$gene"."$cryptic_site".gene.card.pdf 

done
unset IFS 
