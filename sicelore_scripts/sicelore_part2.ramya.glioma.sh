#!/bin/bash

#SBATCH --job-name=sic_part2
#SBATCH --mem=500g
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=20
#SBATCH --output=err_and_out/stdout_%j.log
#SBATCH --error=err_and_out/stderror_%j.log

source /etc/profile.d/modules.sh
module load java/1.9
module load samtools
module load racon

##Executable paths for conensus step 
PATH=$PATH:/gpfs/commons/home/ahawkins/spoa/build/bin/
PATH=$PATH:/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2

workdir=$1
sample_name=$2
temp=$3

### Merge all *.sub.GEUS10xAttributes.bam files 

cd $workdir

#rm 100*

bamlist=$(for f in *.GEUS10xAttributes_umifound_.bam; do echo -n "I=$f " ; done)

java -jar -Xmx44g /gpfs/commons/home/pchamely/software/picard/picard.jar MergeSamFiles $bamlist ASSUME_SORTED=TRUE TMP_DIR=/scratch/tmp/ MAX_RECORDS_IN_RAM=100000000 OUTPUT="$sample_name".GEUS10xAttributes.umifound.bam VALIDATION_STRINGENCY=SILENT

samtools index "$sample_name".GEUS10xAttributes.umifound.bam

echo "all sam files merged" 

mkdir split_files 
mv *.sub* split_files/

mkdir split_files_scripts
mv *_sicelore_2.sh split_files_scripts/

## Generate consensus reads 
java -jar -Xmx500g /gpfs/commons/home/ahawkins/sicelore/Jar/dev_versions/Sicelore-1.0.dev.snp.jar ComputeConsensus I="$sample_name".GEUS10xAttributes.umifound.bam O="$sample_name"_consensus.fastq T=20 TMPDIR=$temp MAXREADS=1000
echo "Consensus sequences generated"

## Re-map consensus sequences to ref genome
/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2 -ax splice -uf --MD --secondary=no --sam-hit-only -t 20 --junc-bed /gpfs/commons/groups/landau_lab/rraviram/Splicing_ONT/MGH105_ACAligned.sortedByCoord.out.bed /gpfs/commons/home/gmullokandov/software/ref_genome/GRCh38.p12.mmi "$sample_name"_consensus.fastq > "$sample_name"_consensus.sam
samtools view -Sb "$sample_name"_consensus.sam -o "$sample_name"_consensus.unsorted.bam
samtools sort "$sample_name"_consensus.unsorted.bam -o "$sample_name"_consensus.sorted.bam
samtools index "$sample_name"_consensus.sorted.bam
echo "consensus sequences aligned" 

## add in cell BC and UMI tags 
java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddBamMoleculeTags I="$sample_name"_consensus.sorted.bam O="$sample_name"_consensus.sorted.tags.bam
samtools index "$sample_name"_consensus.sorted.tags.bam
echo "cellBC and UMI tags added" 

## add gene name tags
java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddGeneNameTag I="$sample_name"_consensus.sorted.tags.bam O="$sample_name"_consensus.sorted.tags.GE.bam REFFLAT=/gpfs/commons/home/gmullokandov/software/ref_genome/gencode.v31.refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index "$sample_name"_consensus.sorted.tags.GE.bam
echo "DONE!"



