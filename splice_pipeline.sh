#!/bin/bash

#SBATCH --job-name=sicelore
#SBATCH --mem=500G
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderror_%j.log

######################################################################################
# This script is used to submit the GoT-ONT splice pipeline for an individual sample 

# The output will include differential junction usage between all MUT and WT cells 
# and between all MUT and WT cells across all specified cell types 

# Any of the following parameters can be invoked at the command line: 

# fastq : Full path to FASTQ file from ONT
# short_read_files : Full path to folder with short read files
# sample_name : Name of sample being processed 
# scripts_dir : Path to where the GoT ONT pipeline scripts are stored 
# genotype_info : Table with genotype calls for each cell. 
# pattern : 
# output_dir : Where to store all output files  
# nperm : Number of permutations to perform for permutation testing. Default is 100000. 
# min_reads : Minimum number of reads a junction must have across all cells to be 
#   included in differential transcript usage. Default is 5.
# sicelore_dir : Directory where sicelore is located. 
# minimap_dir : Directory where minimap2 is located. 

######################################################################################

############ Command line arguments and load modules #################################

# set defaults 
scripts_dir="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/ONT_Splice_Pipeline/"
nperm=100000
min_reads=5
sicelore_dir="/gpfs/commons/home/ahawkins/sicelore/"
minimap_dir="/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2"
ref_junc_bed="/gpfs/commons/home/pchamely/SequencingData/CH_dataset_02_27_20/ONT-GEX-CH305/smart_seq_alljunctions.bed"
ref_genome="/gpfs/commons/home/gmullokandov/software/ref_genome/GRCh38.p12.mmi"

# grab arguments from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# load modules 
source /etc/profile.d/modules.sh
module load java/1.9
module load samtools
module load racon

# add minimap to path 
PATH=$PATH:/gpfs/commons/home/ahawkins/spoa/build/bin/
PATH=$PATH:$minimap_dir

echo Completing GoT-ONT Pipeline for $sample_name...

############ Error Checks #########################################################

##moved simple checks (i.e. directory permissions, etc.) to the beginning
##so that basic mistakes are caught before complex computation
{
cd $output_dir
} || {
echo "invalid directory: $output_dir"
exit
}

# if short read files are not provided, build the short read files directory 
# after checking for the output directory 
short_read_files=${short_read_files:-"$output_dir/input_files/2.short_read_files"}

# build path for sicelore outputs 
sicelore_outputs="$output_dir/output_files/sicelore_outputs"


# check that output_files folder exists and create sicelore_outputs, logs, and temp folders 
if [ ! -d output_files ]
then
    {
    mkdir -p output_files
    echo "output_files folder created"
    } || {
    echo "unable to create output_files folder. Check permissions"
    exit
    }
else
    echo "output_files folder exists"
fi
    
cd $output_dir/output_files
if [ ! -d logs ]
then
    {
    mkdir -p logs
    echo "logs folder created"
    } || {
    echo "unable to create logs folder. Check permissions"
    exit
    }
else
    echo "output_files/logs folder exists"
fi

if [ ! -d $sicelore_outputs ]
then
    {
    mkdir -p $sicelore_outputs
    echo "sicelore_outputs folder created"
    } || {
    echo "unable to create sicelore_outputs folder. Check permissions"
    exit
    }
else
    echo "sicelore_outputs folder exists"
fi

cd $sicelore_outputs
if [ ! -d temp ]
then
    {
    mkdir -p temp
    echo "temp folder created"
    } || {
    echo "unable to create temp folder. Check permissions"
    exit
    }
else
    echo "sicelore_outputs/temp folder exists"
fi

# create 10X parsing folder and check that barcodes/posssorted bam are provided in the input folder 
if [ ! -d  $short_read_files/10X_parsing_files ]
then
    {
    mkdir $short_read_files/10X_parsing_files
    echo "10X_parsing_files folder created"
    } || {
    echo "unable to create 10X_parsing_files folder. Check permissions"
    exit
    }
else
    echo "10X_parsing_files folder exists"
fi

if [ -f $short_read_files/barcodes.tsv ]
then
    barcodes=$short_read_files/barcodes.tsv
else
    echo "file missing: folder $short_read_files needs an unzipped barcodes.tsv file"
    exit
fi

if [ -f $short_read_files/possorted_genome_bam.bam ]
then
    bam="$short_read_files/possorted_genome_bam.bam"
else
    echo "file missing: folder $short_read_files needs a file possorted_genome_bam.bam"
    exit
fi

######### Unzip FASTQ file ######################################################################

# check that fastq file exists
if [ ! -f $fastq ]
then 
    echo $fastq does not exist. 
    exit 
fi 

# check that folder to copy over fastq file exists otherwise create it 
fastqdir=$output_dir/input_files/1.ONT_fastq
if [ ! -d $fastqdir && ! -f "$sample_name".fastq ]
then 
    mkdir -p $fastqdir
else 
    echo $fastqdir already exists.

##check if the fastq file already exists. If it does, don't unzip it again
##and proceed to next step. This is useful, if the pipeline crashes later
##and needs to be rerun
if [ ! -f $fastqdir/$sample_name.fastq ]
then
    {
    gunzip -c $fastq > "$sample_name".fastq
    cp $sample_name.fastq $fastqdir
    } || {
    echo "invalid fastq: $fastq"
    exit
    }
else
    echo "fastq already exists in original folder. running analysis on existing fastq"
fi

########### Poly A Scanning ##################################################################

##polyA scanning, if completed will create a file polyADone.txt. This checks
##if the file exists and skips the step if it does. That way, if the pipeline
##is run multiple times, it can skip the polyA scanning step
if [ ! -f "$sicelore_outputs/polyADone.txt" ]
then
    {
    echo "scanning polyA"
    # run sicelore polyA scanner 
    java -jar $sicelore_dir/Jar/NanoporeReadScanner-0.5.jar \
      -i "$fastqdir"/"$sample_name".fastq \
      -o "$output_dir"/output_files/sicelore_outputs

    echo "polyA scan done" > $sicelore_outputs/polyADone.txt

    } || {

    echo "polyA scan failed"
    exit

    }
else

    echo "polyA scan previously completed. running analyis with existing files"

fi

########### 10X UMI/CB Parsing ##############################################################

##10x parsing , if completed will create a file parsing10x.txt. This checks
##if the file exists and skips the step if it does. That way, if the pipeline
##is run multiple times, it can skip the 10x parsing step
##split this out so that it runs after parallel to other tasks
    
if [ ! -f "$short_read_files"/"$sample_name"_parsed_for_Nanopore.obj ]
then
    echo "parsing 10X object"

    parsedobj="$short_read_files"/"$sample_name"_parsed_for_Nanopore.obj

    ## ----- Parse for cell barcodes and UMIs ------------------
    ##set this as a separate job so other tasks can run simultaneously

    cd "$short_read_files"/10X_parsing_files
    cat > "$sample_name"_10X_parsing.sh <<EOF

#!/bin/bash

#SBATCH --job-name=10x_parsing
#SBATCH --mem-per-cpu=100G
#SBATCH --partition=pe2
#SBATCH --output=stderror_%j.log
#SBATCH --error=stderror_%j.log
#SBATCH --ntasks=5

module load java/1.9

java -Xmx500g -jar $sicelore_dir/Jar/IlluminaParser-1.0.jar \
  --inFileIllumina $bam \
  --tsv $barcodes \
  --outFile $parsedobj \
  --cellBCflag CB \
  --umiFlag UB \
  --geneFlag GN

EOF
    
    sbatch --job-name="$sample_name"_10X_parsing "$sample_name"_10X_parsing.sh

else
    parsedobj="$$short_read_files"/"$sample_name"_parsed_for_Nanopore.obj
    echo "10x object parsing previously completed. running analyis with existing files"
fi

############## FASTQ splitting ##############################################################

##fastq splitting
##fastp outputs a file fastp.json. If this file exists in the folder,
##it has already run, and can be skipped

cd $sicelore_outputs

if [ ! -f "$sicelore_outputs"/fastp.json ]
then

# check to make sure that polyA scanning is complete 
    if [ -f "$sicelore_outputs"/passed/"$sample_name"FWD.fastq ]
    then
        polyAfastq="$sicelore_outputs"/passed/"$sample_name"FWD.fastq ## make sure this includes output_dir/passed/fastqnameFWD.fastq in name
    else
        echo "no passed sicelore outputs"
        exit
    fi

    echo "splitting fastq files"
    {

    /gpfs/commons/home/pchamely/software/fastp \
      -i $polyAfastq \
      -Q \
      -A \
      --thread 1 \
      --split_prefix_digits=3 \
      --out1=sub.fastq \
      --split=100

    echo "fastq files split"
    } || {
    echo "fastq file splitting failed"
    exit
    }
else
    echo "fastq file previously split. running analysis on previous files"
fi

if [ ! -d  "$sicelore_outputs"/err_and_out ]
then
    {
    mkdir err_and_out
    echo "err_and_out folder created"
    } || {
    echo "unable to create err_and_out folder. Check permissions"
    exit
    }
else
    echo "err_and_out folder exists"
fi

####### Alignment, UMI, CB Matching ###############################################################

sicelore_jobids=()

## step 3: Alignment, UMI, and CB matching for all files
## perform this step only for non-empty files
## do not repeat it, if it was completed previously

jobs=$(seq -w 1 100)

for i in ${jobs[@]};
do

if [ -f "$sicelore_outputs"/"$i".sub.fastq ]
then
    if [ ! -f "$sicelore_outputs"/"$i".sub.GEUS10xAttributes.bam ]
    then
        filename="$sicelore_outputs"/"$i".sub.fastq
        size=$(stat -c %s $filename)
        if (( $size>0 ))
        then
            cat > "$sample_name"_"$i"_sicelore_2.sh <<EOF
#!/bin/bash

#SBATCH --job-name=Sicelore_2
#SBATCH --mem=300G
#SBATCH --partition=bigmem,pe2
#SBATCH --output=err_and_out/"$i"_stdout_%j.log
#SBATCH --error=err_and_out/"$i"_stderror_%j.log
#SBATCH --cpus-per-task=5

source /etc/profile.d/modules.sh
module load samtools
module load bedtools
module load java/1.9

sample_name=$sample_name
parsedobj="$short_read_files"/"$sample_name"_parsed_for_Nanopore.obj

$minimap_dir \
  -ax splice \
  -uf \
  --MD \
  --secondary=no \
  --sam-hit-only \
  -t 20 \
  --junc-bed $ref_junc_bed \
  $ref_genome \
  "$i".sub.fastq > "$i".sub.sam

samtools view -Sb "$i".sub.sam -o "$i".sub.unsorted.bam
samtools sort "$i".sub.unsorted.bam -o "$i".sub.bam
samtools index "$i".sub.bam

java -jar -Xmx100g $sicelore_dir/Jar/Sicelore-1.0.jar AddGeneNameTag \
  I="$i".sub.bam \
  O="$i".sub.GE.bam \
  REFFLAT=/gpfs/commons/home/gmullokandov/software/ref_genome/gencode.v31.refFlat \
  GENETAG=GE \
  ALLOW_MULTI_GENE_READS=true \
  USE_STRAND_INFO=true \
  VALIDATION_STRINGENCY=SILENT

samtools index "$i".sub.GE.bam

java -jar -Xmx100g $sicelore_dir/Jar/Sicelore-1.0.jar AddBamReadSequenceTag \
  I="$i".sub.GE.bam \
  O="$i".sub.GEUS.bam \
  FASTQ="$i".sub.fastq

samtools index "$i".sub.GEUS.bam

java -jar -Xmx300g $sicelore_dir/Jar/NanoporeBC_UMI_finder-1.0.jar \
  -i "$i".sub.GEUS.bam -o \
  "$i".sub.GEUS10xAttributes.bam \
  -k $parsedobj \
  --maxUMIfalseMatchPercent 2 \
  --maxBCfalseMatchPercent 5 \
  --ncpu 5 \
  --logFile "$i".sub_NanoporeBC_UMI_finder.log

EOF

        sicelore_jobids+=($(sbatch --job-name="$sample_name" "$sample_name"_"$i"_sicelore_2.sh))
        fi
    else
    echo "analysis previously run on $i.sub.fastq. Using previous outputs"
    fi
fi
done

### after this is complete then merge files and run second script to merge 

########## Consensus Generation ######################################################

## step 4: Merge all files and generate consensus sequences 

sic_part2=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/sicelore_scripts/sicelore_part2.sh "$output_dir"/output_files/sicelore_outputs $sample_name $output_dir/output_files/sicelore_outputs/temp))

cd $output_dir

############# Junction calling #########################################################

## step 5: Junction calling and annotation 
#sbatch --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$output_dir"/output_files "$output_dir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin

junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/junction_annotation/junc_calling_pipeline.w.exon.skip.sh \
  --output_location "$output_dir"/output_files \
  --input_bam "$sicelore_outputs"/"$sample_name"_consensus.sorted.tags.GE.bam \
  --sample $sample_name \
  --scripts_dir "$scripts_dir"/junction_annotation))

############## Diff Transcript Usage ####################################################

## step 6: Differential transcript usage (with permutations within cell types)

diff_transcript=($(sbatch --dependency=singleton --job-name="$sample_name" "$scripts_dir"/diff_transcript_usage/diff_transcript_pipeline.sh \
  "$output_dir"/output_files \
  $run_files \
  "$output_dir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_all.introns.info.w.primary.annotations.txt \
  "$output_dir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  \
  $genotype_info \
  $pattern \
  $sample_name \
  $nperm))
'
