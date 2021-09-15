#!/bin/bash

#SBATCH --job-name=sicelore
#SBATCH --mem=500G
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderror_%j.log

source /etc/profile.d/modules.sh
module load java/1.9
module load samtools

module load racon

PATH=$PATH:/gpfs/commons/home/ahawkins/spoa/build/bin/
PATH=$PATH:/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2
echo $PATH

##### input arguments/ requirements 
## working directory must contain a folder labeled input_files 
## within input_files folder there are 3 folders
## 1.ONT_fastq -starts out empty but unzipped fastq is copied here 
## 2.short_read_files - contains barcodes.tsv (must be unzipped)  and possorted_genome_bam.bam and possorted_genome_bam.bam.bai from cellranger outs folder

workdir=$1 #working directory where all output files will be 
fastq=$2 #full path to fastq
fastqdir=$3
sample_name=$4
run_files=$5 ##location where all scripts are located
#genotype_info=$6
#pattern=$7
#nperm=$8
##barcodes=$8 #full path to barcodes file 

echo $sample_name

##moved simple checks (i.e. directory permissions, etc.) to the beginning
##so that basic mistakes are caught before complex computation
{
cd $workdir
} || {
echo "invalid directory: $workdir"
exit
}

if [ ! -d output_files ]
then
    {
    mkdir output_files
    echo "output_files folder created"
    } || {
    echo "unable to create output_files folder. Check permissions"
    exit
    }
else
    echo "output_files folder exists"
fi
    
cd ./output_files
if [ ! -d logs ]
then
    {
    mkdir logs
    echo "logs folder created"
    } || {
    echo "unable to create logs folder. Check permissions"
    exit
    }
else
    echo "output_files/logs folder exists"
fi

if [ ! -d sicelore_outputs ]
then
    {
    mkdir sicelore_outputs
    echo "sicelore_outputs folder created"
    } || {
    echo "unable to create sicelore_outputs folder. Check permissions"
    exit
    }
else
    echo "sicelore_outputs folder exists"
fi

cd ./sicelore_outputs
if [ ! -d temp ]
then
    {
    mkdir temp
    echo "temp folder created"
    } || {
    echo "unable to create temp folder. Check permissions"
    exit
    }
else
    echo "sicelore_outputs/temp folder exists"
fi

if [ ! -d  "$workdir"/input_files/2.short_read_files/10X_parsing_files ]
then
    {
    mkdir "$workdir"/input_files/2.short_read_files/10X_parsing_files
    echo "10X_parsing_files folder created"
    } || {
    echo "unable to create 10X_parsing_files folder. Check permissions"
    exit
    }
else
    echo "10X_parsing_files folder exists"
fi

if [ -f "$workdir"/input_files/2.short_read_files/barcodes.tsv ]
then
    barcodes=$workdir/input_files/2.short_read_files/barcodes.tsv
else
    echo "file missing: folder 2.short_read_files needs an unzipped barcodes.tsv file"
    exit
fi

if [ -f "$workdir"/input_files/2.short_read_files/possorted_genome_bam.bam ]
then
    bam="$workdir"/input_files/2.short_read_files/possorted_genome_bam.bam
else
    echo "file missing: folder 2.short_read_files needs a file possorted_genome_bam.bam"
    exit
fi

{
cd $fastqdir
} || {
echo "invalid directory: $fastqdir"
exit
}

##check if the fastq file already exists. If it does, don't unzip it again
##and proceed to next step. This is useful, if the pipeline crashes later
##and needs to be rerun
if [ ! -f "$workdir"/input_files/1.ONT_fastq/"$sample_name".fastq ]
then
    {
    gunzip -c $fastq > "$sample_name".fastq
    cp "$sample_name".fastq "$workdir"/input_files/1.ONT_fastq
    } || {
    echo "unable to copy fastq file to: $workdir/input_files/1.ONT_fastq. invalid fastq: $fastq"
    exit
}
else
    echo "fastq already exists in destination folder. running analysis on existing fastq"
fi

##polyA scanning, if completed will create a file polyADone.txt. This checks
##if the file exists and skips the step if it does. That way, if the pipeline
##is run multiple times, it can skip the polyA scanning step
if [ ! -f "$workdir"/output_files/sicelore_outputs/polyADone.txt ]
then
    {
    echo "scanning polyA"
    java -jar /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeReadScanner-0.5.jar -i "$workdir"/input_files/1.ONT_fastq/"$sample_name".fastq -o "$workdir"/output_files/sicelore_outputs
    echo "polyA scan done" > "$workdir"/output_files/sicelore_outputs/polyADone.txt
    } || {
    echo "polyA scan failed"
    exit
    }
else
    echo "polyA scan previously completed. running analyis with existing files"
fi

##10x parsing , if completed will create a file parsing10x.txt. This checks
##if the file exists and skips the step if it does. That way, if the pipeline
##is run multiple times, it can skip the 10x parsing step
##split this out so that it runs after parallel to other tasks
    
if [ ! -f "$workdir"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj ]
then
    echo "parsing 10X object"
    parsedobj="$workdir"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj
    ## ----- Parse for cell barcodes and UMIs ------------------
    ##set this as a separate job so other tasks can run simultaneously

    cd "$workdir"/input_files/2.short_read_files/10X_parsing_files
    cat > "$sample_name"_10X_parsing.sh <<EOF
#!/bin/bash

#SBATCH --job-name=10x_parsing
#SBATCH --mem-per-cpu=100G
#SBATCH --partition=bigmem
#SBATCH --output=stderror_%j.log
#SBATCH --error=stderror_%j.log
#SBATCH --ntasks=10

module load java/1.9

java -Xmx1000g -jar /gpfs/commons/home/ahawkins/sicelore/Jar/IlluminaParser-1.0.jar --inFileIllumina $bam --tsv $barcodes --outFile $parsedobj --cellBCflag CB --umiFlag UB --geneFlag GN

EOF
    
    sbatch --job-name="$sample_name"_10X_parsing "$sample_name"_10X_parsing.sh

else
    parsedobj="$workdir"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj
    echo "10x object parsing previously completed. running analyis with existing files"
fi

##fastq splitting
##fastp outputs a file fastp.json. If this file exists in the folder,
##it has already run, and can be skipped
cd $workdir/output_files/sicelore_outputs

if [ ! -f "$workdir"/output_files/sicelore_outputs/fastp.json ]
then
    if [ -f "$workdir"/output_files/sicelore_outputs/passed/"$sample_name"FWD.fastq ]
    then
        polyAfastq="$workdir"/output_files/sicelore_outputs/passed/"$sample_name"FWD.fastq ## make sure this includes workdir/passed/fastqnameFWD.fastq in name
    else
        echo "no passed sicelore outputs"
        exit
    fi

    echo "splitting fastq files"
    {
    /gpfs/commons/home/pchamely/software/fastp -i $polyAfastq -Q -A --thread 1 --split_prefix_digits=3 --out1=sub.fastq --split=100
    echo "fastq files split"
    } || {
    echo "fastq file splitting failed"
    exit
    }
else
    echo "fastq file previously split. running analysis on previous files"
fi

if [ ! -d  "$workdir"/output_files/sicelore_outputs/err_and_out ]
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

sicelore_jobids=()

## step 3: Alignment, UMI, and CB matching for all files
## perform this step only for non-empty files
## do not repeat it, if it was completed previously
jobs=$(seq -w 1 100)

for i in ${jobs[@]};
do

if [ -f "$workdir"/output_files/sicelore_outputs/"$i".sub.fastq ]
then
    if [ ! -f "$workdir"/output_files/sicelore_outputs/"$i".sub.GEUS10xAttributes.bam ]
    then
        filename="$workdir"/output_files/sicelore_outputs/"$i".sub.fastq
        size=$(stat -c %s $filename)
        if (( $size>0 ))
        then
            cat > "$sample_name"_"$i"_sicelore_2.sh <<EOF
#!/bin/bash

#SBATCH --job-name=Sicelore_2
#SBATCH --mem=200G
#SBATCH --partition=bigmem,pe2
#SBATCH --output=err_and_out/"$i"_stdout_%j.log
#SBATCH --error=err_and_out/"$i"_stderror_%j.log
#SBATCH --cpus-per-task=5

source /etc/profile.d/modules.sh
module load samtools
module load bedtools
module load java/1.9

sample_name=$4

parsedobj="$1"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj
    
/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2 -ax splice -uf --MD --secondary=no --sam-hit-only -t 20 --junc-bed /gpfs/commons/groups/landau_lab/schypertribe/bulk/splice_junctions/H1_MPP2_B3_junctions.bed /gpfs/commons/groups/landau_lab/genomes/Mouse/refdata-gex-mm10-2020-A/fasta/genome.mmi "$i".sub.fastq > "$i".sub.sam
samtools view -Sb "$i".sub.sam -o "$i".sub.unsorted.bam
samtools sort "$i".sub.unsorted.bam -o "$i".sub.bam
samtools index "$i".sub.bam

java -jar -Xmx100g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddGeneNameTag I="$i".sub.bam O="$i".sub.GE.bam REFFLAT=/gpfs/commons/groups/landau_lab/genomes/Mouse/refdata-gex-mm10-2020-A/fasta/gencode.reflat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index "$i".sub.GE.bam

java -jar -Xmx100g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddBamReadSequenceTag I="$i".sub.GE.bam O="$i".sub.GEUS.bam FASTQ="$i".sub.fastq
samtools index "$i".sub.GEUS.bam

java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeBC_UMI_finder-1.0.jar -i "$i".sub.GEUS.bam -o "$i".sub.GEUS10xAttributes.bam -k $parsedobj --maxUMIfalseMatchPercent 2 --maxBCfalseMatchPercent 5 --ncpu 5 --polyawin 150 --logFile "$i".sub_NanoporeBC_UMI_finder.log

EOF

        sicelore_jobids+=($(sbatch --job-name="$sample_name" "$sample_name"_"$i"_sicelore_2.sh))
        fi
    else
    echo "analysis previously run on $i.sub.fastq. Using previous outputs"
    fi
fi
done
echo $sicelore_jobids


### after this is complete then merge files and run second script to merge 

## step 4: Merge all files and generate consensus sequences 

##sbatch  /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/sicelore_scripts/2.sicelore_part2/sicelore_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp

#sic_part2=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/sicelore_scripts/sicelore_mouse_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp))

cd $workdir

## step 5: Junction calling and annotation 
##sbatch --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin

junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline_mouse.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin))

## step 6: Differential transcript usage (with permutations within cell types)

##diff_transcript=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/diff_transcript_usage_scripts/diff_transcript_pipeline.sh "$workdir"/output_files $run_files "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_all.introns.info.w.primary.annotations.txt "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  $genotype_info $pattern $sample_name $nperm))
