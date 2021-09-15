#!/bin/bash

#SBATCH --job-name=sicelore
#SBATCH --mem=500gb
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderror_%j.log

source /etc/profile.d/modules.sh
module load java/1.9
module load samtools

module load racon

PATH=$PATH:~/spoa/build/bin/
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
genotype_info=$6
pattern=$7
nperm=$8
#barcodes=$8 #full path to barcodes file 

echo $sample_name

cd $fastqdir
gunzip -c $fastq > "$sample_name".fastq
cp "$sample_name".fastq "$workdir"/input_files/1.ONT_fastq

cd $workdir
mkdir output_files
cd ./output_files
mkdir logs

mkdir sicelore_outputs
cd ./sicelore_outputs
mkdir temp

#echo "scanning polyA" 
#java -jar /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeReadScanner-0.5.jar -i "$workdir"/input_files/1.ONT_fastq/"$sample_name".fastq -o "$workdir"/output_files/sicelore_outputs
#echo "polyA scan done" 

barcodes="$workdir"/input_files/2.short_read_files/barcodes.tsv
bam="$workdir"/input_files/2.short_read_files/possorted_genome_bam.bam
parsedobj="$workdir"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj


#echo "parsing 10X object" 
 # ----- Parse for cell barcodes and UMIs ------------------
#java -Xmx500g -jar /gpfs/commons/home/ahawkins/sicelore/Jar/IlluminaParser-1.0.jar --inFileIllumina $bam --tsv $barcodes --outFile $parsedobj --cellBCflag CB --umiFlag UB --geneFlag GN

polyAfastq=$workdir/output_files/sicelore_outputs/passed/"$sample_name"FWD.fastq ## make sure this includes workdir/passed/fastqnameFWD.fastq in name 

## split fastq files 
echo "splitting fastq files" 
/gpfs/commons/home/pchamely/software/fastp -i $polyAfastq -Q -A --thread 1 --split_prefix_digits=3 --out1=sub.fastq --split=100
echo "fastq files split" 

jobs=$(seq -w 1 100)

#jobs=('0001' '0002' '0003' '0004' '0005' '0006' '0007' '0008' '0009' '0010' '0011' '0012' '0013' '0014' '0015' '0016' '0017' '0018' '0019' '0020')

mkdir err_and_out

sicelore_jobids=()

# step 3: Alignment, UMI, and CB matching for all files 
for i in ${jobs[@]};
do

cat > "$sample_name"_"$i"_sicelore_2.sh <<EOF
#!/bin/bash

#SBATCH --job-name=Sicelore_2
#SBATCH --mem-per-cpu=75gb
#SBATCH --partition=pe2
#SBATCH --output=err_and_out/"$i"_stdout_%j.log
#SBATCH --error=err_and_out/"$i"_stderror_%j.log
#SBATCH --ntasks=5

source /etc/profile.d/modules.sh
module load samtools
module load bedtools
module load java/1.9

sample_name=$4
parsedobj="$1"/input_files/2.short_read_files/"$sample_name"_parsed_for_Nanopore.obj

/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2 -ax splice -uf --MD --secondary=no --sam-hit-only -t 20 --junc-bed /gpfs/commons/home/pchamely/SequencingData/CH_dataset_02_27_20/ONT-GEX-CH305/smart_seq_alljunctions.bed /gpfs/commons/home/gmullokandov/software/ref_genome/GRCh38.p12.mmi "$i".sub.fastq > "$i".sub.sam
samtools view -Sb "$i".sub.sam -o "$i".sub.unsorted.bam
samtools sort "$i".sub.unsorted.bam -o "$i".sub.bam
samtools index "$i".sub.bam

java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddGeneNameTag I="$i".sub.bam O="$i".sub.GE.bam REFFLAT=/gpfs/commons/home/gmullokandov/software/ref_genome/gencode.v31.refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index "$i".sub.GE.bam

java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddBamReadSequenceTag I="$i".sub.GE.bam O="$i".sub.GEUS.bam FASTQ="$i".sub.fastq
samtools index "$i".sub.GEUS.bam

java -jar -Xmx300g /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeBC_UMI_finder-1.0.jar -i "$i".sub.GEUS.bam -o "$i".sub.GEUS10xAttributes.bam -k $parsedobj --maxUMIfalseMatchPercent 2 --maxBCfalseMatchPercent 5 --ncpu 10 --logFile "$i".sub_NanoporeBC_UMI_finder.log

EOF

sicelore_jobids+=($(sbatch --job-name="$sample_name" "$sample_name"_"$i"_sicelore_2.sh))
done

### after this is complete then merge files and run second script to merge 

## step 4: Merge all files and generate consensus sequences 

#sbatch /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/sicelore_scripts/2.sicelore_part2/sicelore_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp

sic_part2=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/sicelore_scripts/sicelore_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp))

cd $workdir

## step 5: Junction calling and annotation 
#sbatch --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin

junc_calling=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/leafcutter_scripts/junc_calling_pipeline.sh "$workdir"/output_files "$workdir"/output_files/sicelore_outputs/"$sample_name"_consensus.sorted.tags.GE.bam $sample_name "$run_files"/bin))

## step 6: Differential transcript usage (with permutations within cell types)

diff_transcript=($(sbatch --dependency=singleton --job-name="$sample_name" "$run_files"/diff_transcript_usage_scripts/diff_transcript_pipeline.sh "$workdir"/output_files $run_files "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_all.introns.info.w.primary.annotations.txt "$workdir"/output_files/leafcutter_outputs/"$sample_name"_output/"$sample_name"_perind_numbers.counts.txt  $genotype_info $pattern $sample_name $nperm))

