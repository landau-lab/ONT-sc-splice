#!/bin/bash

#SBATCH --job-name=sicelore
#SBATCH --mem-per-cpu=100gb
#SBATCH --ntasks=2
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ahawkins@nygenome.org
#SBATCH --output=stdout_%j.log
#SBATCH --error=stderror_%j.log

source /etc/profile.d/modules.sh
module load java/1.9
module load samtools

module load racon

PATH=$PATH:~/spoa/build/bin/
PATH=$PATH:/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2
echo $PATH

##### input arguments 

workdir=$1 #working directory where all output files will be 
fastq_seq1=$2
fastq_seq2=$3 #full path to fastq post polyA scan and split by barcode if necessary 
fastqdir=$4
sample_name=$5
parsedobj=$6
#barcodes=$8 #full path to barcodes file 

echo $sample_name.sam

#cd $fastqdir
#gunzip -c $fastq > "$sample_name".fastq
#cp "$sample_name".fastq $workdir

cd $workdir
#mkdir output_files
cd ./output_files
#mkdir combine_fastq
cd ./combine_fastq

combine_fastq="$sample_name".merge.fastq
cat $fastq_seq1 $fastq_seq2 > $combine_fastq
fastq_path="$workdir"/output_files/combine_fastq/"$combine_fastq"
cd ..

#mkdir sicelore_outputs
cd ./sicelore_outputs
#mkdir temp

#echo "scanning polyA" 
#java -jar /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeReadScanner-0.5.jar -i $fastq -o "$workdir"/output_files
#echo "polyA scan done" 

#polyAfastq=$workdir/output_files/sicelore_outputs/passed/"$sample_name"FWD.fastq ## make sure this includes workdir/passed/fastqnameFWD.fastq in name 

#echo "starting polyA passed read alignment" 
## step 1: align reads that pass PolyA scan 
#/gpfs/commons/home/ahawkins/software/minimap2/minimap2 -ax splice -uf --secondary=no --sam-hit-only -t 20 --junc-bed /gpfs/commons/home/gmullokandov/software/sicelore/Gencode/gencode.v31.hg38.junctions.bed /gpfs/commons/home/gmullokandov/software/ref_genome/GRCh38.p12.mmi $fastq_path > "$sample_name"_passed_reads.sam
#samtools view -Sb "$sample_name"_passed_reads.sam > "$sample_name"_passed_reads_unsorted.bam
#samtools sort "$sample_name"_passed_reads_unsorted.bam -o "$sample_name"_passed_reads_sorted.bam
#samtools index "$sample_name"_passed_reads_sorted.bam
#echo "alignment done" 

## step 2: split fastq files 
#echo "splitting fastq files" 
#/gpfs/commons/home/pchamely/software/fastp -i $fastq_path -Q -A --thread 1 --split_prefix_digits=3 --out1=sub.fastq --split=100
#echo "fastq files split" 

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
#SBATCH --mem=75g
#SBATCH --partition=pe2
#SBATCH --output=err_and_out/"$i"_stdout_%j.log
#SBATCH --error=err_and_out/"$i"_stderror_%j.log
#SBATCH --ntasks=5

source /etc/profile.d/modules.sh
module load samtools
module load bedtools
module load java/1.9

parsedobj=$6

#/gpfs/commons/home/gmullokandov/software/minimap2-2.17_x64-linux/minimap2 -ax splice -uf --MD --secondary=no --sam-hit-only -t 20 --junc-bed /gpfs/commons/home/pchamely/SequencingData/CH_dataset_02_27_20/ONT-GEX-CH305/smart_seq_alljunctions.bed /gpfs/commons/home/gmullokandov/software/ref_genome/GRCh38.p12.mmi "$i".sub.fastq > "$i".sub.sam
#samtools view -Sb "$i".sub.sam -o "$i".sub.unsorted.bam
#samtools sort "$i".sub.unsorted.bam -o "$i".sub.bam
#samtools index "$i".sub.bam

#java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddGeneNameTag I="$i".sub.bam O="$i".sub.GE.bam REFFLAT=/gpfs/commons/home/gmullokandov/software/ref_genome/gencode.v31.refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
#samtools index "$i".sub.GE.bam

java -jar -Xmx200g /gpfs/commons/home/ahawkins/sicelore/Jar/Sicelore-1.0.jar AddBamReadSequenceTag I="$i".sub.GE.bam O="$i".sub.GEUS.bam FASTQ="$i".sub.fastq
samtools index "$i".sub.GEUS.bam

java -jar -Xmx75g /gpfs/commons/home/ahawkins/sicelore/Jar/NanoporeBC_UMI_finder-1.0.jar -i "$i".sub.GEUS.bam -o "$i".sub.GEUS10xAttributes.bam -k $parsedobj --maxUMIfalseMatchPercent 2 --maxBCfalseMatchPercent 5 --ncpu 10 --logFile "$i".sub_NanoporeBC_UMI_finder.log

EOF

sicelore_jobids+=($(sbatch --job-name="$sample_name" "$sample_name"_"$i"_sicelore_2.sh))
done

### after this is complete then merge files and run second script to merge 

## step 4: Merge all files and generate consensus sequences 

#sbatch /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/sicelore_scripts/2.sicelore_part2/sicelore_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp

sic_part2=($(sbatch --dependency=singleton --job-name="$sample_name" /gpfs/commons/groups/landau_lab/ahawkins/MDS_ONT_splice/ONT_processing_pipeline/sicelore/2.sicelore_part2/sicelore_part2.sh "$workdir"/output_files/sicelore_outputs $sample_name $workdir/output_files/sicelore_outputs/temp))
