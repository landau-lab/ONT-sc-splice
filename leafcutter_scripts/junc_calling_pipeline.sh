#!/usr/bin/env bash
#SBATCH --job-name=junc_annot
#SBATCH --partition=bigmem
#SBATCH --mail-type=NONE
#SBATCH --mem=200g
#SBATCH --output=junc_annot_stdout_%j.log

##Pipeline to go from ONT Bam file (from SiCeLoRe output) --> junction:cell annotated matrix

##Format: bash junc_calling_pipeline.sh path/to/main_output_folder patient_run_name path/to/run_files
##Example: sbatch junc_calling_pipeline.sh /gpfs/commons/groups/landau_lab/SF3B1_splice_project/6.Junction_analysis_files/2.ONT_juncs CH259 run_files 

##Input varibles:
output_location=$1
input_bam=$2
patient_run=$3
run_files=$4

#mkdir "$output_location"/"$patient_run"_output

##### ----------------------- Generate Intron/3p/5p databases from GTF file (only need to do this once)  ----------------------- #####

#mkdir "$run_files"/annotation_reference

#Format: path/to/leafcutter/leafviz/gtf2leafcutter.pl -o path/to/output_files/prefix path/to/reference_gtf

#/gpfs/commons/home/gmullokandov/software/leafcutter/leafviz/gtf2leafcutter.pl -o "$run_files"/annotation_reference/leafviz "$run_files"/gencode.v31.basic.annotation.gtf

annotation_code=""$run_files"/annotation_reference/leafviz"


##### ----------------------- Run the python junction calling script ----------------------- #####

module load python/3.5.1

#Format: python script_to_run path/to/input/bam path_to_output_file
echo "starting junction calling" 

cd $output_location
mkdir leafcutter_outputs
cd ./leafcutter_outputs

python "$run_files"/count_introns_ONT.py $input_bam "$patient_run"_counts_sc_txt.gz
echo "junction calling done" 

##### ----------------------- Run the pre-processing scripts ----------------------- #####

echo "Running pre-processing scripts"

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/pchamely/.conda/envs/RenvForKnowlesPCA
mkdir "$patient_run"_output

#Format: Rscript junc_calling_script.R path/to/*_counts_sc_txt.gz path/to/output_folder patient_ID path/to/annotation_code/prefix"

Rscript "$run_files"/junc_calling_script.R "$patient_run"_counts_sc_txt.gz "$patient_run"_output "$patient_run" "$annotation_code"

echo "Done" 

##### ----------------------- Run Lloyd's annotation script  ----------------------- #####
conda deactivate
scp "$run_files"/annotation_reference/leafviz_all_introns.bed.gz "$patient_run"_output/
gunzip "$output_location"/leafcutter_outputs/"$patient_run"_output/leafviz_all_introns.bed.gz

#Format:python splice_annotator_got.py <path to datafiles> <input datafile> <junction tag .bed file> <output main file name>
python "$run_files"/new_annotator_v2.py "$output_location"/leafcutter_outputs/"$patient_run"_output "$patient_run"_all.introns.info.txt "$output_location"/leafcutter_outputs/"$patient_run"_output/leafviz_all_introns.bed "$patient_run"_all.introns.info.w.primary.annotations.txt                                                                                                                                                                                                          
##### ----------------------- Strand adjust the primary annotations  ----------------------- #####

 
