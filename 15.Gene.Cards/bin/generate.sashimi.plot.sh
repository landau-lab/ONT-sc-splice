module load python
module load R/3.6.1
module load samtools

input_bams=$1
cryptic_junction=$2
output_file=$3

working_directory="/gpfs/commons/groups/landau_lab/SF3B1_splice_project/15.Gene.Cards/bin"
python "$working_directory"/sashimi-plot.py -b "$input_bams" -c "$cryptic_junction" -g "$working_directory"/gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F png -R 350 --base-size=24 --height=3 --width=18 -P "$working_directory"/palette.double.txt -o "$output_file"

