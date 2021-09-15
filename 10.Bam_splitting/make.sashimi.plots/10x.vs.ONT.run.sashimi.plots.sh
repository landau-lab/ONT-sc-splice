module load python
module load R/3.6.1
module load samtools

#ERGIC3
#./sashimi-plot.py -b input.bams/CH305.10x.vs.ONT.bams.tsv -c chr20:35548535-35557544 -M 5000 -C 3 -O 3 --shrink --alpha 0.25 -F pdf --base-size=20 --height=3 --width=18 -o 10x.vs.ONT/CH305.10x.vs.ONT.ERGIC3.sashimi


#Expanded view of gene:
#./sashimi-plot.py -b input.bams/MDSP1.ONT.full.bam.tsv -c chr20:35556000-35557600 --shrink -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 --alpha 1 -F pdf --base-size=20 --height=3 --width=18 -P palette.single.coralblue.txt -o 10x.vs.ONT/MDSP1.ONT.ERGIC3.sashimi

#./sashimi-plot.py -b input.bams/MDSP1.10x.full.bam.tsv -c 20:35556000-35557600 --shrink -M 5000 -C 3 -O 3 --alpha 1 -F pdf --base-size=20 --height=3 --width=18 -P palette.single.orange.txt -o 10x.vs.ONT/MDSP1.10x.ERGIC3.sashimi


#Zoom in on cryptic 3p splice site 
./sashimi-plot.py -b input.bams/MDSP1.ONT.full.bam.tsv -c chr20:35556904-35557004 --shrink -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 --alpha 1 -F pdf --base-size=20 --height=3 --width=12 -P palette.single.coralblue.txt -o 10x.vs.ONT/MDSP1.ONT.ERGIC3.cryptic.site.sashimi

./sashimi-plot.py -b input.bams/MDSP1.10x.full.bam.tsv -c 20:35556904-35557004 --shrink -M 5000 -C 3 -O 3 --alpha 1 -F pdf --base-size=20 --height=3 --width=12 -P palette.single.orange.txt -o 10x.vs.ONT/MDSP1.10x.ERGIC3.cryptic.site.sashimi
