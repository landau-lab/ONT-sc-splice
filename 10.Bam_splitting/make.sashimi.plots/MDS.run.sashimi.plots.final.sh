
module load python
module load R/3.6.1
module load samtools


#----------- DPH5 -----------#
#Try doing 100 basepairs before and after the cryptic site
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr1:100992704-100992804 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/final.plots/MDS.untreated.combined.bulk.mut.wt.DPH5.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.allonly.2UMI.tsv -c chr1:100992700-100992940 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/final.plots/MDS.untreated.combined.bulk.mut.wt.DPH5.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.allonly.2UMI.tsv -c chr1:100994853-100995153 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/final.plots/MDS.untreated.combined.bulk.mut.wt.DPH5.2.sashimi

#----------- METTL17 -----------#
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr14:20990470-20991250 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.METTL17.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr14:20991450-20992221 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.METTL17.2.sashimi

#----------- ATF1 -----------#
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr12:50780200-50780500 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.ATF1.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr12:50795773-50795999 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.ATF1.2.sashimi

#----------- STAU1 -----------#

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.per.ct.2UMI.tsv -c chr20:49124460-49124796 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/final.plots/MDS.untreated.combined.bulk.mut.wt.STAU1.sashimi

#----------- GYPA -----------#
./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.bulk.and.ct.2UMI.tsv -c chr20:49124460-49124796 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/final.plots/MDS.untreated.combined.bulk.mut.wt.GYPA.sashimi
