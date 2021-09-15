
module load python
module load R/3.6.1
module load samtools

#----------- ERGIC3 -----------#
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr20:35556200-35557115 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.ERGIC3.sashimi

#----------- SERBP1 -----------#

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr1:67424880-67425205 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.SERBP1.junc2.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr1:67424880-67425205 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.SERBP1.junc2.sashimi

#----------- RPL22L1 -----------#

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.per.ct.tsv -c chr3:170868000-170868400 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.RPL22L1.sashimi.2

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr3:170868000-170868400 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.RPL22L1.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr3:170868000-170868400 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.RPL22L1.sashimi


#----------- OXA1L -----------#

#1UMI
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.1UMI.tsv -c chr14:22767950-22769900 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.1UMI.OXAIL.sashimi

#2UMI
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr14:22767950-22769900 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.OXAIL.sashimi

#----------- DPH5 -----------#
#Try doing 100 basepairs before and after the cryptic site
./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr1:100992704-100992804 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.DPH5.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr1:100992700-100992940 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.DPH5.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr1:100994853-100995153 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.DPH5.2.sashimi

#----------- METTL17 -----------#
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr14:20990470-20991250 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.METTL17.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr14:20991450-20992221 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.METTL17.2.sashimi

#----------- METTL17 -----------#
#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr12:50780200-50780500 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.ATF1.1.sashimi

#./sashimi-plot.py -b input.bams/MDS.untreated.combined.input.bams.2UMI.tsv -c chr12:50795773-50795999 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.triple.txt -o 2UMI.output/MDS.untreated.combined.bulk.mut.wt.ATF1.2.sashimi
