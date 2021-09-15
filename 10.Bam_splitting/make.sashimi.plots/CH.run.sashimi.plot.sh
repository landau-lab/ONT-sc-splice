
module load python
module load R/3.6.1
module load samtools

#ERGIC3
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr20:35556200-35557115 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.ERGIC3.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr20:35556200-35557115 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.ERGIC3.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.bulk.input.bam.tsv -c chr20:35556200-35557115 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.single.txt -o CH.combined.bulk.ERGIC3.sashimi

#Try doing 50 bp before and after cryptic junction and include gencode reference 
./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr20:35556904-35557004 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F jpeg -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.ERGIC3.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr20:35556200-35557115 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=22 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.ERGIC3.sashimi.full

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr20:35556200-35557115 --shrink -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.ERGIC3.sashimi.full.w.shrink

#SERBP1
#./sashimi-plot.py -b input.bams/CH259.input.bams.per.cttsv -c chr1:67424111-67424959 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.txt -o SERBP1_sashimi
#./sashimi-plot.py -b input.bams/CH259.input.bams.per.ct.tsv -c chr1:67424111-67424959 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o SERBP1_sashimi_2

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.per.ct.tsv -c chr1:67424111-67424959 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.multi.txt -o CH.combined.SERBP1.sashimi.1
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.per.ct.tsv -c chr1:67424111-67424959 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.SERBP1.sashimi.2

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr1:67424111-67424959 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.SERBP1.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr1:67424111-67424959 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.SERBP1.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr1:67424880-67425205 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.SERBP1.junc2.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr1:67424880-67425205 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.SERBP1.junc2.sashimi

#RPL22L1
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.per.ct.tsv -c chr3:170868000-170868400 -M 1000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.RPL22L1.sashimi.2

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.tsv -c chr3:170868000-170868400 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.RPL22L1.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr3:170868000-170868400 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.RPL22L1.sashimi

#OXA1L
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.1UMI.tsv -c chr14:22767950-22769900 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o CH.combined.all.mut.wt.1UMI.OXAIL.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.2UMI.tsv -c chr14:22767950-22769900 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.double.txt -o 2UMI.output/CH.combined.all.mut.wt.2UMI.OXAIL.sashimi

#./sashimi-plot.py -b input.bams/CH.combined.bulk.input.bam.tsv -c chr14:22767950-22769900 -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=16 --height=3 --width=18 -P palette.single.txt -o CH.combined.bulk.OXAIL.sashimi

#SLC25A3
#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr17:44322502-44322602 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.SLC25A3.sashimi


#DPH5

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr1:100992704-100992804 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.DPH5.sashimi

#METTL5

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr2:169812474-169812574 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.METTL5.sashimi

#METTL5

#./sashimi-plot.py -b input.bams/CH.combined.input.bams.bulk.and.ct.1UMI.tsv -c chr1:161167288-161167388 -g gencode.v31.basic.annotation.gtf -M 5000 -C 3 -O 3 -A median --alpha 1 -F pdf -R 350 --base-size=24 --height=3 --width=18 -P palette.double.txt -o 1UMI.output/CH.combined.bulk.and.ct.1UMI.PPOX.sashimi
