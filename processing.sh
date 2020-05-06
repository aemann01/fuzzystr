#########################
#STR FUZZY MATCH PROJECT
#########################
cd /Users/mann/github/NISTB

# got bam files from courtney
cd /home/administrator/NISTB1_alignment
bedtools bamtobed -i NISTB1.sorted.bam > NISTB1.sorted.bed

# how often do sequences have more than one alignment?
awk -F"\t" '{print $4}' NISTB1.sorted.bed | sort | uniq -c | awk '{print $1}' | sort | uniq -c

# first need to convert bam to sam then to paf
samtools view -h NISTB1.sorted.bam > NISTB1.sorted.sam
paftools.js sam2paf NISTB1.sorted.sam > NISTB1.sorted.paf

# first train the model
reorientexpress.py -train -data NISTB1_ps109_75ng_051419.fastq.gz -source mapped --a NISTB1.sorted.paf --v -output NISTB1
# test the model
reorientexpress.py -test -dat NISTB1_ps109_75ng_051419.fastq.gz -source mapped -model NISTB1.model 1> model_test.out 2> model_test.err
# predictions
reorientexpress.py -predict -data NISTB1_ps109_75ng_051419.fastq.gz -source mapped -model NISTB1.model -output NISTB1.predictions
# convert predictions to fasta
awk -F"," '{print ">"$1, "\n", $2}' NISTB1.predictions.csv | sed 's/^ //' > NISTB1.reorient.fa

# filter out low quality mapping reads
samtools view -h -q 21 NISTB1.sorted.bam > NISTB1.qualfilt.bam
bedtools bamtobed -i NISTB1.qualfilt.bam > NISTB1.qualfilt.bed

# filter fasta file - include only quality filter threshold passes
awk -F"\t" '{print $4}' NISTB1.qualfilt.bed | sort | uniq > qualfilt.ids
seqtk subseq NISTB1.reorient.fa qualfilt.ids > NISTB1.reorient.qualfilt.fa

# split fasta file into smaller chunks
python3 fasta_chunk.py NISTB1.reorient.qualfilt.fa 10000

# now can run fuzzy str to get predicted locus information and trimmed sequences
ls group_*fasta | sed 's/.fasta//' | parallel 'python3 fuzzystr.py -f {}.fasta -p 150-300bp_amps_powerseq/Powerseq.config -e 0.20 > {}.out'

# cat together get rid of intermediate files
cat *out > NISTB1.reorient.qualfilt.predictions
rm *out
rm group*

# how many reads with just one hit vs multiple hits?
awk '{print $1}' NISTB1.reorient.qualfilt.predictions | sort | uniq -c | grep "1 " > hit_once.txt
awk '{print $1}' NISTB1.reorient.qualfilt.predictions | sort | uniq -c | grep -v "1 " > hit_morethanonce.txt

# how many reads per locus?
awk -F"\t" '{print $6}' NISTB1.reorient.qualfilt.predictions | sort | uniq -c
 # 170987  Amelogenin
 #   4732  CSF1PO
 #   2908  D10S1248
 #   4334  D12S391
 #  13221  D13S317
 #  13773  D16S539
 #   1757  D18S51
 #  32370  D19S433
 #   4021  D1S1656
 #   3169  D21S11
 #   1549  D22S1045
 #  21147  D2S1338
 #   8396  D2S441
 #   8017  D3S1358
 #   8560  D5S818
 #  15016  D7S820
 #  46468  D8S1179
 #  16045  DYS19
 #   3118  DYS385
 #   9873  DYS389I
 #   9873  DYS389II
 #  17716  DYS390
 #  13829  DYS391
 #   2144  DYS392
 #   5290  DYS393
 #   3726  DYS437
 #   2008  DYS438
 #  39283  DYS439
 #   9998  DYS448
 #   8805  DYS456
 #   3757  DYS458
 #   9158  DYS481
 #  14259  DYS533
 #  25957  DYS549
 #  12479  DYS570
 #   5367  DYS576
 #   8594  DYS635
 #   4937  DYS643
 #   2692  FGA
 #    402  PENTAD
 #  16656  PENTAE
 #   1034  TH01
 #   6177  TPOX # with at least 8 repeats: 1300, 640 11 repeats
 #  17319  vWA
 # 104050  Y-GATA-H4

# get fasta of trimmed reads annotated by predicted locus
awk -F"\t" '{print ">" $6 "_" $1 "\n" $2}' NISTB1.reorient.qualfilt.predictions | sed 's/ //g' > NISTB1.trimmed.fa
sed 's/_/ /' NISTB1.trimmed.fa > temp
mv temp NISTB1.trimmed.fa

# split reads by predicted locus
mkdir split_files
cat query.ids | while read line; do grep -w $line NISTB1.trimmed.fa -A 1 > split_files/$line.fa; done
cd split_files
rename 's/>//' *fa

# modify seq names so that they are uniqu
ls *fa | while read line; do sed -i s'/ /:/' $line; done

# now trim by flanking regions
# only doing those loci that are in the supp table
# got to be a better way to do this, iterate through files etc
mkdir flanktrim
python ../scripts/flanktrim.py -f D1S1656.fa -fw GAACCAAATA \
	-rv GCAACACAGG -e 0.1 > flanktrim/D1S1656.trim &
python ../scripts/flanktrim.py -f D13S317.fa -fw ACAAATACAT \
	-rv ATCATCTAT -e 0.1 > flanktrim/D13S317.trim &
python ../scripts/flanktrim.py -f PENTAE.fa -fw CTTACAATTT \
	-rv GAGACTGAGTCT -e 0.1 > flanktrim/PENTAE.trim &
python ../scripts/flanktrim.py -f TPOX.fa -fw CACTGAATGA \
	-rv AATGTTTGG -e 0.1 > flanktrim/TPOX.trim &
python ../scripts/flanktrim.py -f D2S441.fa -fw TATGAAAACT \
	-rv TATCATAACA -e 0.1 > flanktrim/D2S441.trim &
python ../scripts/flanktrim.py -f D2S1338.fa -fw GAAGGAAGGA \
	-rv AGGCCAAGCC -e 0.1 > flanktrim/D2S1338.trim &
python ../scripts/flanktrim.py -f D3S1358.fa -fw CTTGCATGTA \
	-rv TGAGACAGGG -e 0.1 > flanktrim/D3S1358.trim &
python ../scripts/flanktrim.py -f FGA.fa -fw AAGGAAGAAA \
	-rv CTAGCTTGTA -e 0.1 > flanktrim/FGA.trim &
python ../scripts/flanktrim.py -f D5S818.fa -fw TTATACCTCT \
	-rv TCAAAAT -e 0.1 > flanktrim/D5S818.trim &
python ../scripts/flanktrim.py -f CSF1PO.fa -fw ATCTATCTAT \
	-rv AATCTATCTA -e 0.1 > flanktrim/CSF1PO.trim &
python ../scripts/flanktrim.py -f D7S820.fa -fw TCAATCTGTC \
	-rv GTTAGTTCGT -e 0.1 > flanktrim/D7S820.trim &
python ../scripts/flanktrim.py -f D10S1248.fa -fw ATGAGTGAGT \
	-rv ATGAAGACAA -e 0.1 > flanktrim/D10S1248.trim &
python ../scripts/flanktrim.py -f TH01.fa -fw CTCCATGGTG \
	-rv AGGGAAATAA -e 0.1 > flanktrim/TH01.trim &
python ../scripts/flanktrim.py -f vWA.fa -fw GATAGATGGA \
	-rv TCAAT -e 0.1 > flanktrim/vWA.trim &
python ../scripts/flanktrim.py -f D12S391.fa -fw ATGCATAGGT \
	-rv GAGAGGGGAT -e 0.1 > flanktrim/D12S391.trim &
python ../scripts/flanktrim.py -f D16S539.fa -fw CAGACAGGTG \
	-rv TCATTGAAAG -e 0.1 > flanktrim/D16S539.trim &
python ../scripts/flanktrim.py -f D18S51.fa -fw GTCTCAGAAA \
	-rv AAAGAGAGAG -e 0.1 > flanktrim/D18S51.trim &
python ../scripts/flanktrim.py -f D19S433.fa -fw CTTCCTCTCT \
	-rv CAACAGAATC -e 0.1 > flanktrim/D19S433.trim &
python ../scripts/flanktrim.py -f D21S11.fa -fw TGAATTGCCT \
	-rv TCGTCTATCT -e 0.1 > flanktrim/D21S11.trim &
python ../scripts/flanktrim.py -f PENTAD.fa -fw ACCATCTCTC \
	-rv AAAAACGAA -e 0.1 > flanktrim/PENTAD.trim &
python ../scripts/flanktrim.py -f D22S1045.fa -fw TAGTAGTCTC \
	-rv GTTATAAAAA -e 0.25 > flanktrim/D22S1045.trim &
python ../scripts/flanktrim.py -f D8S1179.fa -fw GATCTATCTA \
	-rv TTCCC -e 0.1 > flanktrim/D8S1179.trim &

# get fasta for each trimmed file
cd flanktrim
ls *trim | while read line; do awk -F"\t" '{print ">" $1 "\n" $2}' $line > $line.fa; done

# sort by length
ls *.fa | sed 's/.trim.fa//' | while read line; do vsearch --sortbylength $line.trim.fa --output $line.sort.fa; done

# dereplication
ls *sort.fa | sed 's/.sort.fa//' | while read line; do vsearch --derep_fulllength $line.sort.fa --output $line.derep.fa --sizeout --minseqlength 10; done

# cluster at 99% identity
ls *derep.fa | sed 's/.derep.fa//' | parallel 'vsearch --cluster_smallmem {}.derep.fa --id 0.99 --centroids {}.cluster.fa --clusterout_sort --consout {}.cons.fa --relabel_keep --sizein --sizeout --uc {}.uc --usersort'

# map to reference sequences
# bowtie2-build ref_seqs.fa ref_seqs
ls *cons* | sed 's/.fa//' | while read line; do bowtie2 --no-unal -x ~/mann/NISTB1_alignment/ref_seqs -U $line.fa -f --end-to-end -S $line.sam; done

# also try to use blast to find best guessed match
cat *cons.fa > cleaned.str.fa
blastn -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send evalue length pident gaps btop' \
	-db ~/mann/NISTB1_alignment/ref_seqs.db \
	-query cleaned.str.fa \
	-out blast.out \
	-evalue 1e-10 \
	-perc_identity 0.90 \
	-gapopen 1 \
	-gapextend 1 \
	-max_target_seqs 1 &

# also get xml version
blastn -outfmt 5 \
	-db ~/mann/NISTB1_alignment/ref_seqs.db \
	-query cleaned.str.fa \
	-out blast.out.xml \
	-evalue 1e-10 \
	-perc_identity 0.90 \
	-gapopen 1 \
	-gapextend 1 \
	-max_target_seqs 1 &

# ok so blast is not great at loci like tpox, need alternative

# get predicted motif counts
# each file visually examined to determine number of max allowed mismatches
# check for each motif IF it repeats, can narrow down later
mkdir motifs
## CSF1P0
Rscript ../../scripts/motif_count.R CSF1PO.cons.fa "ATCT" 0 CSF1P0.motifs
## D10S1248
Rscript ../../scripts/motif_count.R D10S1248.cons.fa "GGAA" 0 D10S1248.motifs
## D12S391
Rscript ../../scripts/motif_count.R D12S391.cons.fa "AGAT" 0 D12S391.motifs_1
Rscript ../../scripts/motif_count.R D12S391.cons.fa "AGAC" 0 D12S391.motifs_2
paste D12S391.motifs_1 D12S391.motifs_2 > D12S391.motifs
## D13S317
Rscript ../../scripts/motif_count.R D13S317.cons.fa "TATC" 0 D13S317.motifs
## D16S539
Rscript ../../scripts/motif_count.R D16S539.cons.fa "GATA" 0 D16S539.motifs
## D18S51
Rscript ../../scripts/motif_count.R D18S51.cons.fa "AGAA" 0 D18S51.motifs
## D19S433
Rscript ../../scripts/motif_count.R D19S433.cons.fa "CCTT" 1 D19S433.motifs
## D1S1656
Rscript ../../scripts/motif_count.R D1S1656.cons.fa "TCTA" 0 D1S1656.motifs
## D21S11
Rscript ../../scripts/motif_count.R D21S11.cons.fa "TCTA" 0 D21S11.motifs_1
Rscript ../../scripts/motif_count.R D21S11.cons.fa "TCTG" 0 D21S11.motifs_2
paste D21S11.motifs_1 D21S11.motifs_2 > D21S11.motifs
## D22S1045
Rscript ../../scripts/motif_count.R D22S1045.cons.fa "ATT" 0 D22S1045.motifs
## D2S1338
Rscript ../../scripts/motif_count.R D2S1338.cons.fa "GGAA" 0 D2S1338.motifs_1
Rscript ../../scripts/motif_count.R D2S1338.cons.fa "GGCA" 0 D2S1338.motifs_2
paste D2S1338.motifs_1 D2S1338.motifs_2 > D2S1338.motifs
## D2S441
Rscript ../../scripts/motif_count.R D2S441.cons.fa "TCTA" 0 D2S441.motifs
## D3S1358
Rscript ../../scripts/motif_count.R D3S1358.cons.fa "TCTG" 0 D3S1358.motifs_1
Rscript ../../scripts/motif_count.R D3S1358.cons.fa "TCTA" 0 D3S1358.motifs_2
paste D3S1358.motifs_1 D3S1358.motifs_2 > D3S1358.motifs
## D5S818
Rscript ../../scripts/motif_count.R D5S818.cons.fa "ATCT" 0 D5S818.motifs
## D7S820
Rscript ../../scripts/motif_count.R D7S820.cons.fa "TATC" 0 D7S820.motifs
## D8S1179
Rscript ../../scripts/motif_count.R D8S1179.cons.fa "TCTA" 0 D8S1179.motifs_1
Rscript ../../scripts/motif_count.R D8S1179.cons.fa "TCTG" 0 D8S1179.motifs_2
paste D8S1179.motifs_1 D8S1179.motifs_2 > D8S1179.motifs
## FGA
Rscript ../../scripts/motif_count.R FGA.cons.fa "GGAA" 0 FGA.motifs_1
Rscript ../../scripts/motif_count.R FGA.cons.fa "AAAG" 0 FGA.motifs_2
Rscript ../../scripts/motif_count.R FGA.cons.fa "GAAA" 0 FGA.motifs_3
Rscript ../../scripts/motif_count.R FGA.cons.fa "GAAG" 0 FGA.motifs_4
paste FGA.motifs_1 FGA.motifs_2 FGA.motifs_3 FGA.motifs_4 > FGA.motifs
## PENTAD
Rscript ../../scripts/motif_count.R PENTAD.cons.fa "AAAGA" 0 PENTAD.motifs
## PENTAE
Rscript ../../scripts/motif_count.R PENTAE.cons.fa "TCTTT" 0 PENTAE.motifs
## TH01
Rscript ../../scripts/motif_count.R TH01.cons.fa "AATG" 0 TH01.motifs
## TPOX
Rscript ../../scripts/motif_count.R TPOX.cons.fa "AATG" 0 TPOX.motifs
## vWA
Rscript ../../scripts/motif_count.R vWA.cons.fa "TAGA" 0 vWA.motifs_1
Rscript ../../scripts/motif_count.R vWA.cons.fa "CAGA" 0 vWA.motifs_2
paste vWA.motifs_1 vWA.motifs_2 >vWA.motifs

## RERUN FROM FLANK TRIM WHEN DONE




