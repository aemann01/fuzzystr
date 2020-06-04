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

# rereplicate for motif counts
ls *cons.fa | sed 's/.cons.fa//' | while read line; do vsearch --rereplicate $line.cons.fa --output $line.rerep.fa; done

# get predicted motif counts
# each file visually examined to determine number of max allowed mismatches
# check for each motif IF it repeats, can narrow down later
mkdir motifs
## CSF1P0
Rscript ../../scripts/motif_count.R CSF1PO.rerep.fa "ATCT" 0 CSF1P0.motifs
## D10S1248
Rscript ../../scripts/motif_count.R D10S1248.rerep.fa "GGAA" 0 D10S1248.motifs
## D12S391
Rscript ../../scripts/motif_count.R D12S391.rerep.fa "AGAT" 0 D12S391.motifs_1
Rscript ../../scripts/motif_count.R D12S391.rerep.fa "AGAC" 0 D12S391.motifs_2
paste D12S391.motifs_1 D12S391.motifs_2 > D12S391.motifs
## D13S317
Rscript ../../scripts/motif_count.R D13S317.rerep.fa "TATC" 0 D13S317.motifs
## D16S539
Rscript ../../scripts/motif_count.R D16S539.rerep.fa "GATA" 0 D16S539.motifs
## D18S51
Rscript ../../scripts/motif_count.R D18S51.rerep.fa "AGAA" 0 D18S51.motifs
## D19S433
Rscript ../../scripts/motif_count.R D19S433.rerep.fa "CCTT" 1 D19S433.motifs
## D1S1656
Rscript ../../scripts/motif_count.R D1S1656.rerep.fa "TCTA" 0 D1S1656.motifs
## D21S11
Rscript ../../scripts/motif_count.R D21S11.rerep.fa "TCTA" 0 D21S11.motifs_1
Rscript ../../scripts/motif_count.R D21S11.rerep.fa "TCTG" 0 D21S11.motifs_2
paste D21S11.motifs_1 D21S11.motifs_2 > D21S11.motifs
## D22S1045
Rscript ../../scripts/motif_count.R D22S1045.rerep.fa "ATT" 0 D22S1045.motifs
## D2S1338
Rscript ../../scripts/motif_count.R D2S1338.rerep.fa "GGAA" 0 D2S1338.motifs_1
Rscript ../../scripts/motif_count.R D2S1338.rerep.fa "GGCA" 0 D2S1338.motifs_2
paste D2S1338.motifs_1 D2S1338.motifs_2 > D2S1338.motifs
## D2S441
Rscript ../../scripts/motif_count.R D2S441.rerep.fa "TCTA" 0 D2S441.motifs
## D3S1358
Rscript ../../scripts/motif_count.R D3S1358.rerep.fa "TCTG" 0 D3S1358.motifs_1
Rscript ../../scripts/motif_count.R D3S1358.rerep.fa "TCTA" 0 D3S1358.motifs_2
paste D3S1358.motifs_1 D3S1358.motifs_2 > D3S1358.motifs
## D5S818
Rscript ../../scripts/motif_count.R D5S818.rerep.fa "ATCT" 0 D5S818.motifs
## D7S820
Rscript ../../scripts/motif_count.R D7S820.rerep.fa "TATC" 0 D7S820.motifs
## D8S1179
Rscript ../../scripts/motif_count.R D8S1179.rerep.fa "TCTA" 0 D8S1179.motifs_1
Rscript ../../scripts/motif_count.R D8S1179.rerep.fa "TCTG" 0 D8S1179.motifs_2
paste D8S1179.motifs_1 D8S1179.motifs_2 > D8S1179.motifs
## FGA
Rscript ../../scripts/motif_count.R FGA.rerep.fa "GGAA" 0 FGA.motifs_1
Rscript ../../scripts/motif_count.R FGA.rerep.fa "AAAG" 0 FGA.motifs_2
Rscript ../../scripts/motif_count.R FGA.rerep.fa "GAAA" 0 FGA.motifs_3
Rscript ../../scripts/motif_count.R FGA.rerep.fa "GAAG" 0 FGA.motifs_4
paste FGA.motifs_1 FGA.motifs_2 FGA.motifs_3 FGA.motifs_4 > FGA.motifs
## PENTAD
Rscript ../../scripts/motif_count.R PENTAD.rerep.fa "AAAGA" 0 PENTAD.motifs
## PENTAE
Rscript ../../scripts/motif_count.R PENTAE.rerep.fa "TCTTT" 0 PENTAE.motifs
## TH01
Rscript ../../scripts/motif_count.R TH01.rerep.fa "AATG" 0 TH01.motifs
## TPOX
Rscript ../../scripts/motif_count.R TPOX.rerep.fa "AATG" 0 TPOX.motifs
## vWA
Rscript ../../scripts/motif_count.R vWA.rerep.fa "TAGA" 0 vWA.motifs_1
Rscript ../../scripts/motif_count.R vWA.rerep.fa "CAGA" 0 vWA.motifs_2
paste vWA.motifs_1 vWA.motifs_2 >vWA.motifs
# clean up
rm *_*

#######################################
# New test, STRique
#######################################
# activate environment
cd /home/administrator/STRique
source bin/activate
# 1 indexing
python3 ~/src/STRique/scripts/STRique.py index \
	~/mann/15cycle_NISTB/ \
	> ~/mann/NISTB1_alignment/strique/reads.fofn

# bam to sam
# samtools view -h -o ~/mann/NISTB1_alignment/NISTB1_ps109_75ng_051419.sam \
# 	~/mann/NISTB1_alignment/NISTB1_ps109_75ng_051419.bam

# filter out unmapped reads
samtools view -F 4 ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.sam \
	> ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam

# repeat counting 150 bp flank
cat ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam | \
	python3 ~/src/STRique/scripts/STRique.py count \
	~/mann/NISTB1_alignment/strique/reads.fofn \
	~/src/STRique/models/r9_4_450bps.model \
	~/src/STRique/configs/strique_config_150bp.txt \
	> ~/mann/NISTB1_alignment/strique/results.tsv

# repeat counting 75 bp flank
cat ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam | \
	python3 ~/src/STRique/scripts/STRique.py count \
	~/mann/NISTB1_alignment/strique/reads.fofn \
	~/src/STRique/models/r9_4_450bps.model \
	~/src/STRique/configs/strique_config_75bp.txt \
	> ~/mann/NISTB1_alignment/strique/results_75bp.tsv

# repeat counting 50 bp flank
cat ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam | \
	python3 ~/src/STRique/scripts/STRique.py count \
	~/mann/NISTB1_alignment/strique/reads.fofn \
	~/src/STRique/models/r9_4_450bps.model \
	~/src/STRique/configs/strique_config_50bp.txt \
	> ~/mann/NISTB1_alignment/strique/results_50bp.tsv

# repeat counting 25 bp flank
cat ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam | \
	python3 ~/src/STRique/scripts/STRique.py count \
	~/mann/NISTB1_alignment/strique/reads.fofn \
	~/src/STRique/models/r9_4_450bps.model \
	~/src/STRique/configs/strique_config_25bp.txt \
	> ~/mann/NISTB1_alignment/strique/results_25bp.tsv

# first need to filter out any reads with repeat count of zero
ls results* | while read line; 
	do awk -F"\t" '$4>0' $line > $line.filt; done

# sanity check
awk -F"\t" '{print $4}' *filt | sort | uniq -c

# round up prefix and suffix score
ls *filt | while read line; 
	do awk -F"\t" '{printf "%.1f\t%.1f\n", $5, $6}' $line > $line.round; done
ls *round | while read line; 
	do sed -i '1 s/^.*$/score_prefix_round\tscore_suffix_round/' $line; done
ls *filt | sed 's/.filt//' | while read line; 
	do paste $line.filt $line.filt.round > $line.paste; done

# next filter by prefix and suffix score (being a little lax here, in pub 4.0 is threshold)
ls *paste | sed 's/.paste//' | while read line; 
	do awk -F"\t" '$11>=3.5 && $12>=3.5' $line.paste > $line.final; done

# pass to hist generation r file

# Why so many reads with 3 repeats?
# first pull read IDs
grep "TPOX" results_25bp.tsv.final | awk -F"\t" '$4==3' | awk '{print $1}' > tpox_3rep_25bp.ids
# actually not that many only 45 reads total?

# filter from fastq
seqtk subseq ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.fastq tpox_3rep_25bp.ids > tpox_3rep_25bp.fq

# fastq to fasta
sed -n '1~4s/^@/>/p;2~4p' tpox_3rep_25bp.fq > tpox_3rep_25bp.fa

# align sequences
mafft tpox_3rep_25bp.fa > tpox_3rep_25bp.align.fa

# only 145 reads were filtered from tpox (25bp)
# 45 are 3 reps
# 20 are 8 etc...

# what about before all of the filtering
# re run strique
cd /home/administrator/STRique
source bin/activate
cd ~/mann/NISTB1_alignment/strique/test/

cat ~/mann/15cycle_NISTB/15cycle_NISTB_hacbasecall/15cycle_NISTB.mapped.sam | \
	python3 ~/src/STRique/scripts/STRique.py count \
	~/mann/NISTB1_alignment/strique/reads.fofn \
	~/src/STRique/models/r9_4_450bps.model \
	~/src/STRique/configs/strique_config_25bp.txt \
	> ~/mann/NISTB1_alignment/strique/test/results_25bp.tsv

# trying to correct nanopore error
# first need to download locally, need to update ubuntu on lab computer to run
scp administrator@10.53.1.235:/home/administrator/mann/NISTB1_alignment/split_files/*fa .
cat *fa > all.fa
lorma.sh all.fa

# dereplicate these
vsearch --derep_fulllength final.fasta --output final.derep.fa --sizeout --minseqlen 100

# what if we bin these by length? Uncorrected ones lormah makes it worse
# read length statistics for each locus
ls *.sort.fa | sed 's/.sort.fa//' | while read line; do ~/bbmap/readlength.sh in=$line.sort.fa out=$line.sum.txt bin=4; done

###################
# LENGTH BINNING
###################
# first need to get rid of reads in the flank trimmed data that fall above or below the min max Courtney sent
ls *.derep.fa | sed 's/.derep.fa//' | while read line; do python ~/github/NISTB/scripts/fasta_len.py $line.derep.fa > $line.len; done

# get summary statistics for each (median +- 1stdev)
# pull reads +- 4bp around median
# CSF1P0 24 19 29
awk '$2 >= 19 && $2 <= 29' CSF1PO.len | awk '{print $1}' > CSF1PO.filt.ids
# D10S1248 52 47 57
awk '$2 >= 47 && $2 <= 57' D10S1248.len | awk '{print $1}' > D10S1248.filt.ids
# D12S391 92 79 105
awk '$2 >= 79 && $2 <= 105' D12S391.len | awk '{print $1}' > D12S391.filt.ids
# D13S317 44 35 53
awk '$2 >= 35 && $2 <= 53' D13S317.len | awk '{print $1}' > D13S317.filt.ids
# D16S539 48 41 55
awk '$2 >= 24 && $2 <= 32' D16S539.len | awk '{print $1}' > D16S539.filt.ids
# D18S51 48 39 57
awk '$2 >= 39 && $2 <= 57' D18S51.len | awk '{print $1}' > D18S51.filt.ids
# D19S433 76 66 86
awk '$2 >= 66 && $2 <= 86' D19S433.len | awk '{print $1}' > D19S433.filt.ids
# D1S1656 104 95 113
awk '$2 >= 95 && $2 <= 113' D1S1656.len | awk '{print $1}' > D1S1656.filt.ids
# D21S11 128 77 179
awk '$2 >= 77 && $2 <= 179' D21S11.len | awk '{print $1}' > D21S11.filt.ids
# D22S1045 48 44 53
awk '$2 >= 44 && $2 <= 53' D22S1045.len | awk '{print $1}' > D22S1045.filt.ids
# D2S1338 40 29 51
awk '$2 >= 29 && $2 <= 51' D2S1338.len | awk '{print $1}' > D2S1338.filt.ids
# D2S441 52 43 61
awk '$2 >= 43 && $2 <= 61' D2S441.len | awk '{print $1}' > D2S441.filt.ids
# D3S1358 72 63 81
awk '$2 >= 63 && $2 <= 81' D3S1358.len | awk '{print $1}' > D3S1358.filt.ids
# D5S818 48 38 58
awk '$2 >= 38 && $2 <= 58' D5S818.len | awk '{print $1}' > D5S818.filt.ids
# D7S820 40 37 43
awk '$2 >= 37 && $2 <= 43' D7S820.len | awk '{print $1}' > D7S820.filt.ids
# D8S1179 32 24 40
awk '$2 >= 24 && $2 <= 40' D8S1179.len | awk '{print $1}' > D8S1179.filt.ids
# FGA 72 60 85
awk '$2 >= 60 && $2 <= 85' FGA.len | awk '{print $1}' > FGA.filt.ids
# PENTAD 52 47 57
awk '$2 >= 47 && $2 <= 57' PENTAD.len | awk '{print $1}' > PENTAD.filt.ids
# PENTAE 72 51 93
awk '$2 >= 51 && $2 <= 93' PENTAE.len | awk '{print $1}' > PENTAE.filt.ids
# TH01 36 28 44
awk '$2 >= 28 && $2 <= 44' TH01.len | awk '{print $1}' > TH01.filt.ids
# TPOX 28 21 35
awk '$2 >= 21 && $2 <= 35' TPOX.len | awk '{print $1}' > TPOX.filt.ids
# vWA 68 52 84
awk '$2 >= 52 && $2 <= 84' vWA.len | awk '{print $1}' > vWA.filt.ids

# filter binned reads from sequences
ls *ids | sed 's/.filt.ids//' | while read line; do seqtk subseq $line.sort.fa $line.filt.ids > $line.bin.fa; done

# dereplicate again and sort by abundance
ls *bin.fa | sed 's/.bin.fa//' | while read line; do vsearch --derep_fulllength $line.bin.fa --output $line.bin.derep.fa --sizeout --minseqlength 10; done
ls *bin.derep.fa | sed 's/.bin.derep.fa//' | while read line; do vsearch --sortbysize $line.bin.derep.fa --output $line.bin.abundance.fa; done







