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
 #   6177  TPOX # with at least 8 repeats: 1940, 640 of which are 11 repeats
 #  17319  vWA
 # 104050  Y-GATA-H4

# get fasta of trimmed reads annotated by predicted locus
awk -F"\t" '{print ">" $6 "_" $1 "\n" $2}' NISTB1.reorient.qualfilt.predictions | sed 's/ //g' > NISTB1.trimmed.fa
sed 's/_/ /' NISTB1.trimmed.fa > temp
mv temp NISTB1.trimmed.fa

# split reads by predicted locus
cat query.ids | while read line; do grep -w $line NISTB1.trimmed.fa -A 1 > split_files/$line.fa; done
cd split_files
rename 's/>//' *fa

# modify seq names so that they are uniqu
ls *fa | while read line; do sed -i s'/ /:/' $line; done

# dereplicate sequences
ls *fa | while read line; do sed -i 's/--//' $line; done
ls *fa | sed 's/.fa//' | parallel -j10 'vsearch --derep_fulllength {}.fa --output {}.uniq.fa --sizeout --minseqlength 25'

# sort by length
ls *uniq.fa | sed 's/.uniq.fa//' | while read line; do vsearch --sortbylength $line.uniq.fa --output $line.sort.fa; done

# cluster at 98% identity
ls *sort.fa | sed 's/.sort.fa//' | while read line; do vsearch --cluster_size $line.sort.fa --id 0.98 --sizein --sizeout --centroids $line.98.fa; done

# alignment
ls *98.fa | sed 's/.98.fa//' | parallel 'mafft {}.98.fa > {}.align.fa' 

# trim sequences and remove spurious alignments
# right now ignoring amelogenin and y chromosome cause slow
cat temp.ids | parallel -j10 'trimal -in {}.align.fa -out {}.trim.fa -gt 0.5 -st 0.001 -resoverlap 0.75 -seqoverlap 80'&



## STOPPED HERE # continued without those three jerk files that won't complete

# remove wordwrap
ls *trim* | while read line; do bash remove_wordwrap.sh $line > $line.fix; done
rm *trim.fa
rename  's/.fix//' *fix

# remove gaps
ls *trim* | sed 's/.trim.fa//' | while read line; do sed 's/-//g' $line.trim.fa > $line.nogap.fa; done

# sort by size
ls *nogap.fa | sed 's/.nogap.fa//' | while read line; do vsearch --sortbylength $line.nogap.fa --output $line.sorta.fa; done

# re cluster at 98%
ls *sorta.fa | sed 's/.sorta.fa//' | while read line; do vsearch --cluster_size $line.sorta.fa --id 0.98 --sizein --sizeout --centroids $line.98a.fa --minseqlength 25; done

# derep again
ls *98a.fa | sed 's/.98a.fa//' | while read line; do vsearch --derep_fulllength $line.98a.fa --output $line.uniq2.fa --sizeout --sizein --minseqlength 25; done

# sort by abundance
ls *uniq2.fa | sed 's/.uniq2.fa//' | while read line; do vsearch --sortbysize $line.uniq2.fa --output $line.final.fa; done

# clean up files
ls *final* | sed 's/.final.fa//' | while read line; do bash remove_wordwrap.sh $line.final.fa > $line.final.fix; done
rm *final.fa
rename 's/.final.fix/.final.fa/' *fix

# re add uniq sequence identifiers
ls *final* | while read line; do bash seq_numbers.sh $line > $line.fix; done
rm *final.fa
rename 's/.final.fix/.final.fa/' *fix

# find motifs
ls *final* | sed 's/.final.fa//' | while read line; do python find_motifs.py -f $line.final.fa > motifs/$line.motifs.txt; done




