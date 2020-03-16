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

# now can run fuzzy str to get predicted locus information and trimmed sequences
python3 fuzzystr.py -f NISTB1.reorient.qualfilt.fa -o NISTB1.reorient.qualfilt.predict -p 150-300bp_amps_powerseq/Powerseq.config -e 0.25

# split reads by predicted locus

# map each batch to the reference dataset

# get consensus for each batch

# motif count and haplotype calling from consensus sequence?


