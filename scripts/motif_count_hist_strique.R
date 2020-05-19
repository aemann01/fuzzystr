#!/usr/bin/env Rscript
library(ggplot2)

# first run with 150bp flank regions
dat <- read.table("results.tsv.final", header=T, sep="\t")

# single repeats
## CSF1P0
dat.sub <- dat[dat$target == "CSF1PO",]
pdf("CSF1P0.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D10S1248
dat.sub <- dat[dat$target == "D10S1248",]
pdf("D10S1248.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D13S317
dat.sub <- dat[dat$target == "D13S317",]
pdf("D13S317.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D16S539
dat.sub <- dat[dat$target == "D16S539",]
pdf("D16S539.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D18S51
dat.sub <- dat[dat$target == "D18S51",]
pdf("D18S51.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D19S433
dat.sub <- dat[dat$target == "D19S433",]
pdf("D19S433.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D1S1656
dat.sub <- dat[dat$target == "D1S1656",]
pdf("D1S1656.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D22S1045
dat.sub <- dat[dat$target == "D22S1045",]
pdf("D22S1045.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D2S441
dat.sub <- dat[dat$target == "D2S441",]
pdf("D2S441.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D5S818
dat.sub <- dat[dat$target == "D5S818",]
pdf("D5S818.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D7S820
dat.sub <- dat[dat$target == "D7S820",]
pdf("D7S820.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAD
dat.sub <- dat[dat$target == "Penta_D",]
pdf("PENTAD.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAE
dat.sub <- dat[dat$target == "Penta_E",]
pdf("PENTAE.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TH01
dat.sub <- dat[dat$target == "TH01",]
pdf("TH01.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TPOX
dat.sub <- dat[dat$target == "TPOX",]
pdf("TPOX.hist.strique.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

##########
# 75bp
##########
dat <- read.table("results_75bp.tsv.final", header=T, sep="\t")
# single repeats
## CSF1P0
dat.sub <- dat[dat$target == "CSF1PO",]
pdf("CSF1P0.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D10S1248
dat.sub <- dat[dat$target == "D10S1248",]
pdf("D10S1248.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D13S317
dat.sub <- dat[dat$target == "D13S317",]
pdf("D13S317.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D16S539
dat.sub <- dat[dat$target == "D16S539",]
pdf("D16S539.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D18S51
dat.sub <- dat[dat$target == "D18S51",]
pdf("D18S51.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D19S433
dat.sub <- dat[dat$target == "D19S433",]
pdf("D19S433.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D1S1656
dat.sub <- dat[dat$target == "D1S1656",]
pdf("D1S1656.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D22S1045
dat.sub <- dat[dat$target == "D22S1045",]
pdf("D22S1045.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D2S441
dat.sub <- dat[dat$target == "D2S441",]
pdf("D2S441.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D5S818
dat.sub <- dat[dat$target == "D5S818",]
pdf("D5S818.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D7S820
dat.sub <- dat[dat$target == "D7S820",]
pdf("D7S820.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAD
dat.sub <- dat[dat$target == "Penta_D",]
pdf("PENTAD.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAE
dat.sub <- dat[dat$target == "Penta_E",]
pdf("PENTAE.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TH01
dat.sub <- dat[dat$target == "TH01",]
pdf("TH01.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TPOX
dat.sub <- dat[dat$target == "TPOX",]
pdf("TPOX.hist.strique75.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

##########
# 50bp
##########
dat <- read.table("results_50bp.tsv.final", header=T, sep="\t")
# single repeats
## CSF1P0
dat.sub <- dat[dat$target == "CSF1PO",]
pdf("CSF1P0.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D10S1248
dat.sub <- dat[dat$target == "D10S1248",]
pdf("D10S1248.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D13S317
dat.sub <- dat[dat$target == "D13S317",]
pdf("D13S317.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D16S539
dat.sub <- dat[dat$target == "D16S539",]
pdf("D16S539.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D18S51
dat.sub <- dat[dat$target == "D18S51",]
pdf("D18S51.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D19S433
dat.sub <- dat[dat$target == "D19S433",]
pdf("D19S433.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D1S1656
dat.sub <- dat[dat$target == "D1S1656",]
pdf("D1S1656.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D22S1045
dat.sub <- dat[dat$target == "D22S1045",]
pdf("D22S1045.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D2S441
dat.sub <- dat[dat$target == "D2S441",]
pdf("D2S441.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D5S818
dat.sub <- dat[dat$target == "D5S818",]
pdf("D5S818.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D7S820
dat.sub <- dat[dat$target == "D7S820",]
pdf("D7S820.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAD
dat.sub <- dat[dat$target == "Penta_D",]
pdf("PENTAD.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAE
dat.sub <- dat[dat$target == "Penta_E",]
pdf("PENTAE.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TH01
dat.sub <- dat[dat$target == "TH01",]
pdf("TH01.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TPOX
dat.sub <- dat[dat$target == "TPOX",]
pdf("TPOX.hist.strique50.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

##########
# 25bp
##########
dat <- read.table("results_25bp.tsv.final", header=T, sep="\t")
# single repeats
## CSF1P0
dat.sub <- dat[dat$target == "CSF1PO",]
pdf("CSF1P0.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D10S1248
dat.sub <- dat[dat$target == "D10S1248",]
pdf("D10S1248.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D13S317
dat.sub <- dat[dat$target == "D13S317",]
pdf("D13S317.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D16S539
dat.sub <- dat[dat$target == "D16S539",]
pdf("D16S539.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D18S51
dat.sub <- dat[dat$target == "D18S51",]
pdf("D18S51.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D19S433
dat.sub <- dat[dat$target == "D19S433",]
pdf("D19S433.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D1S1656
dat.sub <- dat[dat$target == "D1S1656",]
pdf("D1S1656.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D22S1045
dat.sub <- dat[dat$target == "D22S1045",]
pdf("D22S1045.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D2S441
dat.sub <- dat[dat$target == "D2S441",]
pdf("D2S441.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D5S818
dat.sub <- dat[dat$target == "D5S818",]
pdf("D5S818.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## D7S820
dat.sub <- dat[dat$target == "D7S820",]
pdf("D7S820.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAD
dat.sub <- dat[dat$target == "Penta_D",]
pdf("PENTAD.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## PENTAE
dat.sub <- dat[dat$target == "Penta_E",]
pdf("PENTAE.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TH01
dat.sub <- dat[dat$target == "TH01",]
pdf("TH01.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

## TPOX
dat.sub <- dat[dat$target == "TPOX",]
pdf("TPOX.hist.strique25.pdf")
ggplot(data=dat.sub, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat.sub$count)), max(as.integer(dat.sub$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

