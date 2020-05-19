#!/usr/bin/env Rscript
library(ggplot2)

# single repeats
## CSF1P0
dat <- read.table("CSF1P0.motifs", header=T, sep="\t")
pdf("CSF1P0.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D10S1248
dat <- read.table("D10S1248.motifs", header=T, sep="\t")
pdf("D10S1248.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D13S317
dat <- read.table("D13S317.motifs", header=T, sep="\t")
pdf("D13S317.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D16S539
dat <- read.table("D16S539.motifs", header=T, sep="\t")
pdf("D16S539.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D18S51
dat <- read.table("D18S51.motifs", header=T, sep="\t")
pdf("D18S51.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D19S433
dat <- read.table("D19S433.motifs", header=T, sep="\t")
pdf("D19S433.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D1S1656
dat <- read.table("D1S1656.motifs", header=T, sep="\t")
pdf("D1S1656.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D22S1045
dat <- read.table("D22S1045.motifs", header=T, sep="\t")
pdf("D22S1045.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D2S441
dat <- read.table("D2S441.motifs", header=T, sep="\t")
pdf("D2S441.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D5S818
dat <- read.table("D5S818.motifs", header=T, sep="\t")
pdf("D5S818.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## D7S820
dat <- read.table("D7S820.motifs", header=T, sep="\t")
pdf("D7S820.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## PENTAD
dat <- read.table("PENTAD.motifs", header=T, sep="\t")
pdf("PENTAD.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## PENTAE
dat <- read.table("PENTAE.motifs", header=T, sep="\t")
pdf("PENTAE.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## TH01
dat <- read.table("TH01.motifs", header=T, sep="\t")
pdf("TH01.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()
## TPOX
dat <- read.table("TPOX.motifs", header=T, sep="\t")
pdf("TPOX.hist.pdf")
ggplot(data=dat, aes(x=as.integer(count))) + 
	geom_bar(aes(y=..prop..)) + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(dat$count)), max(as.integer(dat$count),1))) + 
	xlab("Predicted repeat count") +
	ylab("Frequency")
dev.off()

# multiple motifs
## D12S391
dat <- read.table("D12S391.motifs", header=T, sep="\t")
AGAT <- as.factor(dat$count)
AGAC <- as.factor(dat$count.1)
reps <- c(rep("AGAT", length(AGAT)), rep("AGAC", length(AGAC)))
com <- c(AGAT, AGAC)
df <- data.frame(reps, com)
pdf("D12S391.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## D21S11
dat <- read.table("D21S11.motifs", header=T, sep="\t")
TCTA <- as.factor(dat$count)
TCTG <- as.factor(dat$count.1)
reps <- c(rep("TCTA", length(TCTA)), rep("TCTG", length(TCTG)))
com <- c(TCTA, TCTG)
df <- data.frame(reps, com)
pdf("D21S11.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## D2S1338
dat <- read.table("D2S1338.motifs", header=T, sep="\t")
GGAA <- as.factor(dat$count)
GGCA <- as.factor(dat$count.1)
reps <- c(rep("GGAA", length(GGAA)), rep("GGCA", length(GGCA)))
com <- c(GGAA, GGCA)
df <- data.frame(reps, com)
pdf("D2S1338.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## D3S1358
dat <- read.table("D3S1358.motifs", header=T, sep="\t")
TCTG <- as.factor(dat$count)
TCTA <- as.factor(dat$count.1)
reps <- c(rep("TCTG", length(TCTG)), rep("TCTA", length(TCTA)))
com <- c(TCTG, TCTA)
df <- data.frame(reps, com)
pdf("D3S1358.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## D8S1179
dat <- read.table("D8S1179.motifs", header=T, sep="\t")
TCTA <- as.factor(dat$count)
TCTG <- as.factor(dat$count.1)
reps <- c(rep("TCTA", length(TCTA)), rep("TCTG", length(TCTG)))
com <- c(TCTA, TCTG)
df <- data.frame(reps, com)
pdf("D8S1179.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## FGA
dat <- read.table("FGA.motifs", header=T, sep="\t")
GGAA <- as.integer(dat$count)
AAAG <- as.integer(dat$count.1)
GAAA <- as.integer(dat$count.2)
GAAG <- as.integer(dat$count.3)
reps <- c(rep("GGAA", length(GGAA)), rep("AAAG", length(AAAG)), rep("GGAA", length(GGAA)), rep("GAAG", length(GAAG)))
com <- c(GGAA, AAAG, GAAA, GAAG)
df <- data.frame(reps, com)
pdf("FGA.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()
## vWA
dat <- read.table("vWA.motifs", header=T, sep="\t")
TAGA <- as.factor(dat$count)
CAGA <- as.factor(dat$count.1)
reps <- c(rep("TAGA", length(TAGA)), rep("CAGA", length(CAGA)))
com <- c(TAGA, CAGA)
df <- data.frame(reps, com)
pdf("vWA.hist.pdf")
ggplot(df,aes(x=com, fill=reps, y=..prop..)) + 
	geom_bar(position='dodge') + 
	theme_minimal() + 
	scale_x_continuous(breaks = seq(min(as.integer(df$com)), max(as.integer(df$com),1))) + 
	xlab("Predicted repeat count") + 
	ylab("Frequency")
dev.off()







