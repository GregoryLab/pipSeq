
library("CSAR")
library("gplots")

args <- commandArgs(T)
if (length(args) < 5) {
	cat("USAGE: ./run_CSAR_saturation.R tag chr_len_fn thresholds_fn out_fn out_counts_fn\n")
	q()
}
tag <- args[1]
chr_len_fn <- args[2]
thresholds_fn <- args[3]
out_fn <- args[4]
out_counts_fn <- args[5]

# read in chromosome lengths
chr_len <- read.table(chr_len_fn, header=F, as.is=T, sep="\t")
colnames(chr_len) <- c("chr", "len")

# generate nhits data structures
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(paste(tag, "protein", sep="_"), "plus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsS.forward <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(paste(tag, "protein", sep="_"), "minus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsS.reverse <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(paste(tag, "noprotein", sep="_"), "plus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsC.forward <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)
##
chr <- chr_len$chr
chrL <- chr_len$len
chrL_0 <- vector("integer",length(chr))
filenames<-vector("character",length(chr))
c1 <- vector("integer",length(chr))
c2 <- vector("integer",length(chr))
for (i in 1:length(chr)) {
	filenames[i] <- paste(paste(tag, "noprotein", sep="_"), "minus", chr[i], "CSARNhits", sep=".")
	values <- LoadBinCSAR(filenames[i])
	chrL_0[i] <- length(which(values>0))
	c1[i] <- sum(as.numeric(values))
	c2[i] <- sum(as.numeric(values)^2)
}
nhitsC.reverse <- list(chr=chr,chrL=chrL,chrL_0=chrL_0,filenames=filenames,c1=c1,c2=c2)

norm.forward <- sum(as.numeric(nhitsC.forward$c1))
norm.reverse <- sum(as.numeric(nhitsC.reverse$c1))
# run CSAR
#score.forward <- ChIPseqScore(control=nhitsC.forward,sample=nhitsS.forward,file=sprintf("%s.plus",tag),norm=norm.forward,times=10000,test="Poisson")
#score.reverse <- ChIPseqScore(control=nhitsC.reverse,sample=nhitsS.reverse,file=sprintf("%s.minus",tag),norm=norm.reverse,times=10000,test="Poisson")
score.forward <- ChIPseqScore(control=nhitsC.forward,sample=nhitsS.forward,norm=norm.forward,times=10000,test="Poisson")
win.forward<-sigWin(score.forward, t=1, g=1)

score.reverse <- ChIPseqScore(control=nhitsC.reverse,sample=nhitsS.reverse,norm=norm.reverse,times=10000,test="Poisson")
win.reverse<-sigWin(score.reverse, t=1, g=1)
win <- c(win.forward, win.reverse)

# read in FDR thresholds (from shuffled data)
data.thresholds <- read.table(thresholds_fn, header=F, as.is=T, sep="\t")
colnames(data.thresholds) <- c("FDR", "threshold", "shuffled_count")

# print out results to BED file
results <- data.frame(as.vector(seqnames(win)))
colnames(results) <- "chrom"
results$start <- start(ranges(win))-1
results$end <- end(ranges(win))
results$id <- unlist(lapply(1:length(win), function(x) sprintf("CSAR_peak_%d",x)))
results$score <- score(win)
results$strand <- c(rep("+",length(win.forward)), rep("-",length(win.reverse)))
write.table(results, file=out_fn, quote=F, sep="\t", row.names=F, col.names=F)

results.fdr05 <- subset(results, score>=data.thresholds$threshold[which(data.thresholds$FDR==0.05)])
write.table(results.fdr05, file=gsub(".bed$", "_FDR05.bed", out_fn), quote=F, sep="\t", row.names=F, col.names=F)

# count how many PPSs are real at each FDR threshold and print results
data.thresholds$PPS_count <- unlist(lapply(data.thresholds$threshold, function(x) length(which(results$score>=x))))
write.table(data.thresholds[,c("FDR","threshold","PPS_count")], file=out_counts_fn, quote=F, sep="\t", row.names=F, col.names=F)




