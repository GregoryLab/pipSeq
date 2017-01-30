
library("CSAR")
library("gplots")

args <- commandArgs(T)
if (length(args) < 3) {
	cat("USAGE: ./make_CSAR_files.R dir tag considerStrand\n")
	q()
}
dir <- args[1]
dir <- gsub("/$", "", dir)
tag <- args[2]
considerStrand <- args[3]

# get list of coverage files
files <- list.files(dir, pattern=sprintf("%s.*.coverage.txt$", tag), full.names=T)
#considerStrand <- "Forward"

# load each file and convert it to a binary *.CSARNhits file
for (i in 1:length(files)) {
	values <- as.integer(ceiling(scan(files[i]))) # wrapping as.integer after scanning for numeric handles overflow better (ie really large values get cropped to .Machine$integer.max)
	outfile <- gsub(".coverage.txt$", ".CSARNhits", files[i])
	properties <- unlist(strsplit(basename(outfile), "\\."))
	
	file1<-file(description=outfile,"wb")
	writeBin("CSARNhits",con=file1)
	writeBin("v.1",con=file1) ###Check how to put version here
	writeBin(considerStrand,con=file1) # considerStrand
	writeBin(as.integer(0),con=file1) # only write a placeholder value for w...otherwise reading from binary files is fucked up
	writeBin(FALSE,con=file1)	# uniquelyMapped
	writeBin(FALSE,con=file1)	# uniquePosition
	writeBin(properties[4],con=file1)	# chromosome
	writeBin(length(values),con=file1)	# length of values
	writeBin(values,con=file1)	# values
	close(file1)
}
