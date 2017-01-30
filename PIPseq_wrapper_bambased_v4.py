#! /usr/bin/python2.7
# Wrapper will take .bam alignments for a PIP-seq tetrad
# must have .bam as follows ./[project_dir]/data/topOut/tophat_output_RD_[sample_id]/NR_hits.bam
# must have scripts in ./[project_dir]/scripts/
# run wrapper in ./[project_dir]/data/topOut
# will generate directories (default: CSAR_output/perms) in ./[project_dir]/data/
# use a tag shared by only the 4 PIP-seq samples to run
# Follow file naming format:'[unique_tag]_(d/s)sRNase_(no/)protein_((un/)trimmed/merged)/NR_hits.bam'
# Examples: 'Pp_dark_ssRNase_protein_trimmed/NR_hits.bam', 'totalArab_nuc_rep4_dsRNase_noprotein_merged/NR_hits.bam'
# Ideally, you will use this in conjunction with 'alignment_wrapper.py', which will reduce naming issues
# This is a wrapper, so you need the other scripts to actually do the work. Here's the list:
# bedtools (suite), 'shuffle_reads_BAM.pl', 'split_coverage_by_chrom.pl', 'generate_empty_coverage_files.rb', 'make_CSAR_files.R', 'run_CSAR_shuffled.R', 'run_CSAR_saturation.R'
# all scripts must be in same dir as wrapper script
# !!! DO NOT perform runs with identical 'alignment_tag' but different 'trimming_status'
# This wrapper is an organized and slightly edited version of the PIP-seq analysis piple developed by Fan Li

## Set up python environment
import sys
import os
import argparse
import subprocess
import re
import glob
import random

## Process arguments
parser = argparse.ArgumentParser(description="process 4 .bam files from a PIP-seq experiment, vi and check comments for more info, Follow file naming format for tophat output directories:'tophat_output_RD_[unique_tag]_(d/s)sRNase_(no/)protein_((un/)trimmed/merged)', Examples: 'Pp_dark_ssRNase_protein_trimmed/', 'totalArab_nuc_rep4_dsRNase_noprotein_merged/'")
parser.add_argument('alignment_tag', help='common tag for 4 .tbl files of PIP-seq')
parser.add_argument('trimming_status', help="specify if using a 'trimmed', 'untrimmed', or 'merged' set of reads")
parser.add_argument('chr_len', help="specify path to chromosome length or contig length file")
parser.add_argument('--dsRNase_only','-do', action='store_true', help='only perform PPS analysis for dsRNase')
parser.add_argument('--ssRNase_only','-so', action='store_true', help='only perform PPS analysis for ssRNase')
parser.add_argument('--specify_out_dir','-spout', help="specify output dir (need path, use trailing '/')")
parser.add_argument('--merge_reads','-mrg',action='store_true',help="call to merge timmed an untrimmed reads using wrapper, will make tophat_output_*_merged directories for you, if not using this must have sorted bam files (named 'NR_hits.bam')")
parser.add_argument('--filter_bam','-f',help="input file to filter out reads in specified regions, should be .bed file of unwanted regions, specify entire path")
#args=parser.parse_args('Pp_dark trimmed /Data03/sgosai/annotation/pp/Ppatens_152_chr_len.txt'.split())
args=parser.parse_args()

## Check organizational requirements for correct data processing

# error out if ds and ss only flags are both raised
if args.dsRNase_only and args.ssRNase_only:
	raise ValueError("either ds_only or ss_only can be raised, if you want both, exclude these args")

# raise error if -spin is called without -spout:
#if args.specify_in_dir and not args.specify_out_dir:
#	raise ValueError("must specify output directory if specifying input directory")

# raise error if -mrg is called when using trimmed or untrimmed reads (i.e. trimming_status != merged)
if args.merge_reads and not args.trimming_status == "merged":
	raise NameError("only merge trimmed and untrimmed bam files if doing analysis on merged")

# check that path is not included in input arg
if "/" in args.alignment_tag:
	raise NameError("do not include .tbl file path or '/' characters, run script in directory containing tophat_output_RD directories")

# if merging reads with wrapper, make tophat_output dir for merged .bam files
if args.merge_reads:
	mrgDir=[]
	if not args.dsRNase_only:
		mrgDir+=['tophat_output_RD_'+args.alignment_tag+'_ssRNase_noprotein_merged','tophat_output_RD_'+args.alignment_tag+'_ssRNase_protein_merged']
	if not args.ssRNase_only:
		mrgDir+=['tophat_output_RD_'+args.alignment_tag+'_dsRNase_noprotein_merged','tophat_output_RD_'+args.alignment_tag+'_dsRNase_protein_merged']
	for make_this in mrgDir:
		subprocess.check_call(['mkdir',make_this])

# check that desired PIP-seq tophat_output directories are present
bamDirList=os.listdir('./')
bamDirStr='\n'.join(os.listdir('./'))
if args.dsRNase_only:
	findER=re.compile('^tophat_output_RD_'+args.alignment_tag+'_dsRNase_.*protein_'+args.trimming_status,re.MULTILINE)
	expected_count=2
elif args.ssRNase_only:
	findER=re.compile('^tophat_output_RD_'+args.alignment_tag+'_ssRNase_.*protein_'+args.trimming_status,re.MULTILINE)
	expected_count=2
elif not args.dsRNase_only and not args.ssRNase_only:
	findER=re.compile('^tophat_output_RD_'+args.alignment_tag+'_[ds]sRNase_.*protein_'+args.trimming_status,re.MULTILINE)
	expected_count=4

inList=findER.findall(bamDirStr)
print inList
if len(inList) != expected_count:
	raise ValueError("unexpected number of tophat_output_RD directories match specifications defined by args, make sure tophat_output_RD directories are named as follows: ./[unique_sample_tag]_(d/s)sRNase_(no/)protein_((un/)trimmed/merged)/, --help for more info")

# set scripts directory
scrPth=sys.argv[0].split('/')
scrDir=('/'.join(scrPth[:(len(scrPth)-1)]))+'/'

# check for required scripts
reqScr=['shuffle_reads_BAM.pl', 'split_coverage_by_chrom.pl', 'generate_empty_coverage_files.rb', 'make_CSAR_files.R', 'run_CSAR_shuffled.R', 'run_CSAR_saturation.R']
havScr=os.listdir(scrDir)
getScr=[]
for needThis in reqScr:
	if needThis not in havScr:
		getScr.append(needThis)

if len(getScr) != 0:
	for getThis in getScr:
		print "(!) Need '%s'" % getThis
	raise ValueError("Missing reqired scripts for analysis, move indicated scripts to same directory as this script")

# generate output directory
if not args.specify_out_dir:
	if not os.path.isdir("../CSAR_output/perms"):
		subprocess.call(["mkdir", "-p", "../CSAR_output/perms"])
	outDir="../CSAR_output/"
	outPrm="../CSAR_output/perms/"
else:
	if args.specify_out_dir[len(args.specify_out_dir)-1] != '/':
		raise NameError("Need to include trailing '/' specified output arg")
	if not os.path.isdir(args.specify_out_dir):
		subprocess.check_call(["mkdir", args.specify_out_dir])
		subprocess.check_call(["mkdir", args.specify_out_dir+'perms'])
	outDir=args.specify_out_dir
	outPrm=args.specify_out_dir+'perms/'

## Begin data processing

# Notify user run is begining
print "Commence Death"

# Get list of leading tages:
bamHash={}
for item in inList:
	k=item.replace("tophat_output_RD_","").replace("_"+args.trimming_status,"")
	bamHash[k]=item

# merge and sort .bam files if prompted
if args.merge_reads:
	for tag,topDir in bamHash.items():
		rTag=str(random.randint(10000,99999))
		subprocess.check_call(['samtools','merge',topDir+'/tmp'+rTag+'.bam',topDir.replace("merged","trimmed")+'/NR_hits.bam',topDir.replace("merged","untrimmed")+'/NR_hits.bam'])
		subprocess.check_call(['samtools','sort',topDir+'/tmp'+rTag+'.bam',topDir+'/NR_hits'])
		os.remove(topDir+'/tmp'+rTag+'.bam')

# set target .bam files and filter if necessary
if not args.filter_bam:
	set_bam='NR_hits'
else:
	set_bam='NR_hits.filtered'
	for tag,topDir in bamHash.items():
		rTag=str(random.randint(10000,99999))
		filtered_bam=open(topDir+"/tmp"+rTag+".bam",'w')
		subprocess.check_call(['intersectBed','-a',topDir+'/NR_hits.bam','-b',args.filter_bam, '-s','-v'],stdout=filtered_bam)
		filtered_bam.close()
		subprocess.check_call(['samtools','sort',topDir+"/tmp"+rTag+".bam",topDir+'/'+set_bam])
		os.remove(topDir+"/tmp"+rTag+".bam")

# Shuffle .tbl files
if not args.ssRNase_only:
	subprocess.check_call(['perl',scrDir+'shuffle_reads_BAM.pl','tophat_output_RD_'+args.alignment_tag+'_dsRNase_protein_'+args.trimming_status+'/'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_dsRNase_noprotein_'+args.trimming_status+'/'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_dsRNase_protein_'+args.trimming_status+'/shuffled_'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_dsRNase_noprotein_'+args.trimming_status+'/shuffled_'+set_bam+'.bam'])

if not args.dsRNase_only:
	subprocess.check_call(['perl',scrDir+'shuffle_reads_BAM.pl','tophat_output_RD_'+args.alignment_tag+'_ssRNase_protein_'+args.trimming_status+'/'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_ssRNase_noprotein_'+args.trimming_status+'/'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_ssRNase_protein_'+args.trimming_status+'/shuffled_'+set_bam+'.bam','tophat_output_RD_'+args.alignment_tag+'_ssRNase_noprotein_'+args.trimming_status+'/shuffled_'+set_bam+'.bam'])

# Calculate genome coverage for shuffled
for tag,topDir in bamHash.items():
	plusCov=open(outPrm+'shuffled_'+tag+'.plus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', topDir+'/shuffled_'+set_bam+'.bam', '-d', '-split','-strand',  '+'], stdout=plusCov)
	plusCov.close()
	subprocess.check_call(['perl',scrDir+'split_coverage_by_chrom.pl', outPrm+'shuffled_'+tag+'.plus.coverage.txt', outPrm+'shuffled_'+tag+'.plus'])
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.plus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+'generate_empty_coverage_files.rb', outPrm+'shuffled_'+tag+'.plus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+'make_CSAR_files.R', outPrm, 'shuffled_'+tag+'.plus', 'Forward'])
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	minusCov=open(outPrm+'shuffled_'+tag+'.minus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', topDir+'/shuffled_'+set_bam+'.bam', '-d', '-split','-strand',  '-'], stdout=minusCov)
	minusCov.close()
	subprocess.check_call(['perl',scrDir+'split_coverage_by_chrom.pl', outPrm+'shuffled_'+tag+'.minus.coverage.txt', outPrm+'shuffled_'+tag+'.minus'])
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.minus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+'generate_empty_coverage_files.rb', outPrm+'shuffled_'+tag+'.minus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+'make_CSAR_files.R', outPrm, 'shuffled_'+tag+'.minus', 'Reverse'])
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR shuffled
if not args.ssRNase_only:
	subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+'run_CSAR_shuffled.R', outPrm+'shuffled_'+args.alignment_tag+'_dsRNase',args.chr_len, outPrm+'shuffled_'+args.alignment_tag+'_dsRNase_'+args.trimming_status+'.CSAR_PPS.bed', outPrm+'shuffled_'+args.alignment_tag+'_dsRNase_'+args.trimming_status+'.CSAR_PPS.threshold.txt'])
	doneFiles=glob.glob(outPrm+'shuffled_'+args.alignment_tag+'_dsRNase'+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+args.alignment_tag+'_dsRNase'+'*CSARScore')
	for target in doneFiles:
		os.remove(target)
if not args.dsRNase_only:
	subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+'run_CSAR_shuffled.R', outPrm+'shuffled_'+args.alignment_tag+'_ssRNase',args.chr_len, outPrm+'shuffled_'+args.alignment_tag+'_ssRNase_'+args.trimming_status+'.CSAR_PPS.bed', outPrm+'shuffled_'+args.alignment_tag+'_ssRNase_'+args.trimming_status+'.CSAR_PPS.threshold.txt'])
	doneFiles=glob.glob(outPrm+'shuffled_'+args.alignment_tag+'_ssRNase'+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+args.alignment_tag+'_ssRNase'+'*CSARScore')
	for target in doneFiles:
		os.remove(target)

# Calculate genome coverage for true
for tag,topDir in bamHash.items():
	plusCov=open(outDir+tag+'.plus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', topDir+'/'+set_bam+'.bam', '-d', '-split', '-strand', '+'], stdout=plusCov)
	plusCov.close()
	subprocess.check_call(['perl',scrDir+'split_coverage_by_chrom.pl', outDir+tag+'.plus.coverage.txt', outDir+tag+'.plus'])
	subprocess.check_call(['rm', outDir+tag+'.plus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+'generate_empty_coverage_files.rb', outDir+tag+'.plus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+'make_CSAR_files.R', outDir, tag+'.plus', 'Forward'])
	doneFiles=glob.glob(outDir+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	minusCov=open(outDir+tag+'.minus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', topDir+'/'+set_bam+'.bam', '-d', '-split', '-strand', '-'], stdout=minusCov)
	minusCov.close()
	subprocess.check_call(['perl',scrDir+'split_coverage_by_chrom.pl', outDir+tag+'.minus.coverage.txt', outDir+tag+'.minus'])
	subprocess.check_call(['rm', outDir+tag+'.minus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+'generate_empty_coverage_files.rb', outDir+tag+'.minus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+'make_CSAR_files.R', outDir, tag+'.minus', 'Reverse'])
	doneFiles=glob.glob(outDir+tag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR on true data
if not args.ssRNase_only:
	subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+'run_CSAR_saturation.R',outDir+args.alignment_tag+'_dsRNase',args.chr_len,outPrm+'shuffled_'+args.alignment_tag+'_dsRNase_'+args.trimming_status+'.CSAR_PPS.threshold.txt',outDir+args.alignment_tag+'_dsRNase_'+args.trimming_status+'.CSAR_PPS.bed',outDir+args.alignment_tag+'_dsRNase_'+args.trimming_status+'.CSAR_PPS_counts.txt'])
	doneFiles=glob.glob(outDir+args.alignment_tag+'_dsRNase'+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+args.alignment_tag+'_dsRNase'+'*CSARScore')
	for target in doneFiles:
		os.remove(target)

if not args.dsRNase_only:
	subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+'run_CSAR_saturation.R',outDir+args.alignment_tag+'_ssRNase',args.chr_len,outPrm+'shuffled_'+args.alignment_tag+'_ssRNase_'+args.trimming_status+'.CSAR_PPS.threshold.txt',outDir+args.alignment_tag+'_ssRNase_'+args.trimming_status+'.CSAR_PPS.bed',outDir+args.alignment_tag+'_ssRNase_'+args.trimming_status+'.CSAR_PPS_counts.txt'])
	doneFiles=glob.glob(outDir+args.alignment_tag+'_ssRNase'+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+args.alignment_tag+'_ssRNase'+'*CSARScore')
	for target in doneFiles:
		os.remove(target)

print "Finished run, output can be found in "+outDir
