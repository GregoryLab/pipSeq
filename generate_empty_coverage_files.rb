#!/usr/bin/env ruby

# generate coverage files for chromosomes that did not get processed by genomeCoverageBed (ie if they have no coverage at all, no values get printed)
if ARGV.size == 2
	prefix, chr_len_fn = ARGV
else
	puts "USAGE: #{$0} prefix chr_len_fn"
	exit -1
end

# store chromosomes
chr_lens = Hash.new
File.open(chr_len_fn).each_line do |line|
	chr, len = line.chomp.split(/\t/)
	chr_lens[chr] = len
end

# get list of existing files
existing = Hash.new
files = Dir.glob("#{prefix}*.coverage.txt")
files.each { |file|
	name = File.basename(file)
#	subset, tag, strand, chrom, junk, junk = name.split(/\./)
	tag, strand, chrom, junk, junk = name.split(/\./)
	existing[chrom] = 1
}

# find the chromosomes that dont have existing files
chr_lens.each { |chr, len|
	len = len.to_i
	if !existing.has_key?(chr)
		outfn = "#{prefix}.#{chr}.coverage.txt"
		puts "missing #{chr} - replacing with empty coverage file #{outfn}"
		fp = File.new(outfn, "w")
		1.upto(len) do
			fp.puts "0"
		end
		fp.close		
	end
}
