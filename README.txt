BSMAP 1.0

1. Introduction

BSMAP is a short reads mapping program for bisulfite sequencing in DNA methylation study.  Bisulfite treatment coupled with next generation sequencing could estimate the methylation ratio of every single Cytosine location in the genome by mapping high throughput bisulfite reads to the reference sequences.

Bisulfite mapping is different from usual sequence mapping in two aspects: 1) The additional C/T mapping is asymmetric, a T in the read could be aligned to C in the reference but not vice versa. 2) The Watson and Crick strand are not complimentary after bisulfite treatment.  Each read need to be compared with 4 reference sequences, namely BSW(bisulfite Watson), BSWC(reverse complimentary of BSW), BSC(bisulfite Crick) and BSCC(reverse complimentary of BS).

BSMAP is designed to be a general-purpose mapping program to handle these special characteristics of bisulfite mapping.  It is based on the open source program SOAP (Short Oligo Alignment Program).  

Main features: 
	read length up to 144nt, allow up to 15 mismatches
	support pair end mapping, support parallel mapping, support SAM format inout/output
    support both whole genome (WGBS) and reduced representation bisulfite sequencing (RRBS)
    support trimming adapters and low quality sequences from 3'end of reads
    allow different running modes with flexible speed/memory/sensitivity to run on different hardware configurations
    include script to extract methylation ratios

BSMAP is under GNU Public License (GPL).


2. Installation

BSMAP is designed for linux64 platform. 
First unpackage the source code:
    $ tar zxfv bamsp-2.4.tgz

Make executable binary:
    $ make

Install the binary into system default path: (optional)
    $ make install
  
The following parameters could be modified in the makefile to achieve better mapping_speed/memory_usage under different situations:

OLIGOLEN: max read length, options: -DREAD_48, -DREAD_80, -DREAD_144(default)
	smaller max readlen will be faster.


3. Usage

bsmap <option>

option:
-a  <str>   query file, FASTA/FASTQ/BAM format.  The input format will be auto-detected. (required)
-b  <str>   query file b for pair end data, FASTA/FASTQ/BAM format.  The input format will be auto-detected. 
            if the input format is in BAM format, it should be the same as the file specified by "-a" option.
            BSMAP will read the two sets of reads w.r.t to the 0x40/0x80 flag in the input BAM file. 
            (required for pair-end mapping)
-d  <str>   reference sequences file, FASTA format. (required)
-o  <str>   output alignment file, if filename has .sam suffix, the output
            will be in SAM format, if the filename has .bam suffix, the output file
            be in sorted BAM file, and a filename.bai index file will be generated, 
            for other filename suffix the output is in BSP format. (required)
-2  <str>   output alignment file for unpaired reads in pair end mapping, only used for BSP format output. 
            If the output format is specified in BAM/SAM format, this option will be ignored, all alignments will be 
            writen to one BAM/SAM output file specified by the "-o" option. 
            (required for pair-end mapping with BSP format output)
-s  <int>   seed size, default=16, min=8, max=16. (WGBS mode)
            For RRBS mode, seed length is fixed to 12 and this command line option is neglected.
            longer seed size is faster, ~1.5 times faster with each additional nt
-v  <int>   max number of mismatches allowed on a read, default=2, max=15, 
            usually this number should be around 10% of the read length.
-w  <int>   max number of equal best hits to count, smaller will be faster, default=MAXHITS in makefile
-q  <int>   quality threshold in trimming 3'end of reads, 0-40, default=0. (no trim)
-z  <int>   base quality, default=33 [Illumina is using 64, Sanger Institute is using 33]
-f  <int>   filter low-quality reads containing >n Ns, default=5
-p  <int>   number of processors to use, default=1. The parallel performance scales well with 8 threads or less.
            For more than 8 threads, there might be no significant overall speed gain.
-x  <int>   max insertion size for pair end mapping, default=500
-m  <int>   min insertion size for pair end mapping, default=28
-L  <int>   mapping the first N nucleotide of the read, default: 0 (map the whole read).
-I  <int>   index interval (1~16), meaning the reference genome will be indexed every Nbp, default=4. (WGBS mode)
            For RRBS mode, index_interval is fixed to 1bp and this command line option is neglected.
            larger index interval uses memory, and slightly reduces mapping sensitivity. (~0.5% difference) 
            for human genome, -I 16 uses ~5GB, compared with ~9GB at the default -I 4.
-A  <str>   set the adapter sequence(s) and trim from 3'end of reads, default=none, requires at least 4nt matched, no mismatch allowed.
            Multiple -A options could be specified to set more than one adapter sequences, i.e. in pair-end sequencing case. 
            default: none (no adapter trimming)
-R          include the reference sequences as the XR:Z:<string> field in SAM output. default=do not include.
-B  <int>   start from the Nth read or read pair, default: 1.
-E  <int>   end at the Nth read or read pair, default: 4,294,967,295.
            Using -B and -E options user can specify part of the input file to be mapped, so that the input file 
            could be divided into several parts and mapped parallely over distributed system, without creating temporary files. 
-D  <str>   set restriction enzyme digestion site and activate reduced representation bisulfite mapping mode (RRBS mode), 
            i.e. reads must be mapped to digestion sites, the digestion site must be palindromic, digestion position is marked by '-', 
            for example: '-D C-CGG' for MspI digestion.
            default: none, meaning whole genome shot gun mapping (WGBS mode).
-S  <int>   seed for random number generation in selecting multiple hits.  default: 0 (seed set from system clock).
            other seed values generate pseudo random number based on read index number, so that mapping results are reproducible. 
-n  [0,1]   set mapping strand information:
            -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+)    (i.e. the "Lister protocol")
            for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --. 
            -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --    (i.e. the "Cokus protocol")
            default: -n 0. Most bisulfite sequencing data is generated only from forward strands.
-M  <str>   set the alignment information for the additional nucleotide transition. <str> is in the form of two different nucleotides, 
            the first one in the reads could be mapped to the second one in the reference sequences.
            default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion.
            example: -M GA could be used to detect to A=>I(G) transition in RNA editing. 
-h          help


4. Output

4.1 BSP format, includes the following tab delimited fields:
id, seq, qual map_flag, ref, ref_loc, strand, ins_size, refseq, #mismatches, mismatches_info

    1) id: read ID
    2) seq: mapped read sequence
    3) map_flag: 
        UM: unique map (unique pair for paired mapping).
        MA: multiple map (multiple pair for paired mapping)
        OF: over map (#multiple map exceeds MAXHITS)
        NM: no map
        QC: low quality reads
    4) ref: reference sequence name
    5) ref_loc: mapping location(1 based, 5'-end coordinates of the mapping region on the Watson strand of reference)
    6) strand: 
        ++: forward strand of Watson of reference (BSW)
        +-: reverse strand of Watson of reference (BSWC)
        -+: forward strand of Crick of reference (BSC)  
        --: reverse strand of Crick of reference (BSCC) 
    7) ins_size: insertion size for pair-end mapping, measured by the total nucleotide of the pair-end segment. 
        (edge to edge size). positive or negative value means mate's coordinate is larger or smaller respectively. 
        0 means single-end or unpaired mapping.
    8) refseq: Waston reference sequence at the mapping location.
    9) #mismatches: number of mismatches of current hit
    10) mismatch_info:  #hits of 0 mismatch to #hits of max_mismatches, separated by ':'

4.2 SAM format
    FLAG field: 
        UM: 0x0
        MA: 0x100
        OF: 0x100
        NM: 0x4
        QC: 0x204
        for mapping on BSC or BSCC: FLAG=FLAG+0x10
	for pair-end mapping:
	    FLAG=FLAG+0x1
            if it's the first read in pair, FLAG=FLAG+0x40
            if it's the second read in pair, FLAG=FLAG+0x80
	    if mappings are paired, FLAG=FLAG+0x2
            if mate is unmapped, FLAG=FLAG+0x8
            if mate is mapped on BSC or BSCC, FLAG=FLAG+0x20
    aux field: 
        ZS:Z:<strand info> same as BSP column 6).
        XR:Z:<reference sequence> same as BSP column 8).
        NM:i:<#mismatches> same as BSP column 9).    
        ZP:i:<int> RRBS fragment start location, only for RRBS mode.
        ZL:i:<int> RRBS fragment size, only for RRBS mode.
	

    for more details, please refer to SAM format specification: 
    http://samtools.sourceforge.net/SAM1.pdf

**Note: all read sequences are recorded as the corresponding sequence following the reference Watson strand direction.

 

5. Speed and sensitivity
    
The longer seed size(option -s), the faster speed. With seed size increase every bp, mapping time reduces by ~1.5-fold. 
On the other hand, the max number of mismatches that could be detected with 100% sensitivity is bounded by the 
seed_size.  

	max_mismatches_with_100%_sensitivity = (read_len+1-index_interval) / seed_size - 1

If the -v option set max mismatches larger than this number, those mappings with larger max mismatches may not be 
guaranteed to be detected. 

In case full sensitivity can not be achieved within feasible time, user will need to make a decision on the trade off 
between the speed and sensitivity by setting the optimal seed size.  


6 Example

for shot gun whole genome BS: 

	single_end: (-A means trimming adapter from 3'end)
	$ bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100  -v 5 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA

	single_end: (map to all 4 possible strands)
	$bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -n 1 -w 100 -v 5

	pair_end: (set -b option)
        $ bsmap -a read1.fq -b read2.fq -d ~/ref/hg19/hg19.fa -o out_pair.bsp -2 out_upair.bsp -p 8 -w 100  -v 5 
	$ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100  -v 5 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA

	using less memory: (set -I option)
        $ bsmap -a SE_read.bam -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -v 5 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA -I 8      

	mapping from read pair #10001 to read pair #20000 in the input file: (set -B and -E option)
	$ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -w 100  -v 5 -B 10001 -E 20000

	using Illumina quality: (set -z option)
	$ bsmap -a PE_read1.fq -b PE_read2.fq -d ~/ref/hg19/hg19.fa -o out.bam -z 64

	trimming low quality 3'end: (set -q option)
	$ bsmap -a PE_reads.bam -b PE_reads.bam -d ~/ref/hg19/hg19.fa -o out.bam -q 2

	detect A=>G editing in RNA_seq instead of C=>T conversion in bisulfite sequencing
	$ bsmap -a reads.bam -d RNA_ref.fa -M GA -o out.bam

for RRBS: (set -D option to specify digestion site information and activate RRBS mode.)

	bsmap -a PE_reads.bam -b PE_reads.bam  -d ~/ref/hg19/hg19.fa -o out.bam -p 8 -w 100 -s 12 -v 5 -D C-CGG


7. Scripts: 

7.1 methratio.py
python script to extract methylation ratios from BSMAP mapping results. Require python 2.X. 
For human genome, methratio.py needs ~26GB memory.  
For systems with limited memory, user can set the -c/--chr option to process specified chromosomes only,
and combine results for all chromosomes afterwards.

Usage: python methratio.py [options] BSMAP_MAPPING_FILES

Options:
  -h, --help            show this help message and exit
  -o FILE, --out=FILE   output file name. (required)
  -d FILE, --ref=FILE   reference genome fasta file. (required)
  -c CHR, --chr=CHR     process only specified chromosomes. [default: all]
                        example: --chr=chr1,chr2 (this uses ~4.5GB compared with ~26GB for the whole genome)
  -s PATH, --sam-path=PATH
                        path to samtools. [default: none]
  -u, --unique          process only unique mappings/pairs.
  -p, --pair            process only properly paired mappings.
  -z, --zero-meth       report loci with zero methylation ratios.
  -q, --quiet           don't print progress on stderr.

Output format: tab delimited txt file with the following columns:
    1) chromorome
    2) coordinate (1-based)
    3) strand
    4) sequence context (2nt upstream to 2nt downstream in Watson strand direction)
    5) methylation ratio
    6) number of reads covering this locus 
    7) number of unconverted Cs in the reads at this locus

Example:
	python methratio.py --chr=chr1,chr2 --ref=hg19.fa --out=methratio.txt rrbsmap_sample*.sam
    python methratio.py -d mm9.fa -o meth.txt -p bsmap_sample1.bsp bsmap_sample2.sam bsmap_sample3.bam 

Note: For overlapping paired hits, nucleotides in the overlapped part should be counted only once instead of twice.
methratio.py can correctly handle such cases for SAM format output, but for BSP format it will still be counted twice,
because the BSP format does not contain mapping information of the mate.


7.2 sam2bam.sh
Shell script to convert SAM format to sorted and indexed BAM format.
The input SAM file will be deleted if the conversion is successful.
This script is automatically called by BSMAP if the output file has .bam suffix.  It can also be used manually.

Usage: ./sam2bam.sh INPUT_SAM_FILE
example: ./sam2bam.sh sample1.sam 
This will generate sorted BAM file sample1.bam and index file sample1.bam.bai. 



8. Citation
    Yuanxin Xi and Wei Li, "BSMAP: whole genome bisulfite sequence MAPping program" (2009) BMC Bioinformatics 2009, 10:232


9. Contact
    Yuanxin Xi
    Bioinformatics Division, 
    Dan L. Duncan Cancer Center,
    Baylor College of Medicince, 
    Houston, TX 77030, USA
    713-798-6254
    yxi@bcm.tmc.edu, xiyuanxin@gmail.com
    
