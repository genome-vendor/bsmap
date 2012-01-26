import sys, time, os, array, optparse
usage = "usage: %prog [options] BSMAP_MAPPING_FILES"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output file name. (required)", default="")
parser.add_option("-d", "--ref", dest="reffile", metavar="FILE", help="reference genome fasta file. (required)", default="")
parser.add_option("-c", "--chr", dest="chroms", metavar="CHR", help="process only specified chromosomes, separated by ','. [default: all]\nexample: --chroms=chr1,chr2", default=[])
parser.add_option("-s", "--sam-path", dest="sam_path", metavar="PATH", help="path to samtools. [default: none]", default='')
parser.add_option("-u", "--unique", action="store_true", dest="unique", help="process only unique mappings/pairs.", default=False)                                    
parser.add_option("-p", "--pair", action="store_true", dest="pair", help="process only properly paired mappings.", default=False)
parser.add_option("-z", "--zero-meth", action="store_true", dest="meth0", help="report loci with zero methylation ratios.", default=False)
parser.add_option("-q", "--quiet", action="store_true", dest="quiet", help="don't print progress on stderr.", default=False)

options, infiles = parser.parse_args()

if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(options.outfile) == 0: parser.error("Missing output file name, use -o or --out option.")
if len(infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if any(options.chroms): options.chroms = options.chroms.split(',')

if any(options.sam_path): 
    if options.sam_path[-1] != '/': options.sam_path += '/'

def disp(txt, nt=0):
    if options.quiet: return 
    print >> sys.stderr, ''.join(['\t' for i in xrange(nt)]+['@ ',time.asctime(),': ',txt])

def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        flag = col[1] 
        if 'u' in flag: return []
        if options.unique and 's' in flag: return []
        if options.pair and 'P' not in flag: return []
        cr, pos, seq, strand, insert = col[2], int(col[3])-1, col[9], '', int(col[8])
        if cr not in options.chroms: return []
        if insert > 0:
            pos2 = int(col[7])-1
            seq = seq[:pos2-pos1]
        for aux in col[11:]:
            if aux[:5] == 'ZS:Z:': 
                strand = aux[5]
                break
        if strand == '': raise ValueError
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if options.unique and flag != 'UM': return []
        if options.pair and col[7] == '0': return []
        seq, strand, cr, pos = col[1], col[6][0], col[4], int(col[5])-1
        if cr not in options.chroms: return []
    return (seq, strand, cr, pos)

meth, depth = {}, {}
ref, cr, seq = {}, '', ''
disp('reading reference %s ...' % options.reffile)
for line in open(options.reffile):
    if line[0] == '>': 
        if any(cr): 
            if len(options.chroms) == 0 or cr in options.chroms: 
                ref[cr] = seq.upper()
                meth[cr] = array.array('I', [0]) * len(seq)
                depth[cr] = array.array('I', [0]) * len(seq)
        cr, seq = line[1:-1], ''
    else: seq += line.strip()

if len(options.chroms) == 0 or cr in options.chroms: 
    ref[cr] = seq.upper()
    meth[cr] = array.array('I', [0]) * len(seq)
    depth[cr] = array.array('I', [0]) * len(seq)

options.chroms = set(ref.keys())
del seq

BS_conversion = {'+': ('C','T'), '-': ('G','A')}
for infile in infiles:
    nline = 0
    disp('reading %s ...' % infile)
    if infile[-4:].upper() == '.SAM': sam_format, fin = 1, os.popen('%ssamtools view -XS %s' % (options.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = 1, os.popen('%ssamtools view -X %s' % (options.sam_path, infile))
    else: sam_format, fin = 0, open(infile)
    for line in fin:
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info
        readlen = len(seq)
        depthcr = depth[cr]
        if pos + readlen > len(depthcr): continue
        methcr = meth[cr]
        refseq = ref[cr][pos:pos+readlen]
        match, convert = BS_conversion[strand]
        index = refseq.find(match)
        while index >= 0:
            depthcr[pos+index] += 1
            if seq[index] == match: methcr[pos+index] += 1
            index = refseq.find(match, index+1)    
    
    fin.close()

disp('writing %s ...' % options.outfile)
ss = {'C': '+', 'G': '-'}
fout = open(options.outfile, 'w')
fout.write('chr\tpos\tstrand\tcontext\tratio\ttotal_C\tmethy_C\n')
for cr in sorted(depth.keys()):
    depthcr = depth[cr]
    methcr = meth[cr]
    refcr = ref[cr]
    for i, d in enumerate(depthcr):
        if d == 0: continue
        m = methcr[i]
        if m == 0 and not options.meth0: continue
        ratio = float(m) / d
        seq = refcr[i-2:i+3]
        strand = ss[refcr[i]]
        fout.write('%s\t%d\t%c\t%s\t%.3f\t%d\t%d\n' % (cr, i+1, strand, seq, ratio, d, m))

fout.close()
disp('Done!')
