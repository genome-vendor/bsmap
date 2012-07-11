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
parser.add_option("-r", "--remove-duplicate", action="store_true", dest="rm_dup", help="remove duplicated reads.", default=False)
parser.add_option("-t", "--trim-fillin", dest="trim_fillin", type="int", metavar='N', help="trim N end-repairing fill-in nucleotides. [default: 2]", default=2)
parser.add_option("-g", "--combine-CpG", action="store_true", dest="combine_CpG", help="combine CpG methylaion ratios on both strands.", default=False)
parser.add_option("-m", "--min-depth", dest="min_depth", type="int", metavar='FOLD', help="report loci with sequencing depth>=FOLD. [default: 1]", default=1)

options, infiles = parser.parse_args()

if len(options.reffile) == 0: parser.error("Missing reference file, use -d or --ref option.")
if len(options.outfile) == 0: parser.error("Missing output file name, use -o or --out option.")
if len(infiles) == 0: parser.error("Require at least one BSMAP_MAPPING_FILE.") 
if any(options.chroms): options.chroms = options.chroms.split(',')

if any(options.sam_path): 
    if options.sam_path[-1] != '/': options.sam_path += '/'

def disp(txt, nt=0):
    if not options.quiet: print >> sys.stderr, ''.join(['\t' for i in xrange(nt)]+['@ ',time.asctime(),': ',txt])

def get_alignment(line):
    col = line.split('\t')
    if sam_format:
        flag = col[1] 
        if 'u' in flag: return []
        if options.unique and 's' in flag: return []
        if options.pair and 'P' not in flag: return []
        cr, pos, seq, strand, insert = col[2], int(col[3])-1, col[9], '', int(col[8])
        if cr not in options.chroms: return []
        for aux in col[11:]:
            if aux[:5] == 'ZS:Z:': 
                strand = aux[5:7]
                break
        if strand == '': raise ValueError
    else:
        flag = col[3][:2]
        if flag == 'NM' or flag == 'QC': return []
        if options.unique and flag != 'UM': return []
        if options.pair and col[7] == '0': return []
        seq, strand, cr, pos, insert = col[1], col[6], col[4], int(col[5])-1, int(col[7])
        if cr not in options.chroms: return []
    if options.rm_dup:  # remove duplicate hits
        if strand == '+-' or strand == '-+': frag_end, direction = pos+len(seq), 2
        else: frag_end, direction = pos, 1
        if coverage[cr][frag_end] & direction: return []
        coverage[cr][frag_end] |= direction
    if options.trim_fillin > 0: # trim fill in nucleotides
        if strand == '+-': seq = seq[:-options.trim_fillin]
        elif strand == '--': seq, pos = seq[options.trim_fillin:], pos+options.trim_fillin
        elif insert != 0 and len(seq) > abs(insert) - options.trim_fillin:
            trim_nt = len(seq) - (abs(insert) - options.trim_fillin)
            if strand == '++': seq = seq[:-trim_nt]
            elif strand == '-+': seq, pos =seq[trim_nt:], pos+trim_nt
    if sam_format and insert > 0: seq = seq[:int(col[7])-1-pos] # remove overlapped regions in paired hits, SAM format only
    return (seq, strand[0], cr, pos)

ref, cr, seq = {}, '', ''
disp('reading reference %s ...' % options.reffile)
for line in open(options.reffile):
    if line[0] == '>': 
        if any(cr): 
            if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

if len(options.chroms) == 0 or cr in options.chroms: ref[cr] = seq.upper()
del seq

meth, depth, coverage = {}, {}, {}
for cr in ref:
    meth[cr] = array.array('I', [0]) * len(ref[cr])
    depth[cr] = array.array('I', [0]) * len(ref[cr])
    if options.rm_dup: coverage[cr] = array.array('B', [0]) * len(ref[cr])

options.chroms = set(ref.keys())

BS_conversion = {'+': ('C','T'), '-': ('G','A')}
nmap = 0
for infile in infiles:
    nline = 0
    disp('reading %s ...' % infile)
    if infile[-4:].upper() == '.SAM': sam_format, fin = True, os.popen('%ssamtools view -XS %s' % (options.sam_path, infile)) 
    elif infile[-4:].upper() == '.BAM': sam_format, fin = True, os.popen('%ssamtools view -X %s' % (options.sam_path, infile))
    else: sam_format, fin = False, open(infile)
    for line in fin:
        nline += 1
        if nline % 10000000 == 0: disp('read %d lines' % nline, nt=1)
        map_info = get_alignment(line)
        if len(map_info) == 0: continue
        seq, strand, cr, pos = map_info
        depthcr = depth[cr]
        if pos + len(seq) > len(depthcr): continue
        nmap += 1
        methcr = meth[cr]
        refseq = ref[cr][pos:pos+len(seq)]
        match, convert = BS_conversion[strand]
        index = refseq.find(match)
        while index >= 0:
            if seq[index] == convert: depthcr[pos+index] += 1
            elif seq[index] == match: 
                methcr[pos+index] += 1
                depthcr[pos+index] += 1
            index = refseq.find(match, index+1)    
    
    fin.close()

if options.combine_CpG:
    disp('combining CpG methylation from both strands ...')
    for cr in depth:
        methcr, depthcr, refcr = depth[cr], meth[cr], ref[cr]
        pos = refcr.find('CG')
        while pos >= 0:
            depthcr[pos] += depthcr[pos+1]
            methcr[pos] += methcr[pos+1]
            depthcr[pos+1] = 0
            methcr[pos+1] = 0
            pos = refcr.find('CG', pos+2)

disp('writing %s ...' % options.outfile)
ss = {'C': '+', 'G': '-'}
fout = open(options.outfile, 'w')
z95, z95sq = 1.96, 1.96 * 1.96
fout.write('chr\tpos\tstrand\tcontext\tratio\ttotal_C\tmethy_C\tCI_lower\tCI_upper\n')
nc, nd, dep0 = 0, 0, options.min_depth
for cr in sorted(depth.keys()):
    depthcr, methcr, refcr = depth[cr], meth[cr], ref[cr]
    for i, d in enumerate(depthcr):
        if d < dep0: continue
        nc += 1
        nd += d
        m = methcr[i]
        if m == 0 and not options.meth0: continue
        ratio = float(m) / d
        seq = refcr[i-2:i+3]
        strand = ss[refcr[i]]
        pmid = ratio + z95sq / (2 * d)
        sd = z95 * ((ratio*(1-ratio)/d + z95sq/(4*d*d)) ** 0.5)
        norminator = 1 + z95sq / d
        CIl, CIu = (pmid - sd) / norminator, (pmid + sd) / norminator
        fout.write('%s\t%d\t%c\t%s\t%.3f\t%d\t%d\t%.3f\t%.3f\n' % (cr, i+1, strand, seq, ratio, d, m, CIl, CIu))

fout.close()
disp('done.')
print 'total %d valid mappings, %d covered cytosines, average coverage: %.2f fold.' % (nmap, nc, float(nd)/nc)
