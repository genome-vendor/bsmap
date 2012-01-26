#include <stdio.h>
#include "sam.h"
/* callback for bam_fetch() */
static int fetch_func(const bam1_t *b)
{
    const bam1_core_t *c = &b->core;
    int i;
    char* read_name=(char*) bam1_qname(b);
    printf("%s\t",read_name);
    char* read_seq=(char*)malloc(c->l_qseq+1);
    char* s=(char*) bam1_seq(b);
    for(i=0;i<c->l_qseq;i++) read_seq[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
    read_seq[i]=0;
    printf("%s\t",read_seq);
    char* read_qual=(char*)malloc(c->l_qseq+1);
    char* t=(char*) bam1_qual(b);    
    for(i=0;i<c->l_qseq;i++) read_qual[i]=t[i]+33;
    read_qual[i]=0;
    printf("%s\n",read_qual);
    free(read_seq); free(read_qual);
    return 0;
}
int main(int argc, char *argv[])
{
    samfile_t *fp;
    if ((fp = samopen(argv[1], "rb", 0)) == 0) {
        fprintf(stderr, "showbam: Fail to open BAM file %s\n", argv[1]);
        return 1;
    }
    bam1_t *b = bam_init1();
    while (samread(fp, b) >= 0) fetch_func(b);
    bam_destroy1(b);
    samclose(fp);
    return 0;
}
