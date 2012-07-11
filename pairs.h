#ifndef _PAIRS_H_
#define _PAIRS_H_

#include "dbseq.h"
#include "reads.h"
#include "align.h"

using namespace std;


extern char chain_flag[];

struct PairHit
{
	bit16_t chain;
	bit8_t na, nb;   //# of snps
	int insert;
	Hit a;
	Hit b;
};

typedef PairHit PairArray[MAXHITS+1];

class PairAlign
{
public:
	PairAlign();
    ~PairAlign();
	void ImportFileFormat(int format1, int format2);
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a1, vector<ReadInf> &a2);
	int GetExactPairs();
	int GetExact2SnpPairs(RefSeq &ref);
	int GetSnp2SnpPairs(RefSeq &ref);
	int GetExact2GapPairs(RefSeq &ref);
	int RunAlign(RefSeq &ref);
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void StringAlign_ClosestUnpair(RefSeq &ref, string &os);
	
	//added by yxi
	int GetPairs(int na, int nb);
    int StringAlignPair(RefSeq &ref, string &os);
	void StringAlignUnpair(int fa, int fb, RefSeq &ref, string &os);
    void s_OutHitPair(PairHit pp, int n, RefSeq &ref, string &os);
    void s_OutHitUnpair(int readinpair, int chain_a, int chain_b, int ma, int na, Hit ha, int mb, Hit hb, RefSeq &ref, string &os);
    int TrimAdapter();
    void FixPairReadName();

public:	
	SingleAlign _sa;
	SingleAlign _sb;
	bit32_t num_reads;
	bit32_t n_aligned_pairs, n_aligned_a, n_aligned_b;	
	string _str_align;
	string _str_align_unpair;
protected:
	bit32_t _cur_n_hits[2*MAXSNPS+1];
	//PairHit pairhits[2*MAXSNPS+1][MAXHITS+1];
    PairArray *pairhits;
    bit32_t rand_rSeed;	//thread safe RNG seed
	//by yxi
    char _mapseq[256];
    char _ch[1024];  
   	SingleAlign * _stmp;  
    int checked_pair_mismatch[MAXSNPS+1][MAXSNPS+1];
};

#endif //_PAIR_ALIGH_H_
