#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<vector>
#include<string>
#include<cstdio>
#include<set>
#include<algorithm>
#include<cmath>

#include "param.h"
#include "reads.h"
#include "dbseq.h"

using namespace std;

const int FIXSIZE=SEGLEN*FIXELEMENT; 
typedef Hit HitArray[MAXHITS+1];

extern Param param;
extern char rev_char[];
extern char chain_flag[];

class SingleAlign
{
public:
	SingleAlign();
	~SingleAlign();
	void ImportFileFormat(int format);
	void ImportBatchReads(bit32_t n, vector<ReadInf> &a);
	int CountNs();
	void set_RRBS_start(); //by yxi
	int TrimLowQual();
	void ConvertBinaySeq();
	inline void GenerateSeeds(int n, int start);
	inline void GenerateCSeeds(int n, int cstart);
    int CountSeeds(RefSeq &ref, int n, int start);
    int CountCSeeds(RefSeq &ref, int n, int cstart);
	
	//inline void GenerateSeeds_2(int n);
	//inline void GenerateSeeds_3(int n);
	inline unsigned int CountMismatch(bit64_t *q, bit64_t *r, bit64_t *s);
	//void SnpAlign_0(RefSeq &ref);
	void SortHits(int n);
	//void SortcHits(int n);
	//void SnpAlign_1(RefSeq &ref);
	inline bool equal_loc(Hit a);
	//void SnpAlign_2(RefSeq &ref);
	void ClearHits();
	int RunAlign(RefSeq &ref);
	int FilterReads();
	//int CountStringMismatch(int offset, string &s1, string s2);
	void Do_Batch(RefSeq &ref);
	void StringAlign(RefSeq &ref, string &os);
	void Reverse_Seq();
	void Reverse_Qual();
	void s_OutHit(int chain, int n, bit8_t nspsn, Hit *hit, int insert_size, RefSeq &ref, string &os);
	//void s_OutGapHit(int chain, size_t n, bit8_t g, Hit *hit, RefSeq &ref, string &os) {};
    
    //by yxi
    void SnpAlign(RefSeq &ref, int mode);
    int TrimAdapter();
    void SortHits4PE(int n);
    void Fix_Unpaired_Short_Fragment(RefSeq &ref);
    void ReorderSeed(RefSeq &ref);
    bit32_t GetTotalSeedLoc(RefSeq &ref, int start);
    bit32_t GetTotalCSeedLoc(RefSeq &ref, int start);
    void AdjustSeedStartArray(RefSeq &ref);
    void AdjustCSeedStartArray(RefSeq &ref);

public:

    //Hot data section
	bit24_t bseq[SEGLEN][10];
	bit24_t reg[SEGLEN][10];
	bit24_t cbseq[SEGLEN][10];
	bit24_t creg[SEGLEN][10];

	bit32_t seeds[MAXSNPS+1][16];
	bit32_t cseeds[MAXSNPS+1][16];
	
	bit32_t seed_array[144];
	bit32_t cseed_array[144];

	int _cur_n_hit[MAXSNPS+1];
	int _cur_n_chit[MAXSNPS+1];

	bit32_t snp_thres;
	bit32_t rand_rSeed; // thread safe RNG seed
	int seed_start_offset, cseed_start_offset;
	int seed_start_array[MAXSNPS+1], cseed_start_array[MAXSNPS+1];
	bit32_t cseed_offset;
    set<ref_loc_t> *hitset; //, *chitset; 
    vector<pair<int,int> > seedindex, cseedindex;
	SeedProfile *_pro;
	bit32_t _seed;
	Hit _hit;
	Hit *_refloc;
    ref_loc_t *_refloc2, *_refchr2;
	bit32_t tmp_snp;
	bit32_t _hitz;
    bit32_t _ref_chr_count;  	
	
	//Hit hits[MAXSNPS+1][MAXHITS+1];
	//Hit chits[MAXSNPS+1][MAXHITS+1];
	HitArray *hits, *chits;

    //cold data section
	int _format;
    bool flag_chain, cflag_chain;	
	vector<ReadInf>::iterator _pread;
	string _ori_read_seq;
	string _ori_read_qual;
	string _revseq;
	string _revqual;
	bit32_t num_reads;
	vector<ReadInf> mreads;
	bit32_t n_aligned;
	string _str_align;   //align results, prepare for output
    int raw_readlen;
    int read_max_snp_num;
    vector<string>::iterator read_motif_iter;
    int seedseg_num;

	//local variables
	string::iterator _sp;
	string::reverse_iterator _sq;
	string _str;
	
	char _ch[1024];
    pair<ref_loc_t,int> seg_info;
    char _mapseq[256];
    string::iterator _readnt, _adapternt;
};

/*n=0: ab; 1: cd; 2: bc; 3: ac; 4: bd; 5: ad*/
//n<3: ab, cd, bc	
inline void SingleAlign::GenerateSeeds(int n, int start)
{
	int i;
 	//cout<<"cseed_offset="<<cseed_offset<<endl;

    if(param.RRBS_flag){
    	_pro=param.profile[n];
    	seeds[n][0]=seed_array[_pro->a+start];
    }
    else{ 
        for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) seeds[n][i]=seed_array[_pro->a+start-i];
    }
}

inline void SingleAlign::GenerateCSeeds(int n, int cstart)
{
	int i;
 	//cout<<"cseed_offset="<<cseed_offset<<endl;

    if(param.RRBS_flag){
    	_pro=param.profile[n];
        cseeds[n][0]=cseed_array[_pro->a+cseed_offset+cstart];
    }
    else{ 
        for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) cseeds[n][i]=cseed_array[_pro->a+cstart-i];
    }
}


inline unsigned int SingleAlign::CountMismatch(register bit64_t *q, register bit64_t *r, register bit64_t *s)
{

#ifdef READ_48
    return param.XM64((*q&param.XC64(*s)^*s)&*r)+param.XM64((*(q+1)&param.XC64(*(s+1))^*(s+1))&*(r+1));

#endif
#ifdef READ_80

    if((tmp_snp=param.XM64((*q&param.XC64(*s)^*s)&*r))>snp_thres)
        return tmp_snp;

    return tmp_snp+param.XM64((*(q+1)&param.XC64(*(s+1))^*(s+1))&*(r+1))
        +param.XM64((*(q+2)&param.XC64(*(s+2))^*(s+2))&*(r+2));

/*
    return tmp_snp+param.XM64X2((*((bit64_t*)q)&param.XC64(*((bit64_t*)s))^*((bit64_t*)s))&*((bit64_t*)r), 
    (*((bit64_t*)q+2)&param.XC64(*((bit64_t*)s+2))^*((bit64_t*)s+2))&*((bit64_t*)r+2));
*/
#endif
#ifdef READ_144

    if((tmp_snp=param.XM64((*q&param.XC64(*s)^*s)&*r))>snp_thres)
        return tmp_snp;

    if((tmp_snp+=param.XM64((*(q+1)&param.XC64(*(s+1))^*(s+1))&*(r+1)))>snp_thres)
        return tmp_snp;

    return tmp_snp+param.XM64((*(q+2)&param.XC64(*(s+2))^*(s+2))&*(r+2))
           +param.XM64((*(q+3)&param.XC64(*(s+3))^*(s+3))&*(r+3))
           +param.XM64((*(q+4)&param.XC64(*(s+4))^*(s+4))&*(r+4));

#endif
}	

inline void SingleAlign::Reverse_Seq()
{
	_revseq=_pread->seq;
	reverse(_revseq.begin(), _revseq.end());
	for(string::iterator p=_revseq.begin(); p!=_revseq.end(); ++p)
		*p=rev_char[*p];
}

inline void SingleAlign::Reverse_Qual()
{
	_revqual=_pread->qual;
	reverse(_revqual.begin(), _revqual.end());
}
#endif //_ALIGN_H_
