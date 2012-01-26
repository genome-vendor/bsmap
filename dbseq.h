#ifndef _DBSEQ_H_
#define _DBSEQ_H_

#include<vector>
#include<set>
#include<string>
#include<iostream>
#include<fstream>
#include "param.h"
#include "utilities.h"

#define PREFETCH_CAL_INDEX 16
#define PREFETCH_CRT_INDEX 8
#define PREFETCH_LOOP 8U
#define REF_MARGIN 400

using namespace std;

struct OneBfa
{
	bit32_t n;  //count
	bit24_t *s;
};
struct RefTitle
{
	string name;
	bit32_t size;
	//added by yxi
	bit32_t rc_offset;
};
struct Block
{
	bit32_t id;
	bit32_t begin;
	bit32_t end;
};
struct KmerLoc
{
	bit32_t n1; //ab, ac, ad seed
	Hit *loc1;
};

struct shortHit{
    ref_id_t chr;
    ref_loc_t loc;
};

struct shorthitcompclass{
    bool operator()(shortHit a, shortHit b) {
    	if(a.loc<b.loc) return 1;
    	else if(a.loc>b.loc) return 0;
    	else if(a.chr<b.chr) return 1;
    	else return 0;
    }
};



class RefSeq
{
public:
	RefSeq();
	ref_loc_t LoadNextSeq(ifstream &fin);
	void BinSeq(OneBfa &a);
	void cBinSeq(OneBfa &a);
	void UnmaskRegion();
	void Run_ConvertBinseq(ifstream &fin);
	inline bit32_t s_MakeSeed_1(bit24_t *_m, int _a);
	inline bit32_t s_MakeSeed_2(bit24_t *_m, int _a);
	//inline bit32_t s_MakeSeed_2(bit24_t *_m, int _a, int _b, int _c, int _d);
	
	void InitialIndex();
	void t_CalKmerFreq_ab();
	//void t_CalKmerFreq_ac();
	//void t_CalKmerFreq_ad();
	void AllocIndex();
	void t_CreateIndex_ab();
	//void t_CreateIndex_ac();
	//void t_CreateIndex_ad();
	void CreateIndex();
	void ReleaseIndex();
	void find_CCGG();
    pair<ref_loc_t,int> CCGG_seglen(ref_id_t chr, ref_loc_t pos, int readlen);
    ref_loc_t hit2int(Hit h);
    Hit int2hit(ref_loc_t p, int c);

public:
	int total_num;
	bit64_t sum_length;
	vector<OneBfa> bfa;
	bit32_t total_kmers;
	KmerLoc *index;
	vector<RefTitle> title;	
protected:
	ref_id_t _count;
	string _name;
	string _seq;
	ref_loc_t _length;
	shortHit tmploc;
public:	
	vector<Block> _blocks;  //unmasked ref region
	//map<shortHit,bit32_t,shorthitcompclass> ccgg_seglen;

    //by yxi
    //int max_seedseg_num;
	vector<vector<ref_loc_t> > CCGG_index[50];
	vector<vector<ref_loc_t> > CCGG_sites;
    vector<ref_loc_t> *CCGG_sites_chr;
    int n_CCGG;
    NewIndex *index2;
    bit32_t *index2_count;
    bit32_t *refcat, *crefcat;
    vector <bit32_t> ref_anchor, cref_anchor;
};

#endif //_DBSEQ_H_
