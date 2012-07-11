#ifndef _PARAM_H_
#define _PARAM_H_

#define SEGLEN 16

#include<string>
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include<vector>
#include<algorithm>

using namespace std;

#ifdef READ_48
const int FIXELEMENT=4;
#endif

#ifdef READ_80
const int FIXELEMENT=6;
#endif

#ifdef READ_144
const int FIXELEMENT=10; // 144/16+1
#endif

const int MAXSNPS=15;

typedef unsigned char bit8_t;
typedef unsigned short bit16_t;
typedef unsigned bit32_t;
typedef unsigned long long bit64_t;
typedef bit32_t bit24_t;

typedef bit32_t ref_id_t;
typedef bit32_t ref_loc_t;
typedef ref_loc_t* NewIndex;

struct SeedProfile
{
	bit8_t a;  //offset of part a on binary seq
	bit8_t b1;  //begin element when creating seed
	bit8_t s1;  //shift when creating seed from binary seq
};

struct Hit
{
	ref_id_t chr;       //index of chr
	ref_loc_t loc;     //location of first bp on reference seq, count from 0

};


class Param
{
public:
	Param();
	void SetSeedSize(int n);
    void InitMapping();
	void BuildMismatchTable();
	void SetAdaptors(int n);
    void SetDigestionSite(const char *a);
    void rc_seq(string &seq);   
    void SetAlign(char readnt, char refnt);

public:
	int num_procs;  //number of parallel processors
	
	int chains;   //0: forward strands only ; 1: forward and reverse strands
	//dbseq
	int max_dbseq_size;
	int append_dbseq_size;
	//read
	int read_size;
	int max_ns;     //throw out reads containning >=max_ns 'N's
	int trim_lowQ;  //trim low-quality at 3'-end, or not?
	//quality
	bit8_t zero_qual;
	bit8_t qual_threshold;
	bit8_t default_qual;
	//pair-end mapping
	int min_insert;
	int max_insert;
	int optimize_output_SV;  //if a pair cannot align with proper orientation and distance, very likely a strctural variation happen here. we prefer to report hit of read 'a' and 'b' with smallest distance, so that to help detect structural variations
	//seed
	int half_seed_size;
	int seed_size;
	bit32_t half_seed_bits;
	bit32_t seed_bits;
	int min_read_size;
	//alignment
	int max_snp_num;   //maximum number of snps on one read allowed
	int max_num_hits;   //maximum number of equal best hits, smaller will be faster
    bit8_t read_nt, ref_nt;
	//report hits
	int report_repeat_hits;   //how report repeat hits? 0: no, 1: pick one randomly, 2: report all
	bool output_id;   //1: output read id, 1: out read index
	
	string adapter[10];
    int n_adapter;
    string useful_nt;
    string nx_nt;
    //added by yxi
    //bit16_t *_C;  // mask of convert T in reads to C at position of C in ref
    bit16_t *_T;  // convert T to C
    //bit16_t *map4to3; //map 3-letter sequence to 3-based number
	int CCGG_min, CCGG_max;
	int out_sam;
    int total_ref_seq;    
    int max_seedseg_num;
    bit32_t read_start, read_end;
    int out_ref;
    int out_unmap;
    string digest_site; 
    int digest_pos; 
    int RRBS_flag;
    int index_interval;
    int randseed;
	SeedProfile profile[MAXSNPS+1][16];
    int pairend;
    int max_readlen;

    inline bit32_t XT(bit32_t tt) {return (bit32_t)_T[tt&0xFFFF]+((bit32_t)_T[tt>>16])*6561UL;}; //convert T to C and map to 3-nt space
    //inline bit32_t XC(bit32_t tt) {return ((bit32_t)_C[tt&0xFFFF]|((bit32_t)_C[tt>>16])<<16);}; //mask T to C for 2 16-bit segments
    inline bit32_t XC(bit32_t tt) {return ((~tt)<<1)|tt|0x55555555U;}  // generate T2C mask according to C locations
    inline bit64_t XC64(bit64_t tt) {return ((~tt)<<1)|tt|0x5555555555555555ULL;}
    //inline bit8_t XM(bit32_t tt) {return num_mismatch[tt&0xFFFF]+num_mismatch[tt>>16];}; //count mismatches for 2 16-bit segments

    inline bit32_t XM(bit32_t tt) {
#ifdef SSE4
        return __builtin_popcount((tt|(tt>>1))&0x55555555); 
#else
        tt=(tt|(tt>>1))&0x55555555;
        tt=(tt+(tt>>2))&0x33333333;
        return (((bit32_t)(tt*0x1111111))>>28)+(tt&0x3);
#endif
    }

    inline bit32_t XM64(bit64_t tt) {
#ifdef SSE4
        return __builtin_popcountl((tt|(tt>>1))&0x5555555555555555ULL);
#else
        tt=(tt|(tt>>1))&0x5555555555555555ULL;
        tt=(tt+(tt>>2))&0x3333333333333333ULL;
        return (((tt+(tt>>4))&0x0F0F0F0F0F0F0F0FULL)*0x0101010101010101ULL)>>56;
#endif
    }

    inline bit32_t XM64X2(bit64_t tt1, bit64_t tt2) {
        tt1=((tt1|(tt1>>1))&0x5555555555555555ULL)+((tt2|(tt2>>1))&0x5555555555555555ULL);
        tt1=(tt1&0x3333333333333333ULL)+((tt1>>2)&0x3333333333333333ULL);
        return (((tt1+(tt1>>4))&0x0F0F0F0F0F0F0F0FULL)*0x0101010101010101ULL)>>56;
    }


    char * StrSeed(bit32_t seed, bit32_t size) { //for debug only 
        char *s = new char[size+1];
        for(int i=size-1; i>=0; i--) {
                s[size-1-i]=useful_nt[(seed>>(i*2))&0x3];
        }
	s[size]=0;
        return s;
    };                                        

    bit32_t map3to4(bit32_t tt){
      	int s=0, i; for(i=0;i<16;i++) {s|=(tt%3)<<i*2; tt/=3;};
	return s;
    }
};

#endif //_PARAM_H_
