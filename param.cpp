#include <unistd.h>
#include "param.h"
#include<iostream>

using namespace std;

Param::Param()
{
	num_procs=sysconf(_SC_NPROCESSORS_ONLN);
	if(num_procs>8) num_procs=8;


/*	
#ifdef DB_CHR  // seqs <256, length <4Gb
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
#endif
#ifdef DB_CONTIG // seqs <65K, length <4Gb
	max_dbseq_size=0x100000; //1Mb
	append_dbseq_size=0x100000;  //1Mb
#endif
#ifdef DB_SHORT // seqs <4G, length <65K
	max_dbseq_size=0x10000; //65Kb
	append_dbseq_size=0x10000;  //65Kb
#endif
#ifdef DB_HUGE // seqs <4G, length <4G
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
#endif
*/
	max_dbseq_size=0x1000000; //16Mb
	append_dbseq_size=0x1000000;  //16Mb
	
//	read_size=30;
	max_ns = 5;
	trim_lowQ=0;
	
	zero_qual= '!';
	qual_threshold= 0;
	default_qual=40;
	
	min_insert= 28;
	max_insert= 500;
	optimize_output_SV=1;
	
	seed_size= 16;
	half_seed_size= seed_size>>1;	
	half_seed_bits= (1<<(half_seed_size*2))-1;
	seed_bits=(1<<(seed_size*2))-1;
	
	max_snp_num = 2;
	max_num_hits = MAXHITS;
	
	min_read_size=seed_size;
	
    n_adapter=0;
		
	report_repeat_hits = 1;
	output_id=1;

	useful_nt="ACGTacgt";
	nx_nt="NXnx";
	
	for(int i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;

	BuildMismatchTable();
	
	CCGG_max=1000000000;
	CCGG_min=0;
	out_sam=0;
	read_start=1;
	read_end=~0;

    //SetDigestionSite("C-CGG");
    out_ref=0;
    out_unmap=0;
    RRBS_flag=0;
    index_interval=4;
    randseed=0;
    chains=0;
    pairend=0;
    max_readlen=(FIXELEMENT-1)*16;
    SetAlign('T','C');

};

void Param::InitMapping(){
    SeedProfile e;
	for(int i=0; i<index_interval; i++) for(int j=0; j<=MAXSNPS; j++) {  //for 4 binary seqs
		e.a=((j*seed_size+i+index_interval-1)/index_interval)*index_interval; //[half_seed_size*2+i, half_seed_size*2+i+4), 4n
		e.b1=e.a/SEGLEN;
		e.s1=SEGLEN*4-(seed_size+e.a%SEGLEN)*2;		//shift s1 to right for 2 elements
		profile[j][i]=e;
	}
}

void Param::SetDigestionSite(const char *a) {
    digest_site=a;
    //cout<<digest_site<<endl;
    if((digest_pos=digest_site.find('-'))<0) {
        cout<<"Digestion position not marked, use \'-\' to mark. example: \'C-CGG\'\n";  
        exit(1);
    }
    digest_site.erase(digest_pos,1);
    RRBS_flag=1;
    index_interval=1;
    SetSeedSize(12);
}

void Param::SetSeedSize(int n)
{
	//seed_size=(n/2)*2;
	seed_size=n;
	half_seed_size= seed_size>>1;	
	half_seed_bits= (1<<(half_seed_size*2))-1;
	//seed_bits=(1<<(seed_size*2))-1;
	min_read_size=seed_size;	
	
	seed_bits=0;
	for(int i=0; i<seed_size; i++) seed_bits|=0x3<<i*2;
}


void Param::BuildMismatchTable()
{
	bit32_t TT, i;
    int n, j;
	//num_mismatch=new bit8_t[0x10000];
	//_C=new bit16_t[0x10000];
	_T=new bit16_t[0x10000];	
	
	for(i=0; i<0x10000; i++) {
	    TT = ((~((i<<1)&i))|0x5555)&i;
	    n=0;
	    for(j=7; j>=0; j--) n=n*3+((TT>>j*2)&0x3);
	    _T[i] = n;
	}

}

bit8_t alphabet[256];
bit8_t rev_alphabet[256];
bit8_t alphabet0[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
0,0,1,0,0,0,2,0,0,0,0,0,0,0,0, /* next is 'p' */
0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

bit8_t reg_alphabet[256]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'A' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0, /* next is 'a' */ 
3,0,3,0,0,0,3,0,0,0,0,0,0,0,0, /* next is 'p' */
0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

bit8_t nv3='N'; //set any unknown char as 'N'
char rev_char[256] = {
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'A' */ 
'T',nv3,'G',nv3,nv3,nv3,'C',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'P' */
nv3,nv3,nv3,nv3,'A',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'a' */ 
't',nv3,'g',nv3,nv3,nv3,'c',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3, /* next is 'p' */
nv3,nv3,nv3,nv3,'a',nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,
nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3,nv3
};

bit8_t bit_nt[4];

char nt_code[4] =
{
	'A', 'C', 'G', 'T'
};


void Param::SetAlign(char readnt, char refnt){
    int i, j; bit8_t tmp=0;
    read_nt=toupper(readnt); ref_nt=toupper(refnt);
    if((!reg_alphabet[read_nt])||(!reg_alphabet[ref_nt])) {
        cerr<<"Unknown nucleotide."<<endl;
        exit(1);
    }
    if(read_nt==ref_nt){
        cerr<<"Must specify different nucleotides for additional alignment."<<endl;
        exit(1);
    }

    for(i=0;i<4;i++) bit_nt[i]=100;    
    bit_nt[alphabet0[read_nt]]=3;
    bit_nt[alphabet0[ref_nt]]=1;
    for(i=0;i<4;i++)
        if(nt_code[i]!=ref_nt&&nt_code[i]!=read_nt) {
            for(j=0;j<4;j++) if(bit_nt[j]==100) break;
            bit_nt[i]=tmp; tmp=2;
        }

    //for(i=0;i<4;i++) cout<<" "<<nt_code[i]<<":"<<(int) bit_nt[i]; cout<<endl;

    for(i=0;i<256;i++) alphabet[i]=bit_nt[0];
    alphabet['c']=bit_nt[1]; alphabet['C']=bit_nt[1]; 
    alphabet['g']=bit_nt[2]; alphabet['G']=bit_nt[2];
    alphabet['t']=bit_nt[3]; alphabet['T']=bit_nt[3];

    for(i=0;i<256;i++) rev_alphabet[i]=bit_nt[3];
    rev_alphabet['c']=bit_nt[2]; rev_alphabet['C']=bit_nt[2]; 
    rev_alphabet['g']=bit_nt[1]; rev_alphabet['G']=bit_nt[1];
    rev_alphabet['t']=bit_nt[0]; rev_alphabet['T']=bit_nt[0];

    for(i=0;i<4;i++) useful_nt[bit_nt[i]]=nt_code[i];
    for(i=0;i<4;i++) useful_nt[bit_nt[i]+4]=tolower(nt_code[i]);

    //cout<<useful_nt<<endl;

    //for(i=0;i<4;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<4;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(alphabet[nt_code[i]]); cout<<endl;

    //for(i=0;i<4;i++) cout<<" "<<(char) nt_code[i]<<":"<<(int) rev_alphabet[nt_code[i]]; cout<<endl;
    //for(i=0;i<4;i++) cout<<" "<<(char) tolower(nt_code[i])<<":"<<(int) tolower(rev_alphabet[nt_code[i]]); cout<<endl;

}


char chain_flag[2] =
{
'+', '-'
};


char revnt_code[4] =
{
	'T', 'G', 'C', 'A'
};


void Param::rc_seq(string &seq){
   	reverse(seq.begin(), seq.end());
   	for(int ii=0;ii<seq.size();ii++) seq[ii]=rev_char[seq[ii]];
}


