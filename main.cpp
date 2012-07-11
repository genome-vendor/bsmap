#include<unistd.h>
#include<cstdio>
#include<iostream>
#include<iomanip>
#include<ostream>
#include<fstream>
#include<string>
#include<vector>
#include "reads.h"
#include "dbseq.h"
#include "align.h"
#include "param.h"
#include "pairs.h"
#include "utilities.h"

#ifdef THREAD
#include<pthread.h>
#endif

using namespace std;

//global variables
Param param;
string query_a_file;
string query_b_file;
string ref_file;
string out_align_file;
string out_align_file_unpair;

ifstream fin_db;
ifstream fin_a;
ifstream fin_b;
ofstream fout;
ofstream fout_unpair;
ReadClass read_a;
ReadClass read_b;
RefSeq ref;

bit32_t n_aligned=0;   //number of reads aligned
bit32_t n_aligned_pairs=0;  //number of pairs aligned
bit32_t n_aligned_a=0;  //number of a reads aligned
bit32_t n_aligned_b=0;  //number of b reads aligned
char version[] = "2.6";

#ifdef THREAD
pthread_mutex_t mutex_fin=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_fout=PTHREAD_MUTEX_INITIALIZER;

void *t_SingleAlign(void *)
{
	SingleAlign a;
	int n;
	bit32_t cur_at;
	a.ImportFileFormat(read_a._file_format);
	while(1)
	{
		pthread_mutex_lock(&mutex_fin);
		n=read_a.LoadBatchReads(fin_a,0);
		cur_at=read_a._index;
		a.ImportBatchReads(read_a.num, read_a.mreads);
		pthread_mutex_unlock(&mutex_fin);
		if(!n)
			break;
		a.Do_Batch(ref);
		pthread_mutex_lock(&mutex_fout);
		fout<<a._str_align;
		cout<<cur_at-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl;
		pthread_mutex_unlock(&mutex_fout);		
	}
	pthread_mutex_lock(&mutex_fout);
	n_aligned+=a.n_aligned;
	pthread_mutex_unlock(&mutex_fout);
};
void Do_SingleAlign()
{
	read_a.CheckFile(fin_a, query_a_file.c_str(), 1);
	vector<pthread_t> pthread_ids(param.num_procs);
	//create
	for(int i=0; i<param.num_procs; i++)
		pthread_create(&pthread_ids[i], NULL, t_SingleAlign, NULL);
	//join
	for (int i=0; i<param.num_procs; i++)
		pthread_join(pthread_ids[i], NULL);
};


void *t_PairAlign(void *)
{
	PairAlign a;
	int n1, n2;
	bit32_t cur_at;
	a.ImportFileFormat(read_a._file_format, read_b._file_format);
	while(1) {
		pthread_mutex_lock(&mutex_fin);
		n1=read_a.LoadBatchReads(fin_a,1);
		n2=read_b.LoadBatchReads(fin_b,2);
		cur_at=read_a._index;
		a.ImportBatchReads(n1, read_a.mreads, read_b.mreads);
		pthread_mutex_unlock(&mutex_fin);
		if(!n1||(n1!=n2))
			break;
		a.Do_Batch(ref);
		pthread_mutex_lock(&mutex_fout);
		fout<<a._str_align;
		fout_unpair<<a._str_align_unpair;
		cout<<cur_at-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl;
		pthread_mutex_unlock(&mutex_fout);		
	}
	pthread_mutex_lock(&mutex_fout);
	n_aligned_pairs+=a.n_aligned_pairs;
	n_aligned_a+=a.n_aligned_a;
	n_aligned_b+=a.n_aligned_b;
	pthread_mutex_unlock(&mutex_fout);		
};

void Do_PairAlign()
{
	//if(param.max_snp_num>0) param.max_snp_num=2;
	read_a.CheckFile(fin_a, query_a_file.c_str(), 2);
	read_b.CheckFile(fin_b, query_b_file.c_str(), 2);
	
	vector<pthread_t> pthread_ids(param.num_procs);
	//create
	//cout <<param.num_procs<<"num_procs\n";
	for(int i=0; i<param.num_procs; i++)
		pthread_create(&pthread_ids[i], NULL, t_PairAlign, NULL);
	//join
	for (int i=0; i<param.num_procs; i++)
		pthread_join(pthread_ids[i], NULL);
	//
};

#else
void Do_SingleAlign()
{
	read_a.CheckFile(fin_a, query_a_file.c_str(), 1);
	SingleAlign a;
	a.ImportFileFormat(read_a._file_format);
	while(read_a.LoadBatchReads(fin_a,0))
	{
		a.ImportBatchReads(read_a.num, read_a.mreads);
		a.Do_Batch(ref);
		fout<<a._str_align;
		cout<<read_a._index-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl;
	}
	n_aligned=a.n_aligned;	
};

void Do_PairAlign()
{
	read_a.CheckFile(fin_a, query_a_file.c_str(), 2);
	read_b.CheckFile(fin_b, query_b_file.c_str(), 2);
	PairAlign a;
	int n1, n2;
	while(1)
	{
		n1=read_a.LoadBatchReads(fin_a,1);
		n2=read_b.LoadBatchReads(fin_b,2);
		if(!n1||(n1!=n2))
			break;
		a.ImportBatchReads(n1, read_a.mreads, read_b.mreads);
		a.Do_Batch(ref);		
		fout<<a._str_align;
		fout_unpair<<a._str_align_unpair;
		cout<<read_a._index-param.read_start+1<<" reads finished. "<<Cal_AllTime()<<" secs passed"<<endl;
	}	
	n_aligned_pairs=a.n_aligned_pairs;
	n_aligned_a=a.n_aligned_a;
	n_aligned_b=a.n_aligned_b;
};

#endif

void Do_Formatdb()
{
	ref.CreateIndex();
	cout<<"Create seed table. "<<Cal_AllTime()<<" secs passed\n";
};


//usage
void usage(void)
{
cout<<"Usage:	bsmap [options]\n"
		<<"       -a  <str>   query a file, FASTA/FASTQ/BAM format\n"
		<<"       -d  <str>   reference sequences file, FASTA format\n"
		<<"       -o  <str>   output alignment file, BSP/SAM/BAM format\n"
		<<"\n  Options for alignment:\n"
		<<"       -s  <int>   seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.\n"
		<<"       -v  <int>   maximum number of mismatches allowed on a read, <="<<MAXSNPS<<". default="<<param.max_snp_num<<".\n"
		<<"       -w  <int>   maximum number of equal best hits to count, <="<<MAXHITS<<"\n"
        <<"       -B  <int>   start from the Nth read or read pair, default: 1\n"
        <<"       -E  <int>   end at the Nth read or read pair, default: 4,294,967,295\n"
       	<<"       -I  <int>   index interval, default="<<param.index_interval<<"\n"        
#ifdef THREAD
		<<"       -p  <int>   number of processors to use, default="<<param.num_procs<<"\n"
#endif	
        <<"       -D  <str>   activating RRBS mapping mode and set restriction enzyme digestion sites. \n"
        <<"                   digestion position marked by \'-\', example: -D C-CGG for MspI digestion.\n"
        <<"                   default: none (whole genome shotgun bisulfite mapping mode)\n"        
    	<<"       -S  <int>   seed for random number generation used in selecting multiple hits\n"
        <<"                   other seed values generate pseudo random number based on read index number, to allow reproducible mapping results. \n"
        <<"                   default="<<param.randseed<<". (get seed from system clock, mapping results not resproducible.)\n"
        <<"       -n  [0,1]   set mapping strand information. default: -n "<<param.chains<<"\n"
        <<"                   -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+), \n" 
        <<"                   for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.\n"
        <<"                   -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, -- \n"
        <<"       -M  <str>   set alignment information for the additional nucleotide transition. \n"
        <<"                   <str> is in the form of two different nucleotides N1N2, \n"
        <<"                   indicating N1 in the reads could be mapped to N2 in the reference sequences.\n"
        <<"                   default: -M TC, corresponds to C=>U(T) transition in bisulfite conversion. \n"
        <<"                   example: -M GA could be used to detect A=>I(G) transition in RNA editing. \n"
		<<"\n  Options for trimming:\n"
		<<"       -q  <int>   quality threshold in trimming, 0-40, default=0 (no trim)\n"
		<<"       -z  <int>   base quality, default="<<(int) param.zero_qual<<" [Illumina is using 64, Sanger Institute is using 33]\n"
		<<"       -f  <int>   filter low-quality reads containing >n Ns, default="<<param.max_ns<<"\n"
        <<"       -A  <str>   3-end adapter sequence, default: none (no trim)\n"
        <<"       -L  <int>   map the first N nucleotides of the read, default:"<<param.max_readlen<<" (map the whole read).\n"
		<<"\n  Options for reporting:\n"
        <<"       -r  [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default:"<<param.report_repeat_hits<<".\n"
        <<"       -R          print corresponding reference sequences in SAM output, default=off\n"
        <<"       -u          report unmapped reads, default=off\n"
		<<"\n  Options for pair-end alignment:\n"
		<<"       -b  <str>   query b file\n"
		<<"       -m  <int>   minimal insert size allowed, default="<<param.min_insert<<"\n"
		<<"       -x  <int>   maximal insert size allowed, default="<<param.max_insert<<"\n"
		<<"       -2  <str>   output file of unpaired alignment hits\n"

		<<"       -h          help\n\n"
		<<endl;
	exit(1);
};

int mGetOptions(int rgc, char *rgv[])
{
	//[options]
	int i;
	for(i=1; i<rgc; i++) {
		if(rgv[i][0]!='-') return i;		
		switch(rgv[i][1]) {
			case 'a': if(rgv[i][2]==0) query_a_file = rgv[++i]; else if(rgv[i][2]=='=') query_a_file=rgv[i]+3; else return i; break;
			case 'b': if(rgv[i][2]==0) query_b_file = rgv[++i]; else if(rgv[i][2]=='=') query_b_file=rgv[i]+3; else return i; 
			    param.pairend=1; break;
			case 'd': if(rgv[i][2]==0) ref_file = rgv[++i]; else if(rgv[i][2]=='=') ref_file=rgv[i]+3; else return i; break;
			case 's': if(rgv[i][2]==0) 
				param.SetSeedSize(atoi(rgv[++i])); else if(rgv[i][2]=='=') param.SetSeedSize(atoi(rgv[i]+3)); else return i;
				if(param.RRBS_flag) param.SetSeedSize(12); 
				break;
			case 'o': if(rgv[i][2]==0) out_align_file = rgv[++i];else if(rgv[i][2]=='=') out_align_file=rgv[i]+3; else return i; break;
			case '2': if(rgv[i][2]==0) out_align_file_unpair = rgv[++i]; else if(rgv[i][2]=='=') out_align_file_unpair=rgv[i]+3; else return i;break;
			case 'm': if(rgv[i][2]==0) param.min_insert = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.min_insert=atoi(rgv[i]+3); else return i; break;
			//case 'n': if(rgv[i][2]==0) param.chains = atoi(rgv[++i])%4; else if(rgv[i][2]=='=') param.chains=atoi(rgv[i]+3)%4; else return i; break;
			case 'n': if(rgv[i][2]==0) param.chains=(atoi(rgv[++i])!=0); else if(rgv[i][2]=='=') param.chains=(atoi(rgv[i]+3)!=0); else return i; break;
			case 'x': if(rgv[i][2]==0) param.max_insert = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_insert=atoi(rgv[i]+3); else return i; break;
			case 'r': if(rgv[i][2]==0) param.report_repeat_hits = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.report_repeat_hits=atoi(rgv[i]+3); else return i; break;
			case 'I': if(rgv[i][2]==0) param.index_interval = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.index_interval=atoi(rgv[i]+3); else return i; 
			    if(param.RRBS_flag) param.index_interval=1;
			    if(param.index_interval>16) {cerr<<"index interval exceeds max value:16\n"; exit(1);}
			    break;			
			case 'v': if(rgv[i][2]==0) param.max_snp_num = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_snp_num=atoi(rgv[i]+3); else return i;
			    if(param.max_snp_num>MAXSNPS) {cerr<<"number of mismatches exceeds max value:"<<MAXSNPS<<endl; exit(1);}
			    break;
			case 'w': if(rgv[i][2]==0) param.max_num_hits = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_num_hits=atoi(rgv[i]+3); else return i;
			    if(param.max_num_hits>MAXHITS) {cerr<<"number of multi-hits exceeds max value:"<<MAXHITS<<endl; exit(1);}
			    break;
			case 'q': if(rgv[i][2]==0) param.qual_threshold = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.qual_threshold = atoi(rgv[i]+3); else return i; break;
			case 'f': if(rgv[i][2]==0) param.max_ns=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_ns=atoi(rgv[i]+3); else return i; break;
			case 'z': if(rgv[i][2]==0) param.zero_qual=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.zero_qual=atoi(rgv[i]+3); else return i; break;
			case 'p': if(rgv[i][2]==0) param.num_procs=atoi(rgv[++i]); else if(rgv[i][2]=='=') param.num_procs=atoi(rgv[i]+3); else return i; break;
			case 'A': if(rgv[i][2]==0) param.adapter[param.n_adapter++]=rgv[++i]; else if(rgv[i][2]=='=') param.adapter[param.n_adapter++]=rgv[i]+3; else return i;   
                break;
            case 'R': if(rgv[i][2]==0) param.out_ref=1; else return i; break;
            case 'u': if(rgv[i][2]==0) param.out_unmap=1; else return i; break;
		    case 'B': if(rgv[i][2]==0) param.read_start = max(atoi(rgv[++i]),1); else if(rgv[i][2]=='=') param.read_start=max(atoi(rgv[i]+3),1); else return i; break;
		    case 'E': if(rgv[i][2]==0) param.read_end = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.read_end=atoi(rgv[i]+3); else return i; break;
	        case 'D': if(rgv[i][2]==0) param.SetDigestionSite(rgv[++i]); else if(rgv[i][2]=='=') param.SetDigestionSite(rgv[i]+3); else return i; break;
	        case 'M': if(rgv[i][2]==0) {i++; param.SetAlign(rgv[i][0], rgv[i][1]);}
                      else if(rgv[i][2]=='=') param.SetAlign(rgv[i][3], rgv[i][4]); 
                      else return i; 
                      break;
            case 'L': if(rgv[i][2]==0) param.max_readlen = atoi(rgv[++i]); else if(rgv[i][2]=='=') param.max_readlen = atoi(rgv[i]+3); else return i; break;
            case 'S': if(rgv[i][2]==0) param.randseed = atoi(rgv[++i]);  else if(rgv[i][2]=='=') param.randseed = atoi(rgv[i]+3); else return i; break;	        
			case 'h':usage();   //usage information
            default: return i;
		}
	}
    param.InitMapping();
	return 0;
};

void RunProcess(void)
{
    if(out_align_file.size()>4){
	    if(out_align_file.compare(out_align_file.size()-4,4,".sam")==0) param.out_sam=1;
	    else if (out_align_file.compare(out_align_file.size()-4,4,".bam")==0) param.out_sam=2;
    }

    cout <<"max mismatches: "<<param.max_snp_num<<"\tmax multi-hits: "<<param.max_num_hits;
    cout <<"\tmax Ns: "<<param.max_ns<<"\tseed size: "<<param.seed_size<<"\tindex interval: "<<param.index_interval<<endl;
    cout<<"quality cutoff: "<<(int)param.qual_threshold<<"\tbase quality char: '"<<param.zero_qual<<"'"<<endl;
    cout<<"min fragment size:"<<param.min_insert<<"\tmax fragemt size:"<<param.max_insert<<endl;
    cout<<"start from read #"<<param.read_start<<"\tend at read #"<<param.read_end<<endl;
    cout<<"additional alignment: "<<param.read_nt<<" in reads => "<<param.ref_nt<<" in reference"<<endl;
    cout<<"param.chains:"<<param.chains<<endl;
    
    if(param.pairend){
        cout<<"mapping strand (read_1): ++,-+";
        if(param.chains) cout<< ",+-,--";
        cout<<endl;
        cout<<"mapping strand (read_2): +-,--";
        if(param.chains) cout<< ",++,-+";
        cout<<endl;    
    }
    else {
        cout<<"mapping strand: ++,-+";
        if(param.chains) cout<< ",+-,--";
        cout<<endl;
    }
    for(int i=0; i<param.n_adapter;i++) cout<<"adapter sequence"<<i+1<<": "<<param.adapter[i]<<endl;

    if(param.RRBS_flag)
        cout <<"RRBS mode. digestion site: "<<param.digest_site.substr(0,param.digest_pos)<<'-'<<param.digest_site.substr(param.digest_pos)<<endl;

	//pair-end alignment    
	if((!query_a_file.empty()) && (!query_b_file.empty()))
	{
		cout<<"Pair-end alignment("<<param.num_procs<<" threads)\n";
		cout<<"Query: "<<query_a_file<<"  "<<query_b_file<<"  Reference: "<<ref_file<<"  Output: "<<out_align_file<<"  "<<out_align_file_unpair<<endl;
		fin_a.open(query_a_file.c_str());
		if(!fin_a) {
			cerr<<"failed to open read file #1 (check -a option): "<<query_a_file<<endl;
			exit(1);
		}
		fin_b.open(query_b_file.c_str());
		if(!fin_b) {
			cerr<<"failed to open read file #2 (check -b option): "<<query_b_file<<endl;
			exit(1);
		}		
		fout.open(out_align_file.c_str());
		if(!fout) {
			cerr<<"failed to open output file (check -o option): "<<out_align_file<<endl;
			exit(1);
		}
		if(param.out_sam){
    		char _ch[1000];
	    	fout<<"@HD\tVN:1.0\n";
	    	for(int i=0;i<ref.total_num;i++){
	    	    sprintf(_ch,"@SQ\tSN:%s\tLN:%u\n",ref.title[i<<1].name.c_str(),ref.title[i<<1].size);
	    	    fout<<_ch;
	    	}
            fout<<"@PG\tID:BSMAP_"<<version<<endl;
        }
        else{
    		fout_unpair.open(out_align_file_unpair.c_str());
	    	if(!fout_unpair) {
	    		cerr<<"failed to open output file for unpaired hits (check -2 option): "<<out_align_file_unpair<<endl;
	    		exit(1);
	    	}
	    }
		n_aligned_pairs=n_aligned_a=n_aligned_b=0;
		read_a.InitialIndex();
		read_b.InitialIndex();
		Do_PairAlign();
		fin_a.close();
		fin_b.close();
		fout.close();
		if(param.out_sam==0) fout_unpair.close();
		if(read_a._file_format==2){
		    bam_destroy1(read_a.SAM_b);
            samclose(read_a.SAM_fp);
        }
   		if(read_b._file_format==2){
		    bam_destroy1(read_b.SAM_b);
            samclose(read_b.SAM_fp);
        }

		cout<<"Total number of aligned reads: \n"
			<<"pairs:       "<<n_aligned_pairs<<" ("<<setprecision(2)<<100.0*n_aligned_pairs/(read_a._index-param.read_start+1)<<"%)\n"
			<<"single a:    "<<n_aligned_a<<" ("<<setprecision(2)<<100.0*n_aligned_a/(read_a._index-param.read_start+1)<<"%)\n"
			<<"single b:    "<<n_aligned_b<<" ("<<setprecision(2)<<100.0*n_aligned_b/(read_b._index-param.read_start+1)<<"%)\n";
	}	
	//single-read alignment
	else
	{
		if(!query_a_file.empty()) {
			fin_a.open(query_a_file.c_str());
			if(!fin_a) {
				cerr<<"failed to open read file (check -a option): "<<query_a_file<<endl;
				exit(1);
			}
		}
		else
		{
			cerr<<"missing query file(s)\n";
			exit(1);
		}
		cout<<"Single read alignment("<<param.num_procs<<" threads)\n";
		cout<<"Query: "<<query_a_file<<"  Reference: "<<ref_file<<"  Output: "<<out_align_file<<endl;
		fout.open(out_align_file.c_str());
		if(!fout) {
			cerr<<"failed to open output file (check -o option): "<<out_align_file<<endl;
			exit(1);
		}
		
		if(param.out_sam){
    		char _ch[1000];
	    	fout<<"@HD\tVN:1.0\n";
	    	for(int i=0;i<ref.total_num;i++){
	    	    sprintf(_ch,"@SQ\tSN:%s\tLN:%u\n",ref.title[i<<1].name.c_str(),ref.title[i<<1].size);
	    	    fout<<_ch;
	    	}
            fout<<"@PG\tID:BSMAP_"<<version<<endl;
	    }
		n_aligned=0;
		read_a.InitialIndex();
		Do_SingleAlign();
		fin_a.close();
		fout.close();
		if(read_a._file_format==2){
		    bam_destroy1(read_a.SAM_b);
            samclose(read_a.SAM_fp);
        }
		cout<<"Total number of aligned reads: "<<n_aligned<<" ("<<setprecision(2)
		<<100.0*n_aligned/(read_a._index-param.read_start+1)<<"%)\n";
	}

	cout<<"Done.\n";
	cout<<"Finished at "<<Curr_Time();
/*
	if(param.out_sam==2){
		char pos[1000];
		if((pos=strrchr(argv[0],'/'))) string exec_path(argv[0],(int)(pos-argv[0]));
		else string exec_path("");
		cout<<"exec_path"<<exec_path<<endl;
	}
*/

	cout<<"Total time consumed:  "<<Cal_AllTime()<<" secs\n";
};

int main(int argc, char *argv[])
{
	//print usage
    cout<<"\nBSMAP v"<<version<<endl; srand(time(NULL));
	if (argc == 1)
	{
		usage();
	}
	Initial_Time();
	cout<<"Start at:  "<<Curr_Time()<<endl;
	int noptions;
	if(noptions=mGetOptions(argc, argv)) {
        cout<<"unknown option: "<<argv[noptions]<<endl;
        exit(noptions);
    }

	fin_db.open(ref_file.c_str());
	if(!fin_db) {
		cerr<<"fatal error: failed to open ref file\n";
		exit(1);
	}
	ref.Run_ConvertBinseq(fin_db);
	cout<<"Load in "<<ref.total_num<<" db seqs, total size "<<ref.sum_length<<" bp. "<<Cal_AllTime()<<" secs passed"<<endl;			
	Do_Formatdb();
	RunProcess();
    if(param.out_sam==2){
		char sys_cmd[1000], exec_path[1000], abs_bam_file[1000];
		//realpath(argv[0], exec_path);
		realpath(out_align_file.c_str(), abs_bam_file);
		//*(strrchr(exec_path,'/')+1)=0;
		sprintf(sys_cmd,"sam2bam.sh %s",abs_bam_file);
		system(sys_cmd);
	}
    ref.ReleaseIndex();
	return 0;
}
