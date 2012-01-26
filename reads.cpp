#include "reads.h"

using namespace std;

extern Param param;

ReadClass::ReadClass()
{
	_index=0;
	mreads.resize(BatchNum);
}

void ReadClass::CheckFile(ifstream &fin, string filename, int mode)
{
	string s1,s2,s3,s4;
	char ch[1000];
	bit32_t i;

	fin>>s1;
	fin.getline(ch, 1000);

	if('>' == s1[0]) {
		_file_format=1;
    	fin>>s2;
    	fin.getline(ch, 1000);
    }
	else if('@' == s1[0]) {
    	fin>>s2;
    	fin.getline(ch, 1000);
		fin>>s3;		
		fin.getline(ch, 1000);
		fin>>s4;
		fin.getline(ch, 1000);
		_file_format=0;
		if(s2.size() != s4.size()) {
			cerr<<"fatal error: fq format, sequence length not equal to quality length\n";
			exit(1);
		}		
	}
	else if((SAM_fp=samopen(filename.c_str(), "rb", 0))!=0) {
    	SAM_b=bam_init1();
	    _file_format=2; //BAM format
	}
	else {
		cerr<<"fatal error: unrecognizable format of reads file.\n";
		exit(1);
	}
	fin.seekg(0);

	switch(_file_format) {
	case 0: //fastq
		for(i=0;i<(param.read_start-1)*4;i++) {
			if(fin.eof()) break;
			fin.getline(ch,1000);
		}
		break;
	case 1: //fasta
		for(i=0;i<(param.read_start-1)*2;i++) {
			if(fin.eof()) break;
			fin.getline(ch,1000);
		}
		break;
	case 2: //bam
		if(mode &1){ //single-end
			for(i=0;i<param.read_start-1;i++) {if(samread(SAM_fp,SAM_b)<0) break;}
		}
		else {//pair-end
			for(i=0;i<(param.read_start-1)*2;i++) {if(samread(SAM_fp,SAM_b)<0) break;}
		}
		break;
	}
}

void ReadClass::InitialIndex()
{
	_index=param.read_start-1;
}

int ReadClass::LoadBatchReads(ifstream &fin, int readset)  // readset {0: single-end, 1:pair-end set1, 2:pair-end set2)  
{
	char ch[1000];
	char c;
	vector<ReadInf>::iterator p=mreads.begin();

	size_t i,l_seq;
	char *s, *t;
	
	if (_file_format<2) //.fa and .fq format
    	for(num=0; num<BatchNum; p++,num++,_index++){
		if(_index>=param.read_end) break;
    		fin>>c;
    		if(fin.eof()) break;
    		p->index=_index;
	    	fin>>p->name;
	    	fin.getline(ch,1000);
	    	fin>>p->seq;
	    	p->readset=readset;
	    	if(!_file_format) {//*.fq
	    		fin>>ch;
	    		fin.getline(ch, 1000);
	    		fin>>p->qual;
	    	}
	    	else p->qual=string(p->seq.size(), param.zero_qual+param.default_qual);
	    	/*
	    	cout<<p->qual<<endl;
	    	if(param.out_sam&&(param.zero_qual!='!')) 
	    	    for(it=p->qual.begin(); it!=p->qual.end();++it) *it-=(param.zero_qual-'!');
	    	cout<<p->qual<<endl;
	    	*/
            if(p->seq.size()>param.max_readlen) {
                p->seq.erase(param.max_readlen); p->qual.erase(param.max_readlen);
            }
	    }
	else //BAM format
    	for(num=0; num<BatchNum; num++,p++,_index++){
		if(_index>=param.read_end) break;
		//cout<< "index:"<<_index<<endl;
            //cout<<"num:"<<num<<"  mode:"<<mode<<endl;
            if(readset==2) if (samread(SAM_fp,SAM_b)<0) break;
    		if (samread(SAM_fp,SAM_b)<0) break;
      		p->index=_index;
       		p->name=string((char*)bam1_qname(SAM_b));
           //  
       		l_seq=min(SAM_b->core.l_qseq,param.max_readlen);
       		p->seq.assign(l_seq,0); p->qual.assign(l_seq,0);
       		s=(char*) bam1_seq(SAM_b);  t=(char*) bam1_qual(SAM_b);
            if(readset) {
                if(SAM_b->core.flag&0x40) p->readset=1;
                else if(SAM_b->core.flag&0x80) p->readset=2;
                else p->readset=readset;
            }
            else p->readset=readset;
       		for(i=0;i<l_seq;i++){
       		    p->seq[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
       		    p->qual[i]=t[i]+33;
       		}
			//cout<<p->name<<" "<<p->seq<<endl;
            if(readset==1) if (samread(SAM_fp,SAM_b)<0) break;
		}

	return num;
}

