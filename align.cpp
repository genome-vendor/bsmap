#include "align.h"


extern bit8_t alphabet[];
extern bit8_t reg_alphabet[];
extern bit8_t rev_alphabet[];
extern char rev_char[];
extern char nt_code[];
extern char revnt_code[];

//create seed profile
SingleAlign::SingleAlign()
{
	//
	if(param.max_snp_num >MAXSNPS) {
		cerr<<"fatal error, set smaller max_snp_num. (<=MAXSNPS)\n";
		exit(1);
	}
	n_aligned=0;

    hits = new HitArray[MAXSNPS+1];
    chits = new HitArray[MAXSNPS+1];
	
	hitset= new set<ref_loc_t>[param.total_ref_seq];
    _str_align.reserve(BatchNum*400);
	//chitset= new set<ref_loc_t>[param.total_ref_seq];
	rand_rSeed=getpid()*time(NULL);
}

SingleAlign::~SingleAlign(){
	delete [] hitset;
    delete [] hits;
    delete [] chits;

	//delete [] chitset;
}
void SingleAlign::ImportFileFormat(int format)
{
	_format=format;
}

void SingleAlign::ImportBatchReads(bit32_t n, vector<ReadInf> &a)
{
	num_reads=n;
	mreads=a;
}

int SingleAlign::CountNs()
{
	int n=0;
	for(_sp=_pread->seq.begin(); _sp!=_pread->seq.end(); _sp++)
		if(!reg_alphabet[*_sp])
			n++;
	return n;	
}


//trim low quality at 3'-end, cut at 3bp continuous high-quality bps
int SingleAlign::TrimLowQual()
{
	int i;
	if(param.qual_threshold==0||_pread->qual.size()==1) return 1;
	bit8_t read_zero_qual=param.zero_qual;
    if(param.out_sam&&(read_zero_qual!='!')){
        for(_sp=_pread->qual.begin(); _sp!=_pread->qual.end();++_sp) *_sp-=(read_zero_qual-'!');
        read_zero_qual='!';
    }

	for(_sq=_pread->qual.rbegin(),i=_pread->qual.size(); _sq!=_pread->qual.rend(); i--,++_sq) {
		if(*_sq >read_zero_qual+param.qual_threshold) {
    		if(i >= param.seed_size) {
				if(_pread->qual.size()>i) _pread->qual.erase(i);
				if(_pread->seq.size()>i) _pread->seq.erase(i);
				return 1;
			}
		}
	}
	return 0;
}

//shift 1bp to left
inline void RightShiftBinSeq(bit24_t *orib, bit24_t *newb)
{
    int i;
	newb[0]=(orib[0]>>2);
    for(i=1;i<FIXELEMENT;i++) newb[i]=(orib[i]>>2)|(orib[i-1]<<30);
}

//convert string seq to binary type
void SingleAlign::ConvertBinaySeq()
{	
	int i,j,h,g; bit32_t _a, _b, s;
    flag_chain=param.chains||(_pread->readset<2);
    cflag_chain=param.chains||(_pread->readset==2);
    if(flag_chain) {	//direct chain
    	h=_a=_b=0;
    	for(_sp=_pread->seq.begin(),i=1; _sp!=_pread->seq.end(); _sp++,i++) {
    		_a<<=2; _b<<=2;
    		_a|=alphabet[*_sp];
    		_b|=reg_alphabet[*_sp];
    		if(i>param.seed_size){
    		    s<<=2; s|=_a&0x3; 
    		    seed_array[i-param.seed_size]=param.XT(s&param.seed_bits);
    		}
    		else if(i==param.seed_size){s=_a; seed_array[0]=param.XT(s);}
    		if(0==i%SEGLEN) {
    			bseq[0][h]=_a;
    			reg[0][h++]=_b;
    			_a=_b=0;
    		}		
    	}
    	for(; i!=FIXSIZE+1; i++) {
    		_a<<=2; _b<<=2;
    		if(0==i%SEGLEN) {
    			bseq[0][h]=_a;
    			reg[0][h++]=_b;
    			_a=_b=0;
    		}		
    	}
    	for(i=1; i!=SEGLEN; i++) {
    		RightShiftBinSeq(bseq[i-1], bseq[i]);
    		RightShiftBinSeq(reg[i-1], reg[i]);
    	}
	}
/*
	for(i=0; i!=SEGLEN; i++) {
		cout<<"bin seq: "<<i<<"  "<<param.StrSeed(bseq[i][0].a,SEGLEN)<<" "<<param.StrSeed(bseq[i][1].a,SEGLEN)<<" "<<param.StrSeed(bseq[i][2].a,SEGLEN)<<" "<<param.StrSeed(bseq[i][3].a,SEGLEN)<<endl;
		cout<<"bin reg: "<<i<<"  "<<param.StrSeed(reg[i][0].a,SEGLEN)<<" "<<param.StrSeed(reg[i][1].a,SEGLEN)<<" "<<param.StrSeed(reg[i][2].a,SEGLEN)<<" "<<param.StrSeed(reg[i][3].a,SEGLEN)<<endl;
	}
*/
    if(cflag_chain) { //reverse seq
    	h=_a=_b=0;
    	for(_sq=_pread->seq.rbegin(),i=1; _sq!=_pread->seq.rend(); _sq++,i++) {
    		_a<<=2; _b<<=2;
    		_a|=rev_alphabet[*_sq];
    		_b|=reg_alphabet[*_sq];
    		if(i>param.seed_size){
    		    s<<=2; s|=_a&0x3; 
    		    cseed_array[i-param.seed_size]=param.XT(s&param.seed_bits);
    		}
    		else if(i==param.seed_size){s=_a; cseed_array[0]=param.XT(s);}
    		if(0==i%SEGLEN) {
    			cbseq[0][h]=_a;
    			creg[0][h++]=_b;
    			_a=_b=0;
    		}		
    	}
    	for(; i!=FIXSIZE+1; i++) {
    		_a<<=2; _b<<=2;
    		if(0==i%SEGLEN) {
    			cbseq[0][h]=_a;
    			creg[0][h++]=_b;
    			_a=_b=0;
    		}		
    	}
    	
    	for(i=1; i!=SEGLEN; i++) {
    		RightShiftBinSeq(cbseq[i-1], cbseq[i]);
    		RightShiftBinSeq(creg[i-1], creg[i]);
    	}
	}			
}




//by yxi
void SingleAlign::SnpAlign(RefSeq &ref, int mode)
{
   	bit32_t i,j,m,w, *ptr_chr; int h, modeindex, mc;
   	Hit prefetch_hit;

    //cout<<"READ:"<<_pread->name<<" "<<_pread->seq<<endl;

    if(param.RRBS_flag){
    	//direct chain
    	if(flag_chain) {
    	    //cout<<"forward chain, mode "<<mode<<endl;
    	    modeindex=seedindex[mode].second;
       		_seed=seeds[modeindex][0];
    		//cout<<"mode:"<<mode<<" m="<<ref.index[_seed].n1<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<' '<<_seed<<endl;
    		if((m=ref.index[_seed].n1)>0){ 
    	    	_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
        		h=param.profile[modeindex][0].a;
        		for(j=0; j!=m; j++) {
    	    		_hit=_refloc[j];
    	    		if((_hit.chr>>16)!=modeindex) continue; // mode or strand not match
                            if(j+PREFETCH_LOOP<m){      
                                  prefetch_hit=_refloc[j+PREFETCH_LOOP];                            
                                  __builtin_prefetch(ref.bfa[prefetch_hit.chr&0xffff].s+(prefetch_hit.loc-h)/SEGLEN,0,0);                                                   
                            }                                                                                        
    	    		_hit.chr=_hit.chr&0xffff;
    	    		if(_hit.loc<h) continue; //underflow the start of refseq
    	    		_hit.loc-=h;
    	    		_hitz=_hit.loc%SEGLEN;
                    w=CountMismatch((bit64_t*)bseq[_hitz], (bit64_t*)reg[_hitz], (bit64_t*)(ref.bfa[_hit.chr].s+_hit.loc/SEGLEN));
    	    		//cout<<" j="<<j<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<endl;
    	    	 	if(w>snp_thres) continue;
    	    	 	if(_hit.chr&1) _hit.loc=ref.title[_hit.chr].rc_offset-_pread->seq.size()-_hit.loc;
    	    	 	if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
            		if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist
                    if(!param.pairend){
                        seg_info=ref.CCGG_seglen(_hit.chr, _hit.loc, _pread->seq.size()); //get fragment information
                        //cout<<"seg1:"<<seg_info.first<<" seg2:"<<seg_info.second<<endl;
                        if(seg_info.second>param.max_insert) continue; // fragent too large
                        if(seg_info.second<param.min_insert) continue; // fragment too small
                    }

    	    		hits[w][_cur_n_hit[w]++]=_hit;
    	    		if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return;
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) 
                        if(w==0) return; else snp_thres=w-1;
                }
    		}
    	}
    	//cout<<"end forward chain\n";
    	//complementary chain
    	if(cflag_chain) {
    	//cout<<"reverse chain, mode "<<mode<<endl;
    	    modeindex=cseedindex[mode].second;
            int cmodeindex=_pread->seq.size()/param.seed_size-1-modeindex;
    		_seed=cseeds[modeindex][0];
    		//cout<<"mode:"<<mode<<" modeindex="<<modeindex<<" m="<<ref.index[_seed].n1<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<endl;
    		if((m=ref.index[_seed].n1)>0){
    	    	_refloc=ref.index[_seed].loc1;  //list of seeded reference locations
    	    	h=param.profile[modeindex][0].a+cseed_offset;
    	    	for(j=0; j!=m; j++) {
    	    		_hit=_refloc[j];
                    if(((_hit.chr^0x1000000)>>16)!=cmodeindex) continue;  //mode or strand not match
                    if(j+PREFETCH_LOOP<m){
                            prefetch_hit=_refloc[j+PREFETCH_LOOP];
                            __builtin_prefetch(ref.bfa[prefetch_hit.chr&0xffff].s+(prefetch_hit.loc-h)/SEGLEN,0,0);
                    }
                    _hit.chr=_hit.chr&0xffff;
                    if(_hit.loc<h) continue; //underflow the start of refseq
    	    		_hit.loc-=h;
    	    		_hitz=_hit.loc%SEGLEN;
    	    		w=CountMismatch((bit64_t*)cbseq[_hitz], (bit64_t*)creg[_hitz], (bit64_t*)(ref.bfa[_hit.chr].s+_hit.loc/SEGLEN));
    	    		//cout<<" j="<<j<<" chr:"<<(int)_hit.chr<<" loc:"<<_hit.loc<<" mis:"<<w<<endl;
    	    		if(w>snp_thres) continue; //too many mismatches 
    	    	 	if(_hit.chr&1) _hit.loc=ref.title[_hit.chr].rc_offset-_pread->seq.size()-_hit.loc;
    	    	 	if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
    	    	 	 
    	    	 	if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist 
    	    		chits[w][_cur_n_chit[w]++]=_hit;
                    if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return;
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) if(w==0) return; else snp_thres=w-1;
                }
    		}
    	}
    	//cout<<"end reverse chain\n";
    }
    else{
    	//direct chain
    	if(flag_chain) {
    	modeindex=seedindex[mode].second;
    	//cout<<"forward chain, mode:"<<mode<<" modeindex:"<<modeindex<<endl;
        	for(i=0; i!=param.index_interval; i++) {
        		_seed=seeds[modeindex][i];
                if(ref.index2[_seed]==NULL) continue; //no match
        		m=ref.index2[_seed][0];  mc=ref.index2[_seed][1];
        		//cout<<" i= "<<i<<" m="<<m<<" mc="<<mc<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<endl;
        		h=-param.profile[modeindex][i].a+i-seed_start_array[modeindex];
        		for(j=2,_refloc2=ref.index2[_seed]+2; j!=mc; j++,_refloc2++) { 
       				__builtin_prefetch(ref.refcat+(*(_refloc2+PREFETCH_LOOP)+h)/SEGLEN,0,0);
        			_hit.loc=(*_refloc2)+h;
                    //cout<<" j="<<j<<" refloc:"<<(*_refloc2)<<" h:"<<h<<" loc:"<<_hit.loc<<endl;
                    w=CountMismatch((bit64_t*)bseq[_hit.loc%SEGLEN], (bit64_t*)reg[_hit.loc%SEGLEN], (bit64_t*)(ref.refcat+(_hit.loc/SEGLEN)));
        			//cout<<"id:"<<_hit.loc%SEGLEN<<" mis:"<<w<<endl;
        		 	if(w>snp_thres) continue;
                    _hit=ref.int2hit(_hit.loc,0);
        			//cout<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<" snp_thres:"<<snp_thres<<endl;                    
        		 	if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
               		if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist			hits[w][_cur_n_hit[w]++]=_hit;
                        hits[w][_cur_n_hit[w]++]=_hit;
                        if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return;
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) 
                        if(w==0) return; else snp_thres=w-1;
        		}
                //cout<<"#############\n";
        		for(; j!=m; j++,_refloc2++) { 
       				__builtin_prefetch(ref.crefcat+(*(_refloc2+PREFETCH_LOOP)+h)/SEGLEN,0,0);
        			_hit.loc=(*_refloc2)+h;
                    //cout<<" j="<<j<<" refloc:"<<(*_refloc2)<<" h:"<<h<<" loc:"<<_hit.loc<<endl;
                    w=CountMismatch((bit64_t*)bseq[_hit.loc%SEGLEN], (bit64_t*)reg[_hit.loc%SEGLEN], (bit64_t*)(ref.crefcat+(_hit.loc/SEGLEN)));
        			//cout<<" mis:"<<w<<endl;
        		 	if(w>snp_thres) continue;
                    _hit=ref.int2hit(_hit.loc,1);
        		 	_hit.loc=ref.title[_hit.chr].rc_offset-_pread->seq.size()-_hit.loc;
        		 	//cout<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<" snp_thres:"<<snp_thres<<endl;
        		 	if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
               		if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist			hits[w][_cur_n_hit[w]++]=_hit;
                        hits[w][_cur_n_hit[w]++]=_hit;
                        if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return; 
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) 
                        if(w==0) return; else snp_thres=w-1;
        		}

        	}
    	}
    	//complementary chain
    	if(cflag_chain) {
            modeindex=cseedindex[mode].second;
        	//cout<<"reverse chain, mode:"<<mode<<" modeindex:"<<modeindex<<endl;
        	for(i=0; i!=param.index_interval; i++) {
        		_seed=cseeds[modeindex][i];
        		//cout<<"mode:"<<mode<<" i= "<<i<<" m="<<ref.index[_seed].n1<<" seed:"<<_seed<<endl;
                if(ref.index2[_seed]==NULL) continue; //no match
        		m=ref.index2[_seed][0];  mc=ref.index2[_seed][1];
        		//cout<<" i="<<i<<" m="<<m<<" mc="<<mc<<" seed:"<<param.StrSeed(param.map3to4(_seed),param.seed_size)<<endl;
        		h=-param.profile[modeindex][i].a+i-cseed_start_array[modeindex];
        		for(j=2,_refloc2=ref.index2[_seed]+2; j!=mc; j++,_refloc2++) { 
       				__builtin_prefetch(ref.refcat+(*(_refloc2+PREFETCH_LOOP)+h)/SEGLEN,0,0);
        			_hit.loc=(*_refloc2)+h;
        			w=CountMismatch((bit64_t*)cbseq[_hit.loc%SEGLEN], (bit64_t*)creg[_hit.loc%SEGLEN], (bit64_t*)(ref.refcat+(_hit.loc/SEGLEN)));
        			//cout<<" i="<<i<<" loc:"<<_hit.loc<<" mis:"<<w<<endl;
        			if(w>snp_thres) continue;
                    _hit=ref.int2hit(_hit.loc,0);
        			//cout<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<" snp_thres:"<<snp_thres<<endl;
        			if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
              		if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist
        			chits[w][_cur_n_chit[w]++]=_hit;
        		if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return;	
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) 
                        if(w==0) return; else snp_thres=w-1;
        		}
        		for(; j!=m; j++,_refloc2++) { 
       				__builtin_prefetch(ref.crefcat+(*(_refloc2+PREFETCH_LOOP)+h)/SEGLEN,0,0);
        			_hit.loc=(*_refloc2)+h;
        			w=CountMismatch((bit64_t*)cbseq[_hit.loc%SEGLEN], (bit64_t*)creg[_hit.loc%SEGLEN], (bit64_t*)(ref.crefcat+(_hit.loc/SEGLEN)));
        			//cout<<" i="<<i<<" loc:"<<_hit.loc<<" mis:"<<w<<endl;
        			if(w>snp_thres) continue;
                    _hit=ref.int2hit(_hit.loc,1);
        			//cout<<" chr"<<(int)_hit.chr<<":"<<_hit.loc<<" mis:"<<w<<" snp_thres:"<<snp_thres<<endl;
        			_hit.loc=ref.title[_hit.chr].rc_offset-_pread->seq.size()-_hit.loc;
        			if(_hit.loc+_pread->seq.size()>ref.title[_hit.chr].size) continue; //overflow the end of refseq
              		if(!hitset[_hit.chr>>1].insert(_hit.loc).second) continue; //hit already exist
        			chits[w][_cur_n_chit[w]++]=_hit;
                    if(w==mode&&!param.pairend&&param.report_repeat_hits==0) if(_cur_n_hit[w]+_cur_n_chit[w]>1) return;
                    if(_cur_n_hit[w]+_cur_n_chit[w]>=param.max_num_hits) 
                        if(w==0) return; else snp_thres=w-1;
        		}

        	}
    	}    
    }
}

void SingleAlign::SortHits(int n)
{
	sort(hits[n], hits[n]+_cur_n_hit[n], HitComp2);
/*
}

void SingleAlign::SortcHits(int n)
{
*/
	sort(chits[n], chits[n]+_cur_n_chit[n], HitComp2);
}



void SingleAlign::SortHits4PE(int n)
{
    //cout<<"SORT "<<_pread->name.c_str()<<" "<<n<<"\t"<<_cur_n_hit[n]<<" "<<_cur_n_chit[n]<<endl;
	sort(hits[n], hits[n]+_cur_n_hit[n], HitComp);
	sort(chits[n], chits[n]+_cur_n_chit[n], HitComp);
}


int SingleAlign::TrimAdapter(){
    int i,j,k,m,m0, pos;
     
    raw_readlen=_pread->seq.size();
    if(param.RRBS_flag){
        for(i=0; i<param.n_adapter; i++){
            for(pos=param.seed_size;pos<_pread->seq.size()-5;pos++){
                m0=0;
                for(k=0,_readnt=_pread->seq.begin()+pos,_adapternt=param.adapter[i].begin();k<param.adapter[i].size()&&k<15&&_readnt!=_pread->seq.end();k++,++_readnt,++_adapternt){
    		        if((m0+=(*_adapternt!=*_readnt))>4) break;
                }
                if(k<m0*5) continue;
    	        m=m0;
                for(_adapternt=param.digest_site.begin(),_readnt=_pread->seq.begin()+pos-param.digest_site.size()+param.digest_pos;
                    _adapternt!=param.digest_site.end()-param.digest_pos;++_adapternt,++_readnt) {
    		        m+=(*_adapternt!=*_readnt)&&(*_adapternt!='C'||*_readnt!='T');
    	        }
                if(k>=m*5) {
                    _pread->seq.erase(pos); 
                    if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                    return 1;
                }
        
                if(param.pairend) { //pair-end
            		m=m0;
                    for(_adapternt=param.digest_site.begin(),_readnt=_pread->seq.begin()+pos-param.digest_site.size()+param.digest_pos;
                        _adapternt!=param.digest_site.end()-param.digest_pos;++_adapternt,++_readnt) {
    	        	    m+=(*_adapternt!=*_readnt)&&(*_adapternt!='G'||*_readnt!='A');
    		        }
                    if(k>=m*5) {
                        _pread->seq.erase(pos); 
                        if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                        return 1;
                    }
                }
            }
        }
    }
    else{
        for(i=0; i<param.n_adapter; i++){
            for(pos=param.seed_size;pos<_pread->seq.size()-4;pos++){
                m0=0;
                for(k=0,_readnt=_pread->seq.begin()+pos,_adapternt=param.adapter[i].begin();k<param.adapter[i].size()&&k<15&&_readnt!=_pread->seq.end();k++,++_readnt,++_adapternt){
                    if((m0+=(*_adapternt!=*_readnt))>4) break;
                }
                if(k>=m0*5&&k>3) {
                    _pread->seq.erase(pos); 
                    if(_pread->qual.size()>pos) _pread->qual.erase(pos);
                    return 1;
                }    
            }
        }    
    }
    return 0;
}


void SingleAlign::ClearHits()
{
    int i;
	for(i=0; i<=param.max_snp_num; i++) _cur_n_hit[i]=_cur_n_chit[i]=0;
	for(i=0; i<param.total_ref_seq; i++) hitset[i].clear(); 
}

int SingleAlign::RunAlign(RefSeq &ref)
{
    int i;
    //cout <<_pread->name.c_str()<< "\t" << _pread->seq.c_str() << _pread->qual.c_str() << endl;
	ClearHits();
    seedseg_num=min((int)((_pread->seq.size()-param.index_interval+1)/param.seed_size),(int)(read_max_snp_num+1));
	ConvertBinaySeq(); 
    snp_thres=read_max_snp_num;
    cseed_offset=_pread->seq.size()%param.seed_size;
    ReorderSeed(ref);
    for(i=0; i<seedseg_num; i++){
        //GenerateSeeds_1(i);
        SnpAlign(ref,i);      
        if(!param.RRBS_flag) for(int ii=0;ii<=i;ii++) if(_cur_n_hit[ii]||_cur_n_chit[ii]) return 1;
	}
    for(i=0; i<=read_max_snp_num; i++) if(_cur_n_hit[i]||_cur_n_chit[i]) return 1;         
    return 0;
}

void SingleAlign::ReorderSeed(RefSeq &ref) {
    bit32_t i,ii,s,cs, total=0xffffffff, ctotal=0xffffffff, tt;
    if(param.RRBS_flag) seed_start_offset=cseed_start_offset=0;
    else {
        ii=(_pread->seq.size()-param.index_interval+1)%param.seed_size;
        for(i=0;i<ii;i++) {
            if(flag_chain) {
                tt=GetTotalSeedLoc(ref, i);
                if(tt<total) {total=tt; seed_start_offset=i;}
            }
            if(cflag_chain) {
                tt=GetTotalCSeedLoc(ref, i);
                if(tt<ctotal) {ctotal=tt; cseed_start_offset=i;}
            }
        }
    }
    //cout<<"start="<<seed_start_offset<<" cstart="<<cseed_start_offset<<endl;
    if(flag_chain) {
        AdjustSeedStartArray(ref);
        seedindex.clear();
        for(i=0; i<seedseg_num; i++){
            GenerateSeeds(i, seed_start_array[i]);  
            s=0;
            if(param.RRBS_flag) s+=ref.index[seeds[i][0]].n1;
            else {
                for(ii=0; ii!=param.index_interval; ii++)
                    if(ref.index2[seeds[i][ii]]!=NULL) s+=ref.index2[seeds[i][ii]][0];
            }
            seedindex.push_back(pair<int,int>(s,i)); 
        }
        //for(vector<pair<int,int> >::iterator sit=seedindex.begin();sit!=seedindex.end();++sit) cout <<sit->first<<" "<<sit->second<<"   "; cout<<endl;
        sort(seedindex.begin(), seedindex.end());
        //for(vector<pair<int,int> >::iterator sit=seedindex.begin();sit!=seedindex.end();++sit) cout <<sit->first<<" "<<sit->second<<"   "; cout<<endl;
    }
    if(cflag_chain) {
        AdjustCSeedStartArray(ref);
        cseedindex.clear();
        for(i=0; i<seedseg_num; i++){
            GenerateCSeeds(i, cseed_start_array[i]);  
            s=0;
            if(param.RRBS_flag) s+=ref.index[cseeds[i][0]].n1;
            else {
                for(ii=0; ii!=param.index_interval; ii++)
                    if(ref.index2[cseeds[i][ii]]!=NULL) s+=ref.index2[cseeds[i][ii]][0];
            }
            cseedindex.push_back(pair<int,int>(s,i));
        }
        sort(cseedindex.begin(), cseedindex.end());    
        //for(vector<sort_index>::iterator sit=cseedindex.begin();sit!=cseedindex.end();++sit) cout <<sit->val<<" "<<sit->index<<"   "; cout<<endl;
    }
}

void SingleAlign::AdjustSeedStartArray(RefSeq &ref) {
    bit32_t i, ii, tt, total;
    int ptr, start, end, max_offset;
    //cout<<"max offset:"<<max_offset<<" seedseg_num="<<seedseg_num<<" seed_start_offset="<<seed_start_offset<<endl;
    for(i=0;i<seedseg_num;i++) seed_start_array[i]=seed_start_offset;
    if(param.RRBS_flag) return;
    max_offset=(_pread->seq.size()-param.index_interval+1)%param.seed_size; 
    for(i=0;i<seedseg_num;i++) {
        if(i%2==0) ptr=i/2; else ptr=seedseg_num-1-i/2;
        total=0xffffffff; 
        //cout<<"i="<<i<<" ptr="<<ptr<<endl;
        if(ptr==0) start=0; else start=seed_start_array[ptr-1];
        if(ptr==seedseg_num-1) end=max_offset; else end=seed_start_array[ptr+1];
        //cout<<"start="<<start<<" end="<<end<<endl;
        seed_start_array[ptr]=start;
        for(ii=start;ii<=end;ii++) {
            tt=CountSeeds(ref,ptr,ii);
            //cout<<"ii="<<ii<<" tt="<<tt<<endl;
            if(tt<total) {total=tt; seed_start_array[ptr]=ii;}
        }
        //for(ii=0;ii<seedseg_num;ii++) cout<<ii<<":"<<seed_start_array[ii]<<" "; cout<<endl;
    }
}    

void SingleAlign::AdjustCSeedStartArray(RefSeq &ref) {
    bit32_t i, ii, tt, total;
    int ptr, start, end, max_offset;
    for(i=0;i<seedseg_num;i++) cseed_start_array[i]=cseed_start_offset;
    if(param.RRBS_flag) return;
    max_offset=(_pread->seq.size()-param.index_interval+1)%param.seed_size;
    for(i=0; i<seedseg_num; i++){
        if(i%2==0) ptr=i/2; else ptr=seedseg_num-1-i/2;
        total=0xffffffff; 
        if(ptr==0) start=0; else start=cseed_start_array[ptr-1];
        if(ptr==seedseg_num-1) end=max_offset; else end=cseed_start_array[ptr+1];
        cseed_start_array[ptr]=start;
        for(ii=start;ii<=end;ii++) {
            tt=CountCSeeds(ref,ptr,ii);
            if(tt<total) {total=tt; cseed_start_array[ptr]=ii;}
        }                                                                     
    }
}

int SingleAlign::CountSeeds(RefSeq &ref, int n, int start) {
	int i, total=0; bit32_t s;
    for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) {
        s=seed_array[_pro->a+start-i];
        if(ref.index2[s]!=NULL) total+=ref.index2[s][0];
    }
    return total;
}

int SingleAlign::CountCSeeds(RefSeq &ref, int n, int cstart) {
	int i, total=0; bit32_t s;
    for(i=0,_pro=param.profile[n];i<param.index_interval;i++,_pro++) {
        s=cseed_array[_pro->a+cstart-i];
        if(ref.index2[s]!=NULL) total+=ref.index2[s][0];
    }
    return total;
}

bit32_t SingleAlign::GetTotalSeedLoc(RefSeq &ref, int start) {
    int i, total=0;
    for(i=0; i<seedseg_num; i++) total+=CountSeeds(ref, i, start);
    return total;
}

bit32_t SingleAlign::GetTotalCSeedLoc(RefSeq &ref, int start) {
    int i, total=0;
    for(i=0; i<seedseg_num; i++) total+=CountCSeeds(ref, i, start);
    return total;
}

int SingleAlign::FilterReads()
{
    //if(_pread->bam_flag&0x200) return 1;
    TrimAdapter();
    if(TrimLowQual()==0) return 1;
	if(_pread->seq.size()<param.min_read_size) return 1;
    if(CountNs()>param.max_ns) return 1;
    read_max_snp_num=(param.max_snp_num+1)*(_pread->seq.size()-1)/raw_readlen;
    //set_RRBS_start();
	return 0;
}

void SingleAlign::Do_Batch(RefSeq &ref)
{
	_str_align.clear(); 
	bit32_t tt;
	//alignment for normal sequences
	for(_pread=mreads.begin(), tt=0; tt<num_reads; _pread++, tt++) {
        //cout <<_pread->name.c_str()<< "\t" << _pread->seq.c_str() << _pread->qual.c_str() << endl;
		if(FilterReads()) {
			if(param.report_repeat_hits) s_OutHit(0, -1, 0, hits[0], 0, ref, _str_align);
        }
        else {
    		RunAlign(ref);
        	StringAlign(ref, _str_align);
        }
	}
}


//output align hits
void SingleAlign::StringAlign(RefSeq &ref, string &os)
{
	//Reverse_Seq();
	//Reverse_Qual();
	int ii, sum=0, j;
	//snp align:
    
    //cout<<_pread->name<<"\t"<<_pread->seq<<"\t";
    for(ii=0; ii<=read_max_snp_num; ii++) if((sum=_cur_n_hit[ii]+_cur_n_chit[ii])>0) break;
    //cout<<sum<<endl;
    
    if(sum==0) s_OutHit(0, sum, ii, hits[0], 0, ref, os);
    else {
        j=myrand(_pread->index,&rand_rSeed)%sum;
   		if(j<_cur_n_hit[ii]) s_OutHit(0, sum, ii, &hits[ii][j], 0, ref, os);
   		else s_OutHit(1, sum, ii, &chits[ii][j-_cur_n_hit[ii]], 0, ref, os);
   	}
}

//write output according to types of hits
/* n: # of hits; chain: 0+/1-; flag: class of read; sig: 1, detect snp sites, 0, not */
void SingleAlign::s_OutHit(int chain, int n, bit8_t nsnps, Hit *hit, int insert_size, RefSeq &ref, string &os)
{
    bit32_t ii, jj, hitloc;
    if(param.out_sam){ //output in .sam format
        int flag,seg_len;
        flag=0x40*_pread->readset;
        ref_loc_t seg_start;
        if(n<0){
            if(!param.out_unmap) return;
            flag|=0x204; //QC
          	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);
        }
        else if(n==0){
            if(!param.out_unmap) return;
            flag|=0x4; //NM
           	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);
        }
        else if(n>1&&param.report_repeat_hits==0){
            if(!param.out_unmap) return;
            flag|=0x104; //NM
           	sprintf(_ch,"%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",_pread->name.c_str(),flag,_pread->seq.c_str(),_pread->qual.c_str()); 
            os.append(_ch);        
        }
        else{
            n_aligned++;
            if(n==1) flag|=0x0; //UM
            else flag|=0x100; //MA & OF
            if((chain^(hit->chr%2))&&n){
                flag|=0x010; //reverse read seq
               	reverse(_pread->seq.begin(), _pread->seq.end());
                for(ii=0;ii<_pread->seq.size();ii++) _pread->seq[ii]=rev_char[_pread->seq[ii]];
    	    	reverse(_pread->qual.begin(), _pread->qual.end());
            }
            sprintf(_ch,"%s\t%d\t%s\t%u\t255\t%dM\t*\t0\t0\t%s\t%s\tNM:i:%d",_pread->name.c_str(),flag,ref.title[hit->chr].name.c_str(),hit->loc+1,_pread->seq.size(),_pread->seq.c_str(),_pread->qual.c_str(),nsnps); 
            
            os.append(_ch);

            if(param.out_ref) {
                int ptr=0;
                for(ii=2;ii>0;ii--,ptr++) {
                    if(hit->loc<ii) continue;
                    _mapseq[ptr]=param.useful_nt[*(ref.bfa[(hit->chr>>1)<<1].s+(hit->loc-ii)/SEGLEN)>>(SEGLEN*2-2-((hit->loc-ii)%SEGLEN)*2)&0x3]+32;
                }
                for(ii=0;ii<_pread->seq.size()+2;ii++,ptr++) {
                    _mapseq[ptr]=param.useful_nt[*(ref.bfa[(hit->chr>>1)<<1].s+(hit->loc+ii)/SEGLEN)>>(SEGLEN*2-2-((hit->loc+ii)%SEGLEN)*2)&0x3];
                }
                _mapseq[ptr]=0; _mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;
                sprintf(_ch, "\tXR:Z:%s",_mapseq);
                os.append(_ch);
            }
            
            if(param.RRBS_flag){
                seg_info=ref.CCGG_seglen(hit->chr, hit->loc, _pread->seq.size());
                sprintf(_ch,"\tZP:i:%d\tZL:i:%d",seg_info.first,seg_info.second); 
                os.append(_ch);
            }

            sprintf(_ch,"\tZS:Z:%c%c\n",chain_flag[hit->chr%2], chain_flag[chain]); 
            os.append(_ch);

            if((chain^(hit->chr%2))&&(n>1)){
                reverse(_pread->seq.begin(), _pread->seq.end());
                for(ii=0;ii<_pread->seq.size();ii++) _pread->seq[ii]=rev_char[_pread->seq[ii]];
                reverse(_pread->qual.begin(), _pread->qual.end());
            }                                                                                                                
        }
    }

    else{ //output in .bsp format    
        if(!param.out_unmap&&(n<=0||(n>1&&param.report_repeat_hits==0))) return;
        sprintf(_ch, "%s\t", _pread->name.c_str());
    	os.append(_ch);
    	
    	if((chain^(hit->chr%2))&&n) {
    	    	reverse(_pread->seq.begin(), _pread->seq.end());
    	    	for(ii=0;ii<_pread->seq.size();ii++) _pread->seq[ii]=rev_char[_pread->seq[ii]];
    	    	reverse(_pread->qual.begin(), _pread->qual.end());
        }

    	sprintf(_ch, "%s\t%s\t", _pread->seq.c_str(), _pread->qual.c_str());
    	os.append(_ch);
    	
    	if (n<0) os.append("QC");
    	else if(n==0) os.append("NM"); 
    	else if(n==1) os.append("UM");
    	else if(n>=param.max_num_hits) os.append("OF");
    	else os.append("MA");
        
        if((n>0&&param.report_repeat_hits==1)||(n==1&&param.report_repeat_hits==0)){
            n_aligned++;
            int ptr=0;
            for(ii=2;ii>0;ii--,ptr++) {
                if(hit->loc<ii) continue;
                _mapseq[ptr]=param.useful_nt[*(ref.bfa[(hit->chr>>1)<<1].s+(hit->loc-ii)/SEGLEN)>>(SEGLEN*2-2-((hit->loc-ii)%SEGLEN)*2)&0x3]+32;
            }
            for(ii=0;ii<_pread->seq.size()+2;ii++,ptr++) {
                _mapseq[ptr]=param.useful_nt[*(ref.bfa[(hit->chr>>1)<<1].s+(hit->loc+ii)/SEGLEN)>>(SEGLEN*2-2-((hit->loc+ii)%SEGLEN)*2)&0x3];
            }
            _mapseq[ptr]=0; _mapseq[ptr-1]+=32; _mapseq[ptr-2]+=32;

        	sprintf(_ch, "\t%s\t%u\t%c%c\t%d\t%s\t%d\t", ref.title[hit->chr].name.c_str(), hit->loc+1, chain_flag[hit->chr%2],  chain_flag[chain], insert_size, _mapseq, nsnps);
            os.append(_ch); 
            for(ii=0; ii<read_max_snp_num; ii++){
                sprintf(_ch, "%d:", _cur_n_hit[ii]+_cur_n_chit[ii]);
                os.append(_ch);
            }
            sprintf(_ch, "%d", _cur_n_hit[ii]+_cur_n_chit[ii]);
            os.append(_ch);
        }
    	os.append("\n");
        if((chain^(hit->chr%2))&&n) {
            reverse(_pread->seq.begin(), _pread->seq.end());
            for(ii=0;ii<_pread->seq.size();ii++) _pread->seq[ii]=rev_char[_pread->seq[ii]];
            reverse(_pread->qual.begin(), _pread->qual.end());
        }
                                                                

    }	

    //1	QNAME	Query (pair) NAME
    //2	FLAG	bitwise FLAG
    //3	RNAME	Reference sequence NAME
    //4	POS	1-based leftmost POSition/coordinate of clipped sequence
    //5	MAPQ	MAPping Quality (Phred-scaled)
    //6	CIAGR	extended CIGAR string
    //7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
    //8	MPOS	1-based Mate POSistion
    //9	ISIZE	Inferred insert SIZE
    //10	SEQ	query SEQuence on the same strand as the reference
    //11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
    //12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE

}


void SingleAlign::Fix_Unpaired_Short_Fragment(RefSeq &ref) {
    int ii, j, k;
    if(_pread->seq.size()>=param.min_insert) return;
    for(ii=0; ii<=read_max_snp_num; ii++) {
        for(j=0;j<_cur_n_hit[ii];j++) {
            seg_info=ref.CCGG_seglen(hits[ii][j].chr, hits[ii][j].loc, _pread->seq.size());
            if(seg_info.second<param.min_insert||seg_info.second>param.max_insert) {
                _cur_n_hit[ii]--; 
                for(k=j;k<_cur_n_hit[ii];k++) hits[ii][k]=hits[ii][k+1];
                j--;
            }
        }
        for(j=0;j<_cur_n_chit[ii];j++) {
            seg_info=ref.CCGG_seglen(chits[ii][j].chr, chits[ii][j].loc, _pread->seq.size());
            if(seg_info.second<param.min_insert||seg_info.second>param.max_insert) {
                _cur_n_chit[ii]--; 
                for(k=j;k<_cur_n_chit[ii];k++) chits[ii][k]=chits[ii][k+1];
                j--;
            }
        }
        if(_cur_n_hit[ii]+_cur_n_chit[ii]>0) break;
    }

}
