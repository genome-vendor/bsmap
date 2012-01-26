#ifndef _READS_H_
#define _READS_H_

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include "param.h"
#include "sam.h"

using namespace std;

const int BatchNum=50000;

struct ReadInf
{
	bit32_t index;
	bit32_t readset; //added by yxi
	string name;
	string seq;
	string qual;
};

class ReadClass
{
public:
	ReadClass();
	void CheckFile(ifstream &fin, string filename, int readset);
	void InitialIndex();
	int LoadBatchReads(ifstream &fin, int mode);
public:
	vector<ReadInf> mreads;
	bit32_t num;
	int _file_format;  //0: fq; 1: fa; 2:BAM
	bit32_t _index;
    
    //added by yxi, for BAM input support
    samfile_t *SAM_fp;
    bam1_t *SAM_b;
};

#endif //_READS_H_
