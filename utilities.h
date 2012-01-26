#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include<time.h>
#include "param.h"

using namespace std;

time_t Initial_Time();
//time used during the past step
time_t Cal_StepTime();
//total time exhaust
time_t Cal_AllTime();
//current time on string format
char * Curr_Time();

bit32_t myrand(int i);
bool HitComp(Hit a, Hit b);
bool HitComp2(Hit a, Hit b);
bool HitCompChr(Hit a, Hit b);
void disp_bfa(bit32_t a);
void disp_bfa64(bit64_t a);

#endif //_UTILITIES_H_
