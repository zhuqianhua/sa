#ifndef SANW_H_
#define SANW_H_

#include <iostream>
#include <sstream>
#include <string.h>
#include "saopt.h"
#include <algorithm>
#include <time.h>

using namespace std;

void *needleman_wunsch(void *args);
int needleman_wunsch_align(string &ref, Alignment* align);

extern pthread_mutex_t mutex;
extern Statistics stas;
extern Options opts;
extern Sequence seqs;
extern map<int, map<string, vector<int> > > mref;
extern vector<Alignment> score;
extern vector<Sequence> vref;
#endif
