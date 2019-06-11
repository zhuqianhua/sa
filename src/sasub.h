#ifndef SASUB_H_
#define SASUB_H_

#include <fstream>
#include <iostream>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include "saopt.h"
#include "sanw.h"
#define THREAD_NUMBER 16

void _usage();
void _getopt(int argc, char* argv[]);
void _loadref();
void _align();
string _revcomp(string &str);

extern ofstream ofile;
extern Statistics stas;
extern Options opts;
extern Sequence seqs;
extern vector<Alignment> score;
extern vector<Sequence> vref;
extern map<int, map<string, vector<int> > > mref;
#endif
