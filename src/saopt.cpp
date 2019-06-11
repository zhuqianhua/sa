#include "saopt.h"

Options opts = {"", "", "./sa", 
                1, 1, 9,
               };

pthread_mutex_t mutex;
ofstream ofile;
Sequence seqs = {"",};
Statistics stas = {0, 0, 0};
vector<Sequence> vref;
map<int, map<string, vector<int> > > mref;
vector<Alignment> score;
