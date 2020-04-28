#ifndef SAOPT_H_
#define SAOPT_H_

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <pthread.h>

using namespace std;

struct Options {
  string typ, seq, ref, prefix;
  int thread, mismatch, seed;
  bool all, debug, help;
};

struct Sequence {
  string nam, ref;
};

struct Alignment {
  string rname, cigar, tag;
  int flag, pos;
  float mapq;
};

struct Statistics {
  float score;
  int total, mapped;
};
#endif
