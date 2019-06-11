#include "sasub.h"

void _usage()
{
  printf ("\nProgram: sa (smRNA Aligner Based on Needleman-Wunsch)\n");
  printf ("Contact: Qianhua ZHU <zhuqianhua@bgi.com>\n\n");
  printf ("Usage  : sa [options]\n\n");
  printf ("Options: -s, --seq    *<s>  small rna sequence in fasta format\n");
  printf ("         -r, --ref    *<s>  reference sequence in fasta format\n");
  printf ("         -p, --prefix  <s>  prefix of output, default ./sa \n");
  printf ("         -m, --mismatch<i>  number of mismatch, default 1\n");
  printf ("         -d, --seed    <i>  length of seed, default 9\n");
  printf ("         -t, --thread  <i>  number of threads, default 1\n");
  printf ("         -a, --all     <b>  report all best alignments, default random one\n");
  printf ("         -b, --debug   <b>  print the debug information\n");
  printf ("         -h, --help    <b>  print this information\n\n");
  exit(0);
}

void _getopt(int argc, char* argv[])
{
  char const * shortOpt = "s:r:p:t:m:d:abh";
  struct option longOpt[] = {
    {"seq", 1, NULL, 's'},
    {"ref", 1, NULL, 'r'},
    {"prefix", 1, NULL, 'p'},
    {"mismatch", 1, NULL, 'm'},
    {"seed", 1, NULL, 'd'},
    {"all", 1, NULL, 'a'},
    {"thread", 1, NULL, 't'},
    {"debug", 0, NULL, 'b'},
    {"help", 0, NULL, 'h'},
    {NULL, 0, NULL, 0},
  };
  int nextOpt;
  while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1) 
  {
    switch (nextOpt) 
    {
      case 's':
        opts.seq = optarg;
        break;
      case 'r':
        opts.ref = optarg;
        break;
      case 'p':
        opts.prefix = optarg;
        break;
      case 't':
        opts.thread = atoi(optarg);
        break;
      case 'm':
        opts.mismatch = atoi(optarg);
        break;
      case 'd':
        opts.seed = atoi(optarg);
        break;
      case 'a':
        opts.all = true;
        break;
      case 'b':
        opts.debug = true;
        break;
      case 'h':
        opts.help = true;
        break;
    }
  }
  if (opts.seq.compare("") == 0 or opts.ref.compare("") == 0)
    _usage();
  if (opts.thread > THREAD_NUMBER)
    opts.thread = THREAD_NUMBER;
}

void _loadref()
{
  ifstream infile(opts.ref.c_str());
  string line;
  int mark = 0;
  Sequence seq = {"", ""};
  while (!infile.eof())
  {
    getline(infile, line);
    if (line.find_first_of('>') != string::npos) {
      int end = line.length();
      if (line.find_first_of(' ') != string::npos)
        end = line.find_first_of(' ');
      if (seq.nam.compare("") != 0) { 
        transform(seq.ref.begin(), seq.ref.end(), seq.ref.begin(), ::toupper);
        vref.push_back(seq);
        for (int i = 0; i <= seq.ref.size() - opts.seed; i ++) 
          mref[mark % opts.thread][seq.ref.substr(i, opts.seed)].push_back(vref.size()-1);
        ofile << "@SQ\tSN:" << seq.nam << "\tLN:" << seq.ref.size() << endl;
        seq.ref = "";
      }
      seq.nam = line.substr(line.find_first_of('>') + 1, end - line.find_first_of('>') - 1);
      mark ++;
    } else {
      seq.ref += line;
    }
  }
  infile.close();
  transform(seq.ref.begin(), seq.ref.end(), seq.ref.begin(), ::toupper);
  vref.push_back(seq);
  for (int i = 0; i <= seq.ref.size() - opts.seed; i ++) 
    mref[mark % opts.thread][seq.ref.substr(i, opts.seed)].push_back(vref.size()-1);
  ofile << "@SQ\tSN:" << seq.nam << "\tLN:" << seq.ref.size() << endl;
}

bool _comp(const Alignment &a, const Alignment &b)
{
  if (a.mapq > b.mapq)
    return true;
  else
    return false;
}

template <class T>
string _itos(T n)
{
  stringstream ss;
  ss << n;
  return ss.str();
}

void _align() 
{
  stas.score = 0.0;
  int THREAD[opts.thread];
  pthread_t thread[opts.thread]; 
  for (int i = 0; i < opts.thread; i++) {
    THREAD[i] = i;
    int res = pthread_create(&thread[i], NULL, needleman_wunsch, (void *) &THREAD[i]);
    if (res != 0) {
      printf("Thread create %d failed\n", i);
      exit(0);
    }
  }
  while (1)
  {
    int mak = 0; 
    void *status;
    for (int i = 0; i < opts.thread; i++)
      if (pthread_join(thread[i], &status) != 0)
        mak ++;
    if (mak == 0)
      break;
  }
  sort(score.begin(), score.end(), _comp);
  int h0 = 0;
  string sa(" SA:Z:"), ln("");
  for (vector<Alignment>::iterator iter = score.begin(); iter != score.end(); ++iter) {
    if (iter->mapq < stas.score)
      break;
    if (ln.compare("") == 0) {
      ln  = seqs.nam + "\t" + _itos(iter->flag) + "\t" + iter->rname + "\t" + _itos(iter->pos) + "\t" + _itos(iter->mapq);
      ln += "\t" + iter->cigar + "\t*\t0\t0\t";
      if (iter->flag != 16) 
        ln += seqs.ref + "\t*\t" + iter->tag;
      else 
        ln += _revcomp(seqs.ref) + "\t*\t" + iter->tag;
      if (iter->flag != 4)
        stas.mapped ++;
      stas.total ++;
    } else if (iter->flag != 4) 
      sa += iter->rname + "," + _itos(iter->pos) + "," + _itos(iter->flag) + "," + _itos(iter->mapq) + "," + iter->tag + ";";
    if (iter->flag != 4)
      h0 ++;
  }
  vector<Alignment>().swap(score);
  ofile << ln << " H0:i:" << h0;
  if (opts.all == true && sa.size() > 6)
    ofile << sa;
  ofile << endl;
}

