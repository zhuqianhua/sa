#include "sa.h"

int main(int argc, char* argv[])
{
  _getopt(argc, argv);
  string oaln = opts.prefix + ".sam", ostt = opts.prefix + "_stat.xls";
  ofile.open(oaln.c_str());
  ofile << "@HD\tVN:1.5\tSO:unsorted" << endl;
  _loadref();
  if (opts.typ.compare("q") == 0)
    _fastq();
  ifstream iseq(opts.seq.c_str());
  string line;
  while (!iseq.eof())
  {
    getline(iseq, line);
    if (line.find_first_of('>') != string::npos) {
      int end = line.length();
      if (line.find_first_of(' ') != string::npos)
        end = line.find_first_of(' ');
      if (seqs.nam.compare("") != 0) {
        transform(seqs.ref.begin(), seqs.ref.end(), seqs.ref.begin(), ::toupper);
        _align(); 
        seqs.ref = ""; 
      }
      seqs.nam = line.substr(line.find_first_of('>') + 1, end - line.find_first_of('>') - 1);
    } else {
      seqs.ref += line;
    }
  }
  iseq.close();
  transform(seqs.ref.begin(), seqs.ref.end(), seqs.ref.begin(), ::toupper);
  _align();
  ofile.close();
  ofstream ostat;
  ostat.open(ostt.c_str());
  ostat << "Total reads\t" << stas.total << endl;
  ostat << "Mapped reads\t" << stas.mapped << endl;
  char rate[10];
  sprintf(rate, "%.2f", float(stas.mapped * 100 / stas.total));
  ostat << "Mapping rate\t" << rate << "%" << endl;
  ostat.close();
  return 0;
}

