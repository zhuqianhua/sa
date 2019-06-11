#include "sanw.h"

int SCORE[3] = {1, -2, -5}; // match mismatch gap

template <class T>
string _itos(T n)
{
  stringstream ss;
  ss << n;
  return ss.str();
}

int _max(int x, int y, int z)
{
  int t = x > y ? x : y;
  return t > z ? t : z;
}

void _cigar(string& cig, Alignment* align)
{
  align->cigar = "";
  char letter = ' ';
  string tag;
  int count = 0, match = 0;
  for (int i = 0; i < cig.size(); i++) {
    if (letter != cig[i]) {
      if (letter != ' ')
        align->cigar += _itos(count) + _itos(letter);
        if (letter == 'M')
          match += count;
        else if (letter == 'X') {
          tag += _itos(match) + _itos(seqs.ref[i]);
          match = 0;
        } else if (letter == 'I') {
          tag += _itos(match) + "^";
          match = 0;
        }
      letter = cig[i];
      count  = 1;
    } else {
      count ++;
    }
  }
  if (letter == 'M') 
    match += count;
  if (letter == 'X') 
    tag += _itos(match) + _itos(seqs.ref[-1]);
  else if (letter == 'I') 
    tag += _itos(match) + "^";
  else 
    tag += _itos(match);
  align->cigar += _itos(count) + _itos(letter);
  align->tag = "MD:Z:" + tag;
}

string _revcomp(string &str)
{
  map<char, char> base;
  base['A'] = 'T'; 
  base['T'] = 'A'; 
  base['C'] = 'G'; 
  base['G'] = 'C'; 
  base['N'] = 'N';
  string revc(str);
  for (string::iterator it = revc.begin(); it != revc.end(); ++it)
    if (base.find(*it) != base.end())
      *it = base[*it];
  return revc;
}

void *needleman_wunsch(void *args)
{
  int p = *(int *) args;
  int avr = (seqs.ref.size() - 4) / (opts.mismatch + 1);
  vector<int> vpos; // position of the reference in vref
  for (int i = 0; i < opts.mismatch+1; i ++) {
    int beg = avr * i + 2;
    if (seqs.ref.size() - 2 - beg < opts.seed)
      beg = seqs.ref.size() - 2 - opts.seed;
    string seed = seqs.ref.substr(beg, opts.seed);
    if (mref[p].find(seed) != mref[p].end())
      vpos.insert(vpos.end(), mref[p][seed].begin(), mref[p][seed].end());
  }
  sort(vpos.begin(), vpos.end());
  vpos.erase(unique(vpos.begin(), vpos.end()), vpos.end());
  srand((unsigned int)(time(NULL)));
  random_shuffle(vpos.begin(), vpos.end());
  vector<Alignment> vnw;
  Alignment align = {"*", "*", "", 4, 0, 0};
  stas.score = align.mapq;
  vnw.push_back(align);
  for (int i = 0; i < vpos.size(); i++) {
    if (stas.score >= seqs.ref.size() * SCORE[0] && opts.all != true)
      break;
    align.rname = "*";
    align.cigar = "*";
    align.tag = "";
    align.flag = 4;
    align.pos = 0;
    align.mapq = 0.0;
    needleman_wunsch_align(vref[vpos[i]].ref, &align);
    if (align.pos != 0) {
      align.rname = vref[vpos[i]].nam;
      align.flag = 0;
    } 
    if (align.mapq >= stas.score && align.pos != 0) {
      vnw.push_back(align);
      stas.score = align.mapq;
    }
    if (align.pos == 0 && stas.score < seqs.ref.size() * SCORE[0]) {
      align.pos = 0;
      align.mapq = 0;
      string rev = _revcomp(vref[vpos[i]].ref);
      needleman_wunsch_align(rev, &align);
      if (align.pos != 0) {
        align.rname = vref[vpos[i]].nam;
        align.flag = 16;
      }
      if (align.mapq >= stas.score && align.pos != 0) {
        vnw.push_back(align);
        stas.score = align.mapq;
      }
    }
  }
  pthread_mutex_lock(&mutex);
  score.insert(score.end(), vnw.begin(), vnw.end());
  pthread_mutex_unlock(&mutex);
}

int needleman_wunsch_align(string &ref, Alignment *align)
{
  int matrix[seqs.ref.size()+1][ref.size()+1];
  // initialize the first row and colunm
  for (int r = 0; r <= seqs.ref.size(); r ++)
    matrix[r][0] = 0;
  for (int c = 0; c <= ref.size(); c ++)
    matrix[0][c] = 0;
  // caculate the scores based on F(i,j) = max{F(i-1,j-1)+s(x,y), F(i-1,j)+g, F(i, j-1)+g}
  int max_score = 0, end_row = 0, end_col = 0; 
  for (int r = 1; r <= seqs.ref.size(); r ++) {
    for (int c = 1; c <= ref.size(); c ++) {
      int val = SCORE[0];
      if (seqs.ref[r-1] != ref[c-1])
        val = SCORE[1];
      matrix[r][c] = _max(matrix[r-1][c-1]+val, matrix[r-1][c]+SCORE[2], matrix[r][c-1]+SCORE[2]);
      if (matrix[r][c] > max_score) {
        max_score = matrix[r][c];
        end_row = r;
        end_col = c;
      }
    }
  }
  if (opts.debug == true) {
    for (int r = 1; r <= seqs.ref.size(); r ++) {
      for (int c = 1; c <= ref.size(); c ++)
        cout << matrix[r][c] << " ";
      cout << endl;
    }
  }
  if (end_row == 0 || end_col == 0 || max_score <= 0)
    return 0;
  // align details
  string cigar = "M";
  int row = end_row, col = end_col, match = 0, mismatch = 0, clip = 0;
  if (seqs.ref[row-1] != ref[col-1]) {
    cigar = "S";
    clip ++;
  } else {
    match = 1;
    align->mapq += SCORE[0];
  }
  int count = 0;
  while (1)
  {
    if (count >= seqs.ref.size()-1)
      break;
    count ++;
    string type;
    if (opts.debug == true)
      cout << "road: " << row << " " << col << " ";
    if (matrix[row-1][col-1] >= matrix[row][col-1] && matrix[row-1][col-1] >= matrix[row-1][col]) {
      if (seqs.ref[row-2] != ref[col-2]) {
        if (row <= 3 || row >= seqs.ref.size()-1) {
          type = 'S';
          clip ++;
        } else {
          type = "X";
          mismatch ++;
          align->mapq += SCORE[1];
        }
      } else {
        type = "M";
        match ++;
        align->mapq += SCORE[0];
      }
      row --;
      col --;
    } else if (matrix[row][col-1] > matrix[row-1][col-1] && matrix[row][col-1] > matrix[row-1][col]) {
      if (row <= 3) {
        type = 'S';
        clip ++;
      } else {
        type = "D";
        mismatch ++;
        align->mapq += SCORE[2];
      }
      col --;
    } else {
      type = "I";
      row --;
      mismatch ++;
      align->mapq += SCORE[2];
    }
    if (opts.debug == true)
      cout << type << endl;
    if (mismatch > opts.mismatch)
      return 0;
    cigar = type + cigar;
    if (row == 1 or col == 1)
      break;
  }
  if (opts.debug == true)
    cout << "end_row: " << end_row << " count: " << count << " seq: " << seqs.ref.size() << endl;
  for (int i = 0; i < seqs.ref.size() - end_row; i++) {
    cigar += "S";
    clip ++;
  }
  for (int i = 0; i < end_row - count - 1; i++) {
    cigar = "S" + cigar;
    clip ++;
  }
  if (clip > 4) 
    return 0;
  align->mapq += float(match/ref.size());
  _cigar(cigar, align);
  align->pos = col;
}

