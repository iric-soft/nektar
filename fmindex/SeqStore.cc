//                              -*- Mode: C++ -*-
// SeqStore.cc
//

#include "SeqStore.h"


// LIFECYCLE ------------------------------------------------------------

SeqStore::SeqStore(BaseRep *rep)
  : mBaseRep(rep), mStore(NULL), mN(0), mOwnStore(false), mUnfolded(false),
    mInsertPoint(0) {}

SeqStore::SeqStore(BaseRep *rep, size_t n) {
  mBaseRep = rep;
  mN = n;
  mStore = new seq_t[mN];
  for (size_t i = 0; i < n; ++i)
    mStore[i] = 0;
  mOwnStore = true;
  mUnfolded = false;
  mInsertPoint = SeqIter(0);
}

SeqStore::SeqStore(SeqStore &o, SeqRange &r) {
  // Creates a store as a "viewport" to o.  The Range r is transformed to be used
  // in the local coordinate system of this new store.

  int dir = r.dir();
  if (dir < 0) r.reverse();
  size_t wbegin = r.begin().getOffset() / mBaseRep->basePerWord;
  size_t wend = (r.end().getOffset() - 1) / mBaseRep->basePerWord + 1;

  mN = wend - wbegin;
  mStore = o.mStore + wbegin;

  r = SeqRange (r.begin().getOffset() - wbegin * mBaseRep->basePerWord,
                r.end().getOffset() - wbegin * mBaseRep->basePerWord);
  if (dir < 0) r.reverse();

  mOwnStore = false;
  mUnfolded = false; // The range might not cover the fold...
}

SeqStore::~SeqStore() {
  //cerr << "!!!!!!!!!!!!!!!!!!!!!!!!! destroy seqstor!!!!!!!!!!!!!!!!!!"  << endl;
  if (mOwnStore) {
    delete[] mStore;
  }
}


// OPERATORS ------------------------------------------------------------

// ACCESS ---------------------------------------------------------------

// METHODS --------------------------------------------------------------

void SeqStore::rebuildIndex() {
  mIndex.clear ();
  for (size_t i = 0; i < mSeqRange.size(); ++i) {
    mIndex.insert (make_pair(mSeqRange[i].end (), i));
  }
}

char SeqStore::getAminoAcid(SeqIter i, int dir) {
  SeqIter iter = i;
  count_t codon = 0;
  for (unsigned int j = 0; j < 3; ++j) {
    char b = getBase(iter, dir);
    iter += dir;
    codon <<= 2;
    switch (b) {
    case 'T': codon += 0; break;
    case 'C': codon += 1; break;
    case 'A': codon += 2; break;
    case 'G': codon += 3; break;
    default: return 'X';
    }
  }
  const char *toAA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  return toAA[codon];
}

char SeqStore::getBase(SeqIter i, int dir) const {
  return mBaseRep->getBase(mStore, i.getOffset(), dir < 0);
}

void SeqStore::setBase(SeqIter i, int dir, char c) {
  if (!mOwnStore)  cerr << "Warning: attempting to modify a SeqStore that isn't owned." << endl;

  mBaseRep->setBase(mStore, i.getOffset(), dir < 0, c);
}

void SeqStore::setBit(SeqIter i, int dir, count_t b) {
  if (!mOwnStore)  cerr << "Warning: attempting to modify a SeqStore that isn't owned." << endl;

  mBaseRep->setBit(mStore, i.getOffset(), dir < 0, b);
}

size_t SeqStore::addSequence(const char *seq_name, const char *seq) {
  size_t n = strlen (seq);
  while (mInsertPoint + n > end())
    extendStore();
  SeqRange r(mInsertPoint, mInsertPoint + n);
  size_t index = mSeqRange.size();
  mIndex.insert(make_pair(mInsertPoint + n, index));
  mSeqRange.push_back(r);
  // mSeqName.newId (seq_name);
  mNames.push_back(mSeqName.getId(seq_name));
  for (size_t i = 0; i < n; ++i) {
    setBase(mInsertPoint + i, +1, seq[i]);
  }
  mInsertPoint += n;
  return index;
}

void SeqStore::extendStore(size_t newSize) {
  // As it is, cannot be used to shrink a store
  if (newSize == 0)
    newSize = max ((size_t)512, wsize() * 2);
  seq_t *new_store = new seq_t[newSize];
  for (size_t i = 0; i < mN; ++i)
    new_store[i] = mStore[i];
  for (size_t i = mN; i < newSize; ++i)
    new_store[i] = 0;
  if (mOwnStore && mStore) delete[] mStore;
  mStore = new_store;
  mN = newSize;
  mOwnStore = true;
}

void SeqStore::clear() {
  mInsertPoint = SeqIter(0);
  mSeqName.clear();
  mSeqRange.clear();
  mNames.clear();
  mIndex.clear();
}


void SeqStore::detach() {
  // make a local copy of the store, if necessary
  if (!mOwnStore) {
    seq_t *old_store = mStore;
    mStore = new seq_t[mN];
    memcpy (mStore, old_store, mN * sizeof(seq_t));
    mOwnStore = true;
  }
}

void SeqStore::dishuffle(SeqRange r) {
  if (!mOwnStore)  cerr << "Warning: attempting to modify a SeqStore that isn't owned." << endl;

  position_t seqSize = r.size();
  char A [seqSize+1]; A[0] = '\0'; int sizeA = 0; bool A_inZ = false; int posLastEdgeA;
  char C [seqSize+1]; C[0] = '\0'; int sizeC = 0; bool C_inZ = false; int posLastEdgeC;
  char G [seqSize+1]; G[0] = '\0'; int sizeG = 0; bool G_inZ = false; int posLastEdgeG;
  char T [seqSize+1]; T[0] = '\0'; int sizeT = 0; bool T_inZ = false; int posLastEdgeT;
  char N [seqSize+1]; N[0] = '\0'; int sizeN = 0; bool N_inZ = false; int posLastEdgeN;
  char Sf = getBase(r.end () - r.dir(), r.dir());
  char S1 = getBase(r.begin(), r.dir());

  for (SeqIter i = r.begin(); i != (r.end () - r.dir()); i += r.dir()) { // ; unsigned int i = 0; i < seqSize-1; i++) {
    char c = getBase(i, r.dir());
    char c_next = getBase(i + r.dir(), r.dir());
    switch (c) {
    case 'A': A[sizeA++] = c_next; break;
    case 'C': C[sizeC++] = c_next; break;
    case 'G': G[sizeG++] = c_next; break;
    case 'T': T[sizeT++] = c_next; break;
    case 'N': N[sizeN++] = c_next; break;
    }
  }
  A[sizeA] = '\0';
  C[sizeC] = '\0';
  G[sizeG] = '\0';
  T[sizeT] = '\0';
  N[sizeN] = '\0';

  char Z [5*2+1];

  do {
    int index = 0;

    if (Sf != 'A' && strlen (A) != 0) {
      A_inZ = true; posLastEdgeA = __random(sizeA-1);
      Z[index] = 'A';
      Z[index+1] = A[posLastEdgeA];
      index += 2;
    }

    if (Sf != 'C' && strlen(C) != 0) {
      C_inZ = true; posLastEdgeC = __random(sizeC-1);
      Z[index] = 'C';
      Z[index+1] = C[posLastEdgeC];
      index += 2;
    }

    if (Sf != 'G' && strlen(G) != 0) {
      G_inZ = true; posLastEdgeG = __random(sizeG-1);
      Z[index] = 'G';
      Z[index+1] = G[posLastEdgeG];
      index += 2;
    }

    if (Sf != 'T' && strlen(T) != 0 ) {
      T_inZ = true; posLastEdgeT = __random(sizeT-1);
      Z[index] = 'T';
      Z[index+1] = T[posLastEdgeT];
      index += 2;
    }

    if (Sf != 'N' && strlen(N) != 0 ) {
      N_inZ = true; posLastEdgeN = __random(sizeN-1);
      Z[index] = 'N';
      Z[index+1] = N[posLastEdgeN];
      index += 2;
    }

    Z[index] = '\0';
  } while (!__connected(Z, Sf));

  // shuffle

  if (A_inZ)  swap (A[posLastEdgeA], A[--sizeA]);  // on place le last edge a la fin et on ajuste le sizeA pour le shuffle
  random_shuffle (A, A + sizeA);

  if (C_inZ)  swap (C[posLastEdgeC], C[--sizeC]);
  random_shuffle (C, C + sizeC);

  if (G_inZ)  swap (G[posLastEdgeG], G[--sizeG]);
  random_shuffle (G, G + sizeG);

  if (T_inZ)  swap (T[posLastEdgeT], T[--sizeT]);
  random_shuffle (T, T + sizeT);

  if (N_inZ)  swap (N[posLastEdgeN], N[--sizeN]);
  random_shuffle (N, N + sizeN);

  // Reconstruction
  setBase(r.begin(), r.dir(), S1);

  int comptA = 0;
  int comptC = 0;
  int comptG = 0;
  int comptT = 0;
  int comptN = 0;

  char c_prev = S1;
  for (SeqIter i = r.begin() + r.dir(); i != r.end(); i += r.dir()) {
    char c = '.';

    switch (c_prev) {
    case 'A':  c = A[comptA++]; break;
    case 'C':  c = C[comptC++]; break;
    case 'G':  c = G[comptG++]; break;
    case 'T':  c = T[comptT++]; break;
    case 'N':  c = N[comptN++]; break;
    }

    setBase(i, r.dir(), c);
    c_prev = c;
  }

}

void SeqStore::revcomp(SeqRange r) {
  // Only use forware ranges
  SeqIter left = r.begin();
  SeqIter right = r.end();  right--;
  while (left <= right) {
    count_t tmp = mBaseRep->complement(getBit(left, +1));
    setBit(left, +1, mBaseRep->complement(getBit(right, +1)));
    setBit(right, +1, tmp);
    left++;
    right--;
  }
}


void SeqStore::unfold() {
  if (mUnfolded)  return;  // Warning: unfolding a SeqStore that is already unfolded.
  if (!mOwnStore)  cerr << "Warning: attempting to modify a SeqStore that isn't owned." << endl;

  // Copies the reverse complement after the current sequence.
  // Beware, padding regions will be TTTT in the second half.
  seq_t *new_store = new seq_t[mN * 2];
  for (size_t i = 0; i < mN; ++i)
    new_store[i] = mStore[i];
  for (size_t i = 0; i < mN; ++i)
    new_store[mN + i] = mBaseRep->word_complement(mStore[mN - i - 1]);

  delete[] mStore;
  mStore = new_store;
  mN = 2 * mN;

  mUnfolded = true;
}

void SeqStore::randomize_n() {
  for (SeqIter i = begin(); i != end(); ++i) {
    if (getBase(i, +1) == 'N') {
      size_t index = rand() % 4;
      switch (index) {
      case 0: setBase(i, +1, 'A'); break;
      case 1: setBase(i, +1, 'C'); break;
      case 2: setBase(i, +1, 'G'); break;
      case 3: setBase(i, +1, 'T'); break;
      }
    }
  }
}

bool SeqStore::forward(SeqIter i) const {
  if (mUnfolded) {
    size_t n = size();
    size_t half = n / 2;
    position_t pos = i.getOffset();

    if (pos < half) {
      return true;
    }
    return false;
  }
  return true;
}

SeqIter SeqStore::unfold(SeqIter i) const {
  if (mUnfolded) {
    if (forward(i)) {
      return i;
    } else {
      size_t n = size();
      return SeqIter((n - 1) - i.getOffset());
    }
  } else return i;
}



// I/O  -----------------------------------------------------------------

SeqIter SeqStore::parseIter (istream &in) {
  char contig_str[4096];
  char pos_str[4096];

  while ( isspace(in.peek()) )  in.get();

  char *p = contig_str;
  while (in.peek() != ':')  *p++ = in.get();
  in.get();
  *p = 0;

  p = pos_str;
  while (isdigit(in.peek()) || in.peek() == ',') {
    *p = in.get();
    if (*p != ',') p++;
  }
  *p = 0;
  position_t pos = atol(pos_str);

  return SeqIter(getRange(contig_str), pos - 1);
}

SeqIter SeqStore::parseIter(const char *chr, size_t pos) {
  // pos is 1-based
  return SeqIter(getRange(chr), pos - 1);
}

SeqRange SeqStore::parseRangeGB(istream &in) {
  // Adhoc format for ranges: tab-sep, 1-based, inclusive, strand F/R (must be reversed)
  char buffer[4096];
  in.getline (buffer, 4096);
  char *ptr = buffer;

  char *contig_str = strsep(&ptr, ",");
  char *strand_str = strsep(&ptr, ",");
  char *start_str = strsep(&ptr, ",");
  char *end_str = strsep(&ptr, ",");

  if (end_str == NULL) return SeqRange();

  SeqRange res = SeqRange(getRange(contig_str), atoi(start_str) - 1, atoi(end_str));
  if (strand_str[0] == 'F') res.reverse();
  return res;
}

SeqRange SeqStore::parseRange(const char *chr, const char *start, const char *end, bool bed_format) {
  if (bed_format)
    return SeqRange(getRange(chr), atoi(start), atoi(end)); // already zero-base and past-the-end
  return SeqRange(getRange(chr), atoi(start) - 1, atoi(end));
}

string SeqStore::getID(SeqIter iter) const {
  //cerr << "SeqStore -> getID" << endl;
  if (mUnfolded)
    iter = unfold(iter);

  if (mIndex.upper_bound(iter) == mIndex.end())
    return "";
  unsigned int id = mIndex.upper_bound(iter)->second;

  ostringstream outf;
  outf << mSeqName.getStr(id);
  return outf.str();
}

string SeqStore::output(SeqIter iter) const {
  if (mUnfolded)
    iter = unfold(iter);

  unsigned int id = mIndex.upper_bound(iter)->second;

  string seq_name = mSeqName.getStr(id);
  size_t space = seq_name.find(" " , 0);
  seq_name = seq_name.substr(0, space);

  ostringstream outf;
  outf <<  seq_name << ':' << mSeqRange[id].getLocalOffset(iter) + 1;
  return outf.str();
}

string SeqStore::output(SeqRange range, bool tab_delim_format) const {
  if (mUnfolded)
    range = unfold(range);

  unsigned int id = mIndex.upper_bound(range.begin())->second;

  // ostringstream outf;
  char buffer[4096];

  string seq_name = mSeqName.getStr(id);
  size_t space = seq_name.find(" " , 0);
  seq_name = seq_name.substr(0, space);

  if (tab_delim_format) {
    sprintf (buffer, "%s\t%ld\t%ld",
             seq_name.c_str(),
             mSeqRange[id].getLocalOffset(range.begin()) + 1,
             mSeqRange[id].getLocalOffset(range.end()) + 1 - range.dir());
  } else {
    sprintf (buffer, "%s:%ld-%ld",
             seq_name.c_str(),
             mSeqRange[id].getLocalOffset(range.begin()) + 1,
             mSeqRange[id].getLocalOffset(range.end()) + 1 - range.dir());
  }

  return string (buffer);
}

void SeqStore::convert(SeqIter iter, const char **chr, position_t *pos) const {
  if (mUnfolded)
    iter = unfold (iter);
  unsigned int id = mIndex.upper_bound(iter)->second;
  *chr = mSeqName.getStr(id);
  *pos = mSeqRange[id].getLocalOffset(iter) + 1;
}

void SeqStore::readFasta(istream &in, bool truncate) {
  char buffer[8192];
  vector< seq_t > data;
  position_t current_pos = 0;
  position_t current_start = 0;
  seq_t cur = 0;
  size_t count = 0;

  while (!in.eof()) {
    in.getline(buffer, 8192);

    if (buffer[0] == '>' || in.eof()) {

      if (current_pos > 0) mSeqRange.push_back(SeqRange (current_start, current_pos));

      // pad and flush the current word
      // In most case it might not be necessary to keep sequences aligned...
      if (count > 0) {
        cur <<= mBaseRep->bitPerBase * (mBaseRep->basePerWord - count);
      }
      data.push_back(cur);
      current_pos += mBaseRep->basePerWord - count;
      cur = 0;
      count = 0;

      if (buffer[0] == '>') {
        current_start = current_pos;
        if (truncate) {
          char *p = &(buffer[1]);
          while (*p != '\0' && *p != ' ') p++;
          *p = '\0';
        }
        mSeqName.newId (&(buffer[1]));
      }

    } else {
      if (strlen (buffer) >= 8192) {
        cerr << strlen (buffer) << endl;
        abort ();
      }
      char *s = buffer;
      while (*s) {
        cur <<= mBaseRep->bitPerBase;
        cur += mBaseRep->baseToBit (*s);
        current_pos++;
        ++s;
        if (++count == mBaseRep->basePerWord) {
          data.push_back(cur);
          cur = 0;
          count = 0;
        }
      }

    }
  }

  mStore = new seq_t[data.size()];
  mOwnStore = true;
  mN = data.size();
  for (size_t i = 0; i < mN; ++i)
    mStore[i] = data[i];
}

void SeqStore::readBinary(const char *fn) {
  int fd = open(fn, O_RDONLY);
  struct stat sb;
  fstat(fd, &sb);  // File size

  char *addr = (char*) mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  readBinaryMM(addr);
}

void SeqStore::readBinaryOrFasta(const char *fn, bool truncate) {

  char fn_bin[4096];
  sprintf (fn_bin, "%s.bin", fn);
  ifstream test_bin(fn_bin);
  if (test_bin) {
    test_bin.close();
    readBinary(fn_bin);

  } else {
    test_bin.close();
    ifstream in(fn);
    readFasta(in, truncate);

    ofstream bin_out(fn_bin);
    if (bin_out)
      writeBinary(bin_out);

  }
  // mIndex.init (&mSeqRange);
  rebuildIndex();
}

char *SeqStore::readBinaryMM(char *base) {
  base = mSeqName.readBinary(base);
  unsigned int n = *(unsigned int *)base;
  base += sizeof(unsigned int);
  for (unsigned int i = 0; i < n; ++i) {
    SeqRange r;
    base = r.readBinary(base);
    mSeqRange.push_back(r);
  }

  mN = *(size_t*)base;  base += sizeof(size_t);
  mStore = (seq_t*)base;  base += (size_t)mN * sizeof(seq_t);
  mUnfolded = *(bool*)base;  base += sizeof(bool);
  mOwnStore = false;
  // mIndex.init (&mSeqRange);
  rebuildIndex();
  return base;
}

void SeqStore::writeBinary(ostream &out) const {
  mSeqName.writeBinary(out);
  unsigned int n = mSeqRange.size();
  out.write ((char *)&n, sizeof(unsigned int));
  for (unsigned int i = 0; i < n; ++i) {
    mSeqRange[i].writeBinary(out);
  }
  out.write ((char *)&mN, sizeof(size_t));
  long_write (out, (char *)mStore, mN * sizeof(seq_t));
  out.write ((char *)&mUnfolded, sizeof(bool));

  //cout << "Size of the saved store (B): " << (mN * sizeof(seq_t)) << endl;
}

void SeqStore::writeFastaSingle(ostream &out, SeqRange r, const char *title, unsigned int width) const {
  unsigned int count = 0;
  out << '>' << title;
  for (SeqIter i = r.begin(); i != r.end (); i += r.dir()) {
    if (count++ % width == 0) out << endl;
    out << getBase(i, r.dir());
  }
  out << endl;
}

ostream &SeqStore::output(ostream &out) const {
  for (size_t i = 0; i < mSeqName.size(); ++i)
    writeFastaSingle(out, mSeqRange[i], mSeqName.getStr(i), 200);
  return out;
}


// PRIVATE --------------------------------------------------------------

int SeqStore::__random(int n) {
  return int( float(rand()) / RAND_MAX * (n + 1) );
}

bool SeqStore::__connected(char* Z, char Sf) {
  stack< char > terminal;  terminal.push(Sf);
  int n = strlen(Z);
  int count = n / 2;

  while (!terminal.empty()) {
    char dest = terminal.top(); terminal.pop();
    for (int i = 0; i < n; i += 2) {
      if (Z[i+1] == dest) {
        terminal.push(Z[i]);
        Z[i+1] = 'X';
        count--;
      }
    }
  }
  return count == 0;
}

ostream &operator<< (ostream &out, const SeqStore &o) {
  return o.output(out);
}
