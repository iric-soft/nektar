//                              -*- Mode: C++ -*-
// SeqStore.h

#ifndef __SeqStore_h_
#define __SeqStore_h_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <stack>
#include <sstream>

#include <sys/mman.h>  // for memory-mapped files
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

class SeqStore;

#include "types.h"
#include "util.h"
#include "SeqIter.h"
#include "SeqRange.h"
#include "MultiRange.h"
#include "BaseRep.h"
#include "StrDict.h"

class SeqStore {

  BaseRep *mBaseRep;  // Base representation object

  seq_t *mStore;
  size_t mN;     // Number of seq_t in mSeq

  StrDict mSeqName;
  vector< SeqRange > mSeqRange;
  vector< unsigned int > mNames;  // This keeps the names even when duplicate names are present...
  map< SeqIter, size_t > mIndex;


  bool mOwnStore;
  bool mUnfolded;
  SeqIter mInsertPoint;

public:

  // LIFECYCLE ------------------------------------------------------------

  SeqStore(BaseRep *rep);
  SeqStore(SeqStore &o, SeqRange &r);
  SeqStore(BaseRep *rep, size_t n);
  ~SeqStore();

  // OPERATORS ------------------------------------------------------------

  // ACCESS ---------------------------------------------------------------

  size_t size() const { return mN * mBaseRep->basePerWord; }
  size_t wsize() const  { return mN; }
  SeqIter begin() const { return SeqIter(0); }
  SeqIter end() const { return SeqIter(size()); }

  // METHODS --------------------------------------------------------------

  void rebuildIndex();

  count_t getBit(SeqIter i, int dir = +1) const { return mBaseRep->getBit(mStore, i.getOffset(), dir < 0); }
  char getBase(SeqIter i, int dir = +1) const ;
  void setBase(SeqIter i, int dir, char c);
  void setBit(SeqIter i, int dir, count_t b);

  // watch your step!
  seq_t *getStore() { return mStore; }
  const BaseRep *getBaseRep() const { return mBaseRep; }

  size_t addSequence(const char *seq_name, const char *seq);
  void extendStore(size_t newSize = 0);
  void clear();

  const StrDict &names() const { return mSeqName; }
  vector< SeqRange > &ranges() { return mSeqRange; }

  char getAminoAcid(SeqIter i, int dir = +1);
  SeqRange getRange(size_t i) const { return mSeqRange[i]; }
  SeqRange getRange(const char *seq_name) { return getRange(mSeqName.getId(seq_name)); }
  bool hasRange(const char *seq_name) { return mSeqName.hasId(seq_name); }
  const char *getName(size_t i) const { return mSeqName.getStr(mNames[i]); }
  void detach();

  void dishuffle(SeqRange r);
  void revcomp(SeqRange r);
  void unfold();
  void randomize_n();

  bool forward(SeqIter i) const ;
  SeqIter unfold(SeqIter i) const ;
  SeqRange unfold(SeqRange i) const { return SeqRange(unfold(i.begin()), unfold(i.end())); }
  MultiRange unfold(MultiRange i) const ;

  // I/O  -----------------------------------------------------------------

  SeqIter parseIter(istream &in);
  SeqIter parseIter(const char *chr, size_t pos); // pos is 1-based
  SeqRange parseRangeGB(istream &in);
  SeqRange parseRange(const char *chr, const char *start, const char *end, bool bed_format = false);

  string getID(SeqIter iter) const ;
  string output(SeqIter iter) const ;
  string output(SeqRange range, bool tab_delim_format = false) const ;

  void convert(SeqIter iter, const char **chr, position_t *pos) const ;

  void readFasta(istream &in, bool truncate = false);
  void readBinary(const char *fn);
  void readBinaryOrFasta(const char *fn, bool truncate = false);
  char *readBinaryMM(char *base);
  void writeBinary(ostream &out) const ;

  void writeFastaSingle(ostream &out, SeqRange r, const char *title, unsigned int width = 80) const ;

  ostream &output(ostream &out) const;

  // PRIVATE --------------------------------------------------------------

private:
  int __random(int n);
  bool __connected(char* Z, char Sf);

};

ostream &operator<<(ostream &out, const SeqStore &o);

#endif
