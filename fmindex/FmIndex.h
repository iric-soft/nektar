//                              -*- Mode: C++ -*-
// FmIndex.h
//

#ifndef __FmIndex_h_
#define __FmIndex_h_

// Includes for mmap stuff...
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>


#include "util.h"
#include "SeqStore.h"
#include "SeqIter.h"
#include "BaseRep.h"

class FmIndex
{

  BaseRep *mBaseRep;
  SeqStore mSeq;
  position_t mSeqLength; // size of the store in bp

  count_t mBucketSize;
  count_t mNbBucket;
  count_t mMarkFreq;
  count_t mNbMark;
  count_t mCacheSize; // as 4^n
  seq_t *mL;
  count_t *mOcc;  // These counts should not exceed 4G (32-bit)
  position_t *mC;
  position_t *mPos;
  position_t *mO;

  bool mNew;

public:

  // LIFECYCLE ------------------------------------------------------------

  FmIndex ();
  FmIndex (const char *base);
  FmIndex (const char *genome_fn, const char *base, const bool rand=true, const bool reverse=true);
  ~FmIndex ();

  // OPERATORS ------------------------------------------------------------

  // ACCESS ---------------------------------------------------------------

  SeqStore &Seq () { return mSeq; }
  //const BaseRep &getBaseRep () const { return mBaseRep; }
  //Tariq
  const BaseRep &getBaseRep () const { return *mBaseRep; }
  SeqStore *getStore () { return &mSeq; }

  // METHODS --------------------------------------------------------------

  void CreateIndex (unsigned int nBucketSize, unsigned int nMarkFreq);
  void CreateSortedArray ();
  void RebuildOcc (count_t bucketSize);

  position_t Occ (count_t c_bit, count_t c_index, position_t pos) const ;
  position_t Locate (position_t i) const ;

  void InitLookup (count_t c_bit, position_t *out_sp, position_t *out_ep) const ;
  void Lookup (position_t cur_sp, position_t cur_ep, count_t c_bit, position_t *out_sp, position_t *out_ep) const ;
  const char* SearchSeq (const char *seq, int max_match=50) const;
  const char* CheckSeq (const char *seq) const;


  // I/O  -----------------------------------------------------------------

  void Save (const char *fn);
  void Load (const char *fn);
  string GetLocateChr(position_t i) const;
  bool GetLocateStrand(position_t i) const;

  // Private --------------------------------------------------------------

private:

  struct comp_seq {
    // Compare two suffixes, the sequence is considered circular

    position_t begin, end;
    seq_t *base;
    const BaseRep *rep;

    comp_seq (SeqStore *store_ptr);
    bool operator() (position_t a, position_t b);
  };
  void DeleteArray();

};

#endif
