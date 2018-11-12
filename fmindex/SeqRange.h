//                              -*- Mode: C++ -*-
// SeqRange.h
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale
//                Universite de Montreal.
// Author           : Sebastien Lemieux <s.lemieux@umontreal.ca>
// Created On       : Wed Nov 10 14:35:51 2010
//

#ifndef __SeqRange_h_
#define __SeqRange_h_

#include <iostream>

using namespace std;

#include "types.h"
#include "SeqIter.h"

class SeqStore; //#include "SeqStore.h"

class SeqRange
{

  position_t mBegin;
  position_t mEnd;

public:

  // LIFECYCLE ------------------------------------------------------------

  SeqRange ();
  SeqRange (position_t begin, position_t end);
  SeqRange (SeqIter begin, SeqIter end);
  SeqRange (const SeqRange &parent, position_t begin, position_t end);
  SeqRange (const SeqRange &o);
  ~SeqRange () {};

  // OPERATORS ------------------------------------------------------------

  SeqRange &operator= (const SeqRange &o);
  bool operator< (const SeqRange &o) const ;
  bool operator== (const SeqRange &o) const ;
  // Although valid, a range with begin: 0, end: 0 is considered null
  operator bool () const { return mBegin != 0 || mEnd != 0; }

  // ACCESS ---------------------------------------------------------------

  //TODO avoid cast position_t and seqIter
  void setBegin (const SeqIter begin) { mBegin = begin.getOffset(); }
  void setEnd (const SeqIter end) { mEnd = end.getOffset(); }
  SeqIter begin () { return SeqIter (mBegin); }
  SeqIter end () { return SeqIter (mEnd); }
  position_t getBegin () const { return mBegin; }
  position_t getEnd () const { return mEnd; }
  const SeqIter begin () const { return SeqIter (mBegin); }
  const SeqIter end () const { return SeqIter (mEnd); }
  position_t size () const { return (mBegin > mEnd)?(mBegin - mEnd):(mEnd - mBegin); }
  int dir () const { return (mBegin > mEnd)?(-1):(+1); }
  inline bool empty() const {return (mBegin==0 && mEnd ==0); }

  // METHODS --------------------------------------------------------------

  void reverse () { mEnd -= dir (); swap (mBegin, mEnd); mEnd += dir(); }
  void merge (const SeqRange &r);
  void extend (position_t left, position_t right, const SeqStore* refStore);  // assumes both +1 dir
  void clip (const SeqRange &r);         // assumes both +1 dir
  bool overlap (const SeqRange &r) const ;  // includes juxtaposition
  bool contains (SeqIter i) const ;


  position_t getLocalOffset (SeqIter i) const ;
  SeqIter getIter (int diff) { return SeqIter ((dir () > 0)?(mBegin + diff):(mBegin - diff)); }

  // I/O  -----------------------------------------------------------------

  char *readBinary (char *base);
  void writeBinary (ostream &out) const ;

  ostream &output (ostream &out) const ;
};

ostream &operator<< (ostream &out, const SeqRange &o);

#endif
