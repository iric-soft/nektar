//                              -*- Mode: C++ -*-
// SeqIter.h
// Copyright 2009 Laboratoire de Bioinformatique Fonctionnelle et Structurale
//                Universite de Montreal.
// Author           : Sebastien Lemieux <s.lemieux@umontreal.ca>
// Created On       : Mon Oct  4 12:55:07 2010
//

#ifndef __SeqIter_h_
#define __SeqIter_h_

#include <iostream>

using namespace std;

#include "types.h"

class SeqRange;

class SeqIter
{

  position_t mOffset;

public:

  // LIFECYCLE ------------------------------------------------------------

  SeqIter () : mOffset (0) {}
  SeqIter (position_t offset)
    : mOffset (offset) {}
  SeqIter (const SeqRange &parent, position_t p);

  // OPERATORS ------------------------------------------------------------

  SeqIter &operator= (const SeqIter &o);
  bool operator< (const SeqIter &o) const { return mOffset < o.mOffset; }
  bool operator> (const SeqIter &o) const { return mOffset > o.mOffset; }
  bool operator>= (const SeqIter &o) const { return mOffset >= o.mOffset; }
  bool operator<= (const SeqIter &o) const { return mOffset <= o.mOffset; }

  SeqIter &operator++ () { mOffset++; return *this;}
  SeqIter &operator++ (int) { mOffset++; return *this;}
  SeqIter &operator-- () { mOffset--; return *this;}
  SeqIter &operator-- (int) { mOffset--; return *this;}
  SeqIter &operator+= (int i) { mOffset += i; return *this; }
  SeqIter &operator-= (int i) { mOffset -= i; return *this; }
  bool operator== (const SeqIter &o) const ;
  bool operator!= (const SeqIter &o) const ;
  position_t operator- (const SeqIter &o) const ;
  SeqIter operator+ (int i) const ;
  SeqIter operator- (int i) const ;

  // ACCESS ---------------------------------------------------------------

  position_t getOffset () const { return mOffset; }

  // METHODS --------------------------------------------------------------

  // I/O  -----------------------------------------------------------------

  ostream &output (ostream &out) const;

};

ostream &operator<< (ostream &out, const SeqIter &o);

#endif
