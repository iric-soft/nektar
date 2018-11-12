//                              -*- Mode: C++ -*-
// SeqRange.cc
// Copyright 2010 Laboratoire de Bioinformatique Fonctionnelle et Structurale
//                Universite de Montreal.
// Author           : Sebastien Lemieux <s.lemieux@umontreal.ca>
// Created On       : Wed Nov 10 14:35:51 2010
//

#include "SeqRange.h"
#include "SeqStore.h"

// LIFECYCLE ------------------------------------------------------------

SeqRange::SeqRange ()
    : mBegin (0), mEnd (0) {}

SeqRange::SeqRange (position_t begin, position_t end)
  : mBegin (begin), mEnd (end) {}

SeqRange::SeqRange (SeqIter begin, SeqIter end)
  : mBegin (begin.getOffset ()), mEnd (end.getOffset ()) {}

SeqRange::SeqRange (const SeqRange &parent, position_t begin, position_t end) {
  mBegin = parent.mBegin + parent.dir() * begin;
  mEnd = parent.mBegin + parent.dir() * end;
}

SeqRange::SeqRange (const SeqRange &o)
  : mBegin (o.mBegin), mEnd (o.mEnd) {}


// OPERATORS ------------------------------------------------------------

SeqRange &SeqRange::operator= (const SeqRange &o)
{
  if (this != &o) {
    mBegin = o.mBegin;
    mEnd = o.mEnd;
  }
  return *this;
}

bool SeqRange::operator< (const SeqRange &o) const {
  if (dir() < o.dir()) return true;
  if (dir() > o.dir()) return false;
  if (dir() > 0) {
    return (mBegin < o.mBegin);
  } else {
    return (mBegin > o.mBegin);
  }
}

bool SeqRange::operator== (const SeqRange &o) const {
  return (mBegin == o.mBegin) && (mEnd == o.mEnd);
}

// ACCESS ---------------------------------------------------------------

// METHODS --------------------------------------------------------------

void SeqRange::merge(const SeqRange &r) {
  if (dir() != r.dir()) {
    cerr << "SeqRange::merge: incompatible directions." << endl;
    cerr << *this << endl;
    cerr << r << endl;
  }
  mBegin = (dir() > 0)?min (mBegin, r.mBegin):max (mBegin, r.mBegin);
  mEnd = (dir() > 0)?max (mEnd, r.mEnd):min (mEnd, r.mEnd);
}

void SeqRange::extend(position_t left, position_t right, const SeqStore* refStore) {
  //cerr << "SeqRange -> extend" << endl;
  string id = refStore->getID(mBegin);
  position_t extend_begin = mBegin;
  position_t extend_end = mEnd;

  if (extend_begin < left) extend_begin = 0;
  else extend_begin -= left;
  extend_end += right;

  while ( id.compare(refStore->getID(extend_begin))!= 0 ) {
    --left;
    extend_begin = mBegin;
    if (extend_begin < left) extend_begin = 0;
    else extend_begin -= left;
  }

  while ( id.compare(refStore->getID(extend_end))!= 0 ) {
    --right;
    extend_end = mEnd;
    extend_end += right;
  }

  mBegin = extend_begin;
  mEnd = extend_end;
}

void SeqRange::clip(const SeqRange &r) {
  if (mBegin < r.mBegin) mBegin = r.mBegin;
  if (mEnd > r.mEnd) mEnd = r.mEnd;
}

bool SeqRange::overlap(const SeqRange &r) const {  // includes juxtaposition
  if (dir() != r.dir()) {
    cerr << "SeqRange::overlap: incompatible directions." << endl;
    cerr << *this << endl;
    cerr << r << endl;
  }
  if (mBegin > r.mBegin)
    if (r.mEnd >= mBegin) return true;
    else return false;
  else
    if (r.mBegin <= mEnd) return true;
    else return false;
}

bool SeqRange::contains(SeqIter i) const {
  position_t pos = i.getOffset();
  if (dir() > 0) {
    return pos >= mBegin && pos < mEnd;
  } else {
    return pos <= mBegin && pos > mEnd;
  }
}

position_t SeqRange::getLocalOffset(SeqIter i) const {
  // Assumes that i shares the same SeqStore
  // Assumes that i is after mBegin (return value is unsigned)
  position_t res = i.getOffset();
  return (dir() > 0)?(res - mBegin):(mBegin - res);
}

// I/O  -----------------------------------------------------------------

char *SeqRange::readBinary (char *base) {
  mBegin = *(position_t *)base;  base += sizeof (position_t);
  mEnd = *(position_t *)base;    base += sizeof (position_t);
  return base;
}

void SeqRange::writeBinary (ostream &out) const {
  out.write ((char *)&mBegin, sizeof (position_t));
  out.write ((char *)&mEnd, sizeof (position_t));
}

ostream &SeqRange::output (ostream &out) const {
  out << "(" << mBegin << ", " << mEnd << ")";
  return out;
}

ostream &operator<< (ostream &out, const SeqRange &o) {
  return o.output (out);
}
