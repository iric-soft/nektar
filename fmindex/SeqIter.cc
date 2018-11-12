//                              -*- Mode: C++ -*-
// SeqIter.cc
// Copyright 2009 Laboratoire de Bioinformatique Fonctionnelle et Structurale
//                Universite de Montreal.
// Author           : Sebastien Lemieux <s.lemieux@umontreal.ca>
// Created On       : Mon Oct  4 12:55:07 2010
//

#include "SeqIter.h"
#include "SeqRange.h"

// LIFECYCLE ------------------------------------------------------------

SeqIter::SeqIter (const SeqRange &parent, position_t p) {
  mOffset = parent.begin().mOffset + parent.dir() * p;
}

// OPERATORS ------------------------------------------------------------

SeqIter &SeqIter::operator= (const SeqIter &o) {
  if (this != &o) {
    mOffset = o.mOffset;
  }
  return *this;
}

bool SeqIter::operator== (const SeqIter &o) const {
  return (mOffset == o.mOffset);
}

bool SeqIter::operator!= (const SeqIter &o) const {
  return (mOffset != o.mOffset);
}

position_t SeqIter::operator- (const SeqIter &o) const {
  return mOffset - o.mOffset;
}

SeqIter SeqIter::operator+ (int i) const {
  SeqIter res;
  res.mOffset = mOffset + i;
  return res;
}

SeqIter SeqIter::operator- (int i) const {
  SeqIter res;
  res.mOffset = mOffset - i;
  return res;
}

// ACCESS ---------------------------------------------------------------

// METHODS --------------------------------------------------------------

// I/O  -----------------------------------------------------------------

ostream &SeqIter::output (ostream &out) const
{
  out << "(pos: " << mOffset << ")";
  return out;
}

ostream &operator<< (ostream &out, const SeqIter &o)
{
  return o.output (out);
}

