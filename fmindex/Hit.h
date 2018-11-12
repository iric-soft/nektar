//                              -*- Mode: C++ -*-
// Hit.h
//

#ifndef __Hit_h_
#define __Hit_h_

#include <iostream>
#include <algorithm>

using namespace std;

#include "types.h"
#include "SeqRange.h"
#include "SeqIter.h"


class Hit
{

  size_t num_diag;
  size_t pos_i;
  size_t pos_j;
  SeqRange seq_r;

public:

  // LIFECYCLE ------------------------------------------------------------

  Hit() {}
  Hit(const size_t d, const size_t pi, const size_t pj, const SeqRange range);
  Hit(const Hit& h);

  // OPERATORS ------------------------------------------------------------
  bool operator== (const Hit &h) const ;

  // ACCESS ---------------------------------------------------------------

  SeqRange seq_range() const { return seq_r;}
  size_t diag() const { return num_diag;}
  size_t i() const { return pos_i;}
  size_t j() const { return pos_j;}
  SeqIter range_begin() const { return seq_r.begin(); }
  SeqIter range_end() const { return seq_r.end(); }
  size_t range_size() const { return seq_r.size(); }

  // METHODS --------------------------------------------------------------

  void range_extend(size_t extension, const SeqStore* refStore);
  bool range_overlap(Hit& h) const;
  void merge(const Hit& h);

  // I/O  -----------------------------------------------------------------

  ostream &output(ostream &out) const;

};

ostream &operator<< (ostream &out, const Hit &h);

#endif
