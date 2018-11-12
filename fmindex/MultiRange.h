//                              -*- Mode: C++ -*-
// MultiRange.h
//

#ifndef __MultiRange_h_
#define __MultiRange_h_

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

#include "SeqRange.h"

class SeqStore; //#include "SeqStore.h"

class MultiRange
{

  typedef set< SeqRange > set_t;
  set_t mRange;

public:

  struct iterator {
    MultiRange *parent;
    SeqIter i;
    set_t::iterator r;
    bool operator!= (const MultiRange::iterator &o);
    bool operator== (const MultiRange::iterator &o);
    MultiRange::iterator &operator++ ();
    MultiRange::iterator &operator-- ();
    operator SeqIter() { return i; }
  };

  typedef set_t::const_iterator range_iterator;

  // LIFECYCLE ------------------------------------------------------------

  MultiRange() {}
  MultiRange(SeqRange r) { mRange.insert(r); }
  MultiRange(const MultiRange &o);

  // OPERATORS ------------------------------------------------------------

  MultiRange &operator= (const MultiRange &o);
  operator bool() const { return !empty(); }

  // ACCESS ---------------------------------------------------------------

  iterator begin();
  iterator end();
  iterator rbegin();
  iterator rend();
  size_t size() const ;  // This cannot be used to iterate
  range_iterator range_begin() const { return mRange.begin(); }
  range_iterator range_end() const { return mRange.end(); }
  size_t range_size() const { return mRange.size(); }
  bool empty() const { return mRange.empty(); }

  iterator find(const SeqIter pos);

  // METHODS --------------------------------------------------------------

  void insert(const SeqRange o) { mRange.insert(o); }
  void insert(const MultiRange o);
  void extend_merge(size_t extension, const SeqStore* refStore);  // also merges

  // I/O  -----------------------------------------------------------------

  ostream &output(ostream &out) const;

};

ostream &operator<< (ostream &out, const MultiRange &o);


#endif
