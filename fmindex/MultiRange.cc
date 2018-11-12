//                              -*- Mode: C++ -*-
// MultiRange.cc
//

#include "MultiRange.h"
#include "SeqStore.h"


bool MultiRange::iterator::operator!= (const iterator &o) {
  return (r != o.r) || (i != o.i);
}

bool MultiRange::iterator::operator== (const iterator &o) {
  return (r == o.r) && (i == o.i);
}

MultiRange::iterator &MultiRange::iterator::operator++() {
  i += r->dir();
  if (i == r->end()) {
    ++r;
    if (r == parent->mRange.end()) {
      i = parent->mRange.rbegin()->end();
    } else {
      i = r->begin();
    }
  }
  return *this;
}

MultiRange::iterator &MultiRange::iterator::operator--() {
  if (*this == parent->end()) {
    *this = parent->rbegin();
  } else {
    if (r != parent->mRange.begin() && i == r->begin()) {
      --r;
      i = r->end();
    }
    i-= r->dir();
  }

  return *this;
}


// LIFECYCLE ------------------------------------------------------------

MultiRange::MultiRange(const MultiRange &o)
  : mRange(o.mRange) {}


// OPERATORS ------------------------------------------------------------

MultiRange &MultiRange::operator= (const MultiRange &o)
{
  if (this != &o) {
    mRange = o.mRange;
  }
  return *this;
}

// ACCESS ---------------------------------------------------------------

MultiRange::iterator MultiRange::begin() {
  iterator i;
  i.parent = this;
  i.r = mRange.begin();
  i.i = i.r->begin();
  return i;
}

MultiRange::iterator MultiRange::end() {
  iterator i;
  i.parent = this;
  i.r = mRange.end();
  i.i = mRange.rbegin()->end ();
  return i;
}

MultiRange::iterator MultiRange::rbegin() {
  iterator i;
  i.parent = this;
  i.r = mRange.end();
  --i.r;
  i.i = i.r->end();
  i.i -= i.r->dir();
  return i;
}

MultiRange::iterator MultiRange::rend() {
  iterator i;
  i.parent = this;
  i.r = mRange.begin();
  i.i = i.r->begin();
  i.i -= i.r->dir();
  return i;
}

size_t MultiRange::size() const {
  size_t s = 0;
  for (range_iterator i = range_begin(); i != range_end(); ++i)
    s += i->size();
  return s;
}

MultiRange::iterator MultiRange::find(const SeqIter pos)
{
  // This could be optimized by skipping ranges
  for (iterator i=this->begin(); i!=this->end(); ++i) {
    if (i.i == pos) {
      return i;
    }
  }
  return this->end();
}

// METHODS --------------------------------------------------------------

void MultiRange::insert(const MultiRange o) {
  for(MultiRange::range_iterator i = o.range_begin(); i != o.range_end(); ++i)
    insert(*i);
}

void MultiRange::extend_merge(size_t extension, const SeqStore* refStore) {
  if (mRange.empty()) return; // Nothing to do

  set_t new_range;
  SeqRange last_range;
  for (set_t::iterator i = mRange.begin(); i != mRange.end(); ++i) {
    SeqRange tmp = *i;
    tmp.extend(extension, extension, refStore);
    if (i == mRange.begin()) {
      last_range = tmp;
    } else {
      if (last_range.overlap(tmp)) {
	       last_range.merge(tmp);
      } else {
	       new_range.insert(last_range);
         last_range = tmp;
      }
    }
  }
  new_range.insert(last_range);
  mRange = new_range;
}


// I/O  -----------------------------------------------------------------

ostream &MultiRange::output(ostream &out) const {
  for (set_t::const_iterator i = mRange.begin(); i != mRange.end(); ++i) {
    out << *i << i->dir() << endl;
  }

  return out;
}

ostream &operator<< (ostream &out, const MultiRange &o) {
  return o.output(out);
}
