//                              -*- Mode: C++ -*-
// Hits.h
//

#ifndef __Hits_h_
#define __Hits_h_

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>

using namespace std;

#include "./types.h"
#include "Hit.h"
#include "SeqRange.h"
#include "SeqIter.h"
#include "MultiRange.h"


class Hits {
  set<size_t> all_diags;
  vector<Hit> all_hits;
  double** adj_matrix;
  size_t len_matrix;

 public:
  // LIFECYCLE ------------------------------------------------------------

  Hits();
  Hits(const Hits &h);
  ~Hits();

  // ACCESS ---------------------------------------------------------------
  inline size_t size() const { return all_hits.size(); }
  inline const vector<Hit>& get_hits() const { return all_hits; }
  MultiRange get_genome_range() const;

  // METHODS --------------------------------------------------------------

  void add_hit(const size_t& d, const size_t& i, const size_t& j, const SeqRange& seq_r);
  
  inline void sort_on_range() {sort(all_hits.begin(), all_hits.end(), [](const Hit& first, const Hit& sec){return first.seq_range() < sec.seq_range();});}
  void sort_on_pos();
  
  void extend_merge(size_t extension, const SeqStore* refStore);
    
  void init_adj_matrix(const size_t& len_seq, const int max_overlap=31);
  void run_warshall(const size_t& len_seq, const int max_overlap=31);

  string matrix_to_string() const;
  void delete_matrix();

  // I/O  -----------------------------------------------------------------
  ostream &output(ostream &out) const;
};

ostream &operator<< (ostream &out, const Hits &h);

#endif //__Hits_h_
