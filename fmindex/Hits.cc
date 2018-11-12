//                              -*- Mode: C++ -*-
// Hits.cc
//

#include "Hits.h"

// LIFECYCLE ------------------------------------------------------------

Hits::Hits()
  : all_diags(set<size_t>()), all_hits(vector<Hit>()), adj_matrix(NULL), len_matrix(0) {}

Hits::Hits(const Hits &h)
  : all_diags(h.all_diags), all_hits(h.all_hits), adj_matrix(NULL), len_matrix(0) {}

Hits::~Hits() {
  //cerr << "!!!!!!!!!!!!!!!!!!!!!!!!! destroy seqstor!!!!!!!!!!!!!!!!!!"  << endl;
  this->delete_matrix();
}

// OPERATORS ------------------------------------------------------------

// ACCESS ---------------------------------------------------------------

MultiRange Hits::get_genome_range() const {
  MultiRange res;

  for (auto & hit : all_hits) {
    res.insert(SeqRange(hit.range_begin(),hit.range_end()));
  }

  return res;
}

// METHODS --------------------------------------------------------------

void Hits::add_hit(const size_t& d, const size_t& i, const size_t& j,
                   const SeqRange& seq_r) {
  all_diags.insert(d);
  all_hits.push_back(Hit(d, i, j, seq_r));
}

void Hits::sort_on_pos() {
  sort(all_hits.begin(), all_hits.end(), [](const Hit& first, const Hit& sec) {
    if (first.i() == sec.i()) {
      if (first.j() == sec.j())
        return first.seq_range() < sec.seq_range();
      else
        return first.j() < sec.j();
    }

    return first.i() < sec.i();
  });
}

void Hits::extend_merge(size_t extension, const SeqStore* refStore) {
  //cerr << "hits -> extend_merge" << endl;
  if (all_hits.empty()) return;  // Nothing to do

  this->sort_on_range();

  vector<Hit> new_hits;
  Hit last_hit;
  for (auto & hit : all_hits) {
    Hit tmp_hit = hit;
    tmp_hit.range_extend(extension, refStore);
    if (hit == *all_hits.begin()) {
      last_hit = hit;
    } else {
      if (last_hit.range_overlap(tmp_hit)) {
         last_hit.merge(tmp_hit);
      } else {
         new_hits.push_back(last_hit);
         last_hit = tmp_hit;
      }
    }
  }
  new_hits.push_back(last_hit);
  all_hits = new_hits;
}

void Hits::init_adj_matrix(const size_t& len_seq, const int max_overlap) {
  //cerr << "init_adj_matrix" << endl;
  if (all_hits.empty()) return;  // Nothing to do
  this->sort_on_pos();

  double infinity = numeric_limits<double>::max();
  len_matrix = all_hits.size() + 2;
  double diff = 0;

  adj_matrix = new double* [len_matrix];

  for (size_t x = 0; x < len_matrix; ++x) {
    adj_matrix[x] = new double[len_matrix];
    for (size_t y = 0; y < len_matrix; ++y) {
      if (x == 0) {
        if (y == len_matrix - 1)
          adj_matrix[x][y] = static_cast<double>(len_seq);
        else if (y > 0)
          adj_matrix[x][y] = all_hits[y-1].i();
        else
          adj_matrix[x][y] = infinity;
      } else {
        if (y == len_matrix - 1) {
          if (x ==len_matrix - 1)
            adj_matrix[x][y] = infinity;  
          else 
            adj_matrix[x][y] = static_cast<double>(len_seq-1) - all_hits[x-1].j();
        } else {
          if (y <= x)
            adj_matrix[x][y] = infinity;
          else {
            diff = static_cast<double>(all_hits[y-1].i()) - static_cast<double>(all_hits[x-1].j());
            if (all_hits[y-1].j() > all_hits[x-1].j() && diff > -max_overlap)
              adj_matrix[x][y] = max(0.0, diff);
            else 
              adj_matrix[x][y] = infinity;
          }
        }
      }
    }
  }
  //cerr << this->matrix_to_string() << endl;
}

void Hits::run_warshall(const size_t& len_seq, const int max_overlap) {
  //cerr << "run_warshall" << endl;
  if (all_hits.empty()) return;  // Nothing to do
  this->delete_matrix();
  this->init_adj_matrix(len_seq, max_overlap);
  //cerr << "end init" << endl;
  //cerr << matrix_to_string();

  double infinity = numeric_limits<double>::max();
  double cost_change_hit = 0.0;
  len_matrix = all_hits.size() + 2;
  // cerr << "size\t" << len_matrix << endl;

  for (size_t k = 0; k < len_matrix; ++k) {
    for (size_t i = 0; i < len_matrix; ++i) {
      for (size_t j = k; j < len_matrix; ++j) {
        if (adj_matrix[i][k] != infinity && adj_matrix[k][j] != infinity &&
            adj_matrix[i][j] > (adj_matrix[i][k] + adj_matrix[k][j] + cost_change_hit))
          adj_matrix[i][j] = adj_matrix[i][k] + adj_matrix[k][j] + cost_change_hit;
      }
    }
  }

  double best_score = adj_matrix[0][len_matrix -1];
  //cerr << "best score: " << best_score << endl;
  vector<Hit> new_hits;
  for (size_t i = 1; i < len_matrix - 1; ++i) {
    //cerr << i << " score: " << adj_matrix[i][len_matrix - 1] + adj_matrix[0][i] + 0.01 << endl;
    if (adj_matrix[i][len_matrix - 1] + adj_matrix[0][i] + cost_change_hit == best_score) {
      new_hits.push_back(all_hits[i-1]);
    }
  }

  all_hits = new_hits;
}

string Hits::matrix_to_string() const {
  string res = "Adjacente matrix:\n" ;

  if (adj_matrix != NULL) {

    for (size_t i = 0; i < len_matrix; ++i) {
      if (i == 0)
        res += "   S       ";
      else if (i == len_matrix - 1)
        res += "    P       ";
      else if (i < 10)
        res += "    " + to_string(i) + "       ";
      else if (i < 100)
        res += "   " + to_string(i) + "       ";
      else if (i < 1000)
        res += "  " + to_string(i) + "       ";
      else
        res += " " + to_string(i) + "       ";
    }
    res += "\n";

    for (size_t i = 0; i < len_matrix; ++i) {
      for (size_t j = 0; j < len_matrix; ++j) {
        if (adj_matrix[i][j] > 10000) {
          res += "   *        ";
        } else if (adj_matrix[i][j] < 10) {
          res += "   " + to_string(adj_matrix[i][j]) + " ";
        } else if (adj_matrix[i][j] < 100) {
          res += "  " + to_string(adj_matrix[i][j]) + " ";
        } else if (adj_matrix[i][j] < 1000) {
          res += " " + to_string(adj_matrix[i][j]) + " ";
        } else {
          res += to_string(adj_matrix[i][j]) + " ";
        }
      }
      res += "\n";
    }
  }

  return res;
}

void Hits::delete_matrix() {
  if (adj_matrix != NULL) {
    for (size_t i=0; i < len_matrix; ++i) {
      delete[] adj_matrix[i];
    }
    delete[] adj_matrix;
  }
}


// I/O  -----------------------------------------------------------------

ostream &Hits::output(ostream &out) const {
  out << this->matrix_to_string();
  int cpt = 1;
  for (auto & hit : all_hits) {
    out << "\t" << cpt << ") " << hit;
    ++cpt;
  }

  return out;
}

ostream &operator<< (ostream &out, const Hits &h) {
  return h.output(out);
}
