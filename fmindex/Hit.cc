//                              -*- Mode: C++ -*-
// Hit.cc
//

#include "Hit.h"

// LIFECYCLE ------------------------------------------------------------

Hit::Hit(const size_t d, const size_t pi, const size_t pj, const SeqRange range)
  : num_diag(d), pos_i(pi), pos_j(pj), seq_r(range) {}

Hit::Hit(const Hit &h)
  : num_diag(h.num_diag), pos_i(h.pos_i), pos_j(h.pos_j), seq_r(h.seq_r) {}


// OPERATORS ------------------------------------------------------------
bool Hit::operator== (const Hit &h) const {
  return this->seq_r.begin() == h.seq_range().begin() && 
  		 this->seq_r.end() == h.seq_range().end() &&
  		 this->pos_i == h.i() && this->pos_j == h.j();
}


// ACCESS ---------------------------------------------------------------

// METHODS --------------------------------------------------------------

void Hit::range_extend(size_t extension, const SeqStore* refStore) {
	//cerr << "Hit -> range_extend" << endl;
  this->seq_r.extend(extension, extension, refStore);
}

bool Hit::range_overlap(Hit& h) const {
	return this->seq_r.overlap(h.seq_range());
}

void Hit::merge(const Hit& h) {
  if (this->seq_r.dir() != h.seq_range().dir()) {
    cerr << "SeqRange::merge: incompatible directions." << endl;
    cerr << this->seq_r << endl;
    cerr << h.seq_range() << endl;
  }
  if (this->seq_r.dir() > 0) {
  	this->seq_r.setBegin(min(this->seq_r.begin(), h.seq_range().begin()));
  	this->seq_r.setEnd(max(this->seq_r.end(), h.seq_range().end()));
  } else {
  	this->seq_r.setBegin(max(this->seq_r.begin(), h.seq_range().begin()));
  	this->seq_r.setEnd(min(this->seq_r.end(), h.seq_range().end()));
  }

  this->pos_i = min(this->pos_i, h.i());
  this->pos_j = max(this->pos_j, h.j());

  this->num_diag = min(this->num_diag, h.diag());
}

// I/O  -----------------------------------------------------------------

ostream &Hit::output(ostream &out) const {
  out << "diag: " << num_diag;
  out << ", i: " << pos_i;
  out << ", j: " << pos_j;
  out << ", range: " << seq_r << endl;
  return out;
}

ostream &operator<< (ostream &out, const Hit &h) {
  return h.output(out);
}
