#ifndef NEKTAR_CPP_UTIL_KMER_H_
#define NEKTAR_CPP_UTIL_KMER_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <sstream>

// using std::map;
// using std::vector;
// using std::string;
// using std::stringstream;
// using std::getline;

using namespace std;

namespace util {

class Kmer {
  string s_seq;
  int i_group;
  int i_vertex;
  int i_dbsnp;
  vector<int> v_count;
  map<string, vector<double> > m_score;

 public:
  Kmer();
  Kmer(const vector<string>& k_header, const vector<string>& columns,
       const size_t& id_score);
  Kmer(const vector<string>& k_header, const vector<string>& columns,
       const size_t& id_score, const int& num_dbsnp);

  inline size_t len() const     { return s_seq.size();  }
  // inline int vertex() const  { return i_vertex;      }
  inline string seq() const     { return s_seq;         }
  inline int group() const      { return i_group;       }

  //inline double score(const string& h)
  //  const { return m_score.find(h)->second;  }

  inline int get_dbsnp() const { return i_dbsnp; }
  inline map<string, vector<double> > get_map_score()
    const { return m_score; }

  inline vector<double> get_score(const string& key)
    const { return m_score.find(key)->second; }

  inline map<string, vector<double> >::const_iterator get_score_begin()
    const { return m_score.begin(); }

  inline map<string, vector<double> >::const_iterator get_score_end()
    const { return m_score.end(); }


  vector<string> get_header_score() const;
  string print() const;


  inline void set_vertex(int i) { i_vertex =  i;    }
  inline void set_dbsnp(int i) { i_dbsnp =  i;    }
};

}  // namespace util

#endif  // NEKTAR_CPP_UTIL_KMER_H_
