#include "Kmer.h"
namespace util {

Kmer::Kmer() {
  s_seq = "";
  i_vertex = -1;
  i_group = 0;
}

Kmer::Kmer(const vector<string>& k_header, const vector<string>& columns,
           const size_t& id_score) {
  //cerr << "Create kmer" << endl;
  s_seq = columns[0];
  i_vertex = -1;
  i_group = 0;
  i_dbsnp = -1;

  string column = "";
  for (size_t i=1; i< columns.size(); ++i) {
    if (i >= id_score) {
      m_score.insert(pair<string, vector<double> >(k_header[i],
                     vector<double>(1,stod(columns[i]))));
    } 
    /*
    else {
      m_count.push_back(atoi(columns[i].c_str())));
    }
    */
  }
  if (m_score.find("P_1") != m_score.end()) {
    if (m_score["P_1"][0] > m_score["P_2"][0])
      i_group = 1;
  }
}

Kmer::Kmer(const vector<string>& k_header, const vector<string>& columns,
           const size_t& id_score, const int& num_dbsnp) {
  //cerr << "Create kmer" << endl;
  s_seq = columns[0];
  i_vertex = -1;
  i_group = 0;
  i_dbsnp = num_dbsnp;

  string column = "";
  for (size_t i=1; i< columns.size(); ++i) {
    if (i >= id_score) {
      m_score.insert(pair<string, vector<double> >(k_header[i],
                     vector<double>(1,stod(columns[i]))));
    } 
    /*
    else {
      m_count.push_back(atoi(columns[i].c_str())));
    }
    */
  }
  if (m_score.find("P_1") != m_score.end()) {
    if (m_score["P_1"][0] > m_score["P_2"][0])
      i_group = 1;
  }
}

vector<string> Kmer::get_header_score() const {
  vector<string> res;
  for (map<string, vector<double>>::const_iterator it=m_score.begin();
       it != m_score.end(); ++it) {
    res.push_back(it->first);
  }
  return(res);
}

string Kmer::print() const {
  string res = s_seq;
  for (vector<int>::const_iterator it=v_count.begin();
       it != v_count.end(); ++it) {
    res += "\t" + to_string(*it);
  }
  res += "\n";
  for (map<string, vector<double>>::const_iterator it=m_score.begin();
       it != m_score.end(); ++it) {
    res += "\t" + it->first;
  }
  return(res);
}

}   // namespace util
