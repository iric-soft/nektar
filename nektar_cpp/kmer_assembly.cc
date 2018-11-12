// #include "./conf_nektar.h"
#include "./util/Kmer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <numeric>

#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <sstream>
#include <algorithm> // for min and max
//#include <boost/algorithm/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/connected_components.hpp>

#include "../fmindex/FmIndex.h"

using namespace std;
using namespace boost;
using namespace util;


typedef map<string, Kmer> map_kmer;
typedef map<string, Kmer*> map_pkmer;
typedef vector<Kmer> vec_kmer;
typedef pair<string, string> key_project;
typedef map<key_project, string> map_project;
typedef map<string, unsigned int> map_uint;
typedef map<string, double> map_double;
typedef multimap<double, int> map_assembly;
typedef vector<string> vec_str;
typedef array<string, 4> four_string;
typedef vector<uint64_t> vec_int;
typedef vector<string> vec_str;
typedef vector<double> vec_double;
typedef vector<map<string, vector<double> > > vec_score;
typedef uint64_t two_int[2];
typedef double two_double[2];
typedef unsigned int uint;

typedef string val_param[8];
typedef bool bool_param[8];
const int F_PROJECT = 0;
const int F_KMER = 1;
const int F_OUT = 2;
const int K_LEN = 3;
const int REV = 4;
const int ASSEMBLY = 5;
const int DBSNP = 6;
const int PVALUE = 7;

const char* NUC[5] = {"A", "C", "G", "T", "N"};
const char* GROUP[2] = {"Ref", "Query"};

// typedef for graph
// typedef pair<int, int> edge;
typedef vector<pair<unsigned int, unsigned int> > vec_edge;
typedef vector<string> Vertex;
// create a tag for our new property
enum vertex_seq_t { vertex_seq };
namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, seq);
}
typedef adjacency_list <vecS, vecS, undirectedS> UndirectedGraph;
typedef adjacency_list<vecS, vecS, bidirectionalS,
  no_property, property < edge_weight_t, int > > GraphAssembly;
typedef graph_traits < GraphAssembly >::vertex_descriptor vertex_descriptor;

// ///////////////////////////////////////////////////////
// util function
// ///////////////////////////////////////////////////////
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vec_str split(const string &s, char delim) {
    vec_str elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

char complement(const char c) {
  switch (c) {
    case 'A': return('T');
    case 'T': return('A'); // case 'U':
    case 'C': return('G');
    case 'G': return('C');
    case 'N': return('N');
    /*
    case 'B': return 'V';
    case 'D': return 'H';
    case 'G': return 'C';
    case 'H': return 'D';
    case 'K': return 'M';
    case 'M': return 'K';
    case 'N': return 'N';
    case 'R': return 'Y';
    case 'S': return 'S';
    case 'V': return 'B';
    case 'W': return 'W';
    case 'Y': return 'R';
    */
    default:
      cerr << "ERROR in complement" << endl;
      cerr << c << endl;
      exit(EXIT_FAILURE);
  }
  //throw "parse error";
}

string reverse_complement(const string& seq) {
  //cerr << "reverse_complement: " << seq << endl;
  string rev_seq (seq.size(), 'A');
  string::iterator it_rc = rev_seq.begin();
  for (string::const_reverse_iterator c = seq.rbegin(); c != seq.rend(); ++c){
    *it_rc = complement(*c);
    ++it_rc;
  }
  return(rev_seq);
}

int get_inuc(const string& nuc) {
  if (nuc.size() == 0) return(-1);
  switch (nuc[0]) {
    case 'A': return(0);
    case 'C': return(1);
    case 'G': return(2);
    case 'T': return(3); //case 'U':
    case 'N': return(4);
    default:
      cerr << "ERROR in get_inuc" << endl;
      cerr << nuc << endl;
      exit(EXIT_FAILURE);
  }
  cerr << "ERROR in get_inuc(), dont kmow: " << nuc[0] << endl;
  exit(EXIT_FAILURE);
  return(-1);
}

void calculate_scores(vector<double> vec_all, array<double,4>& res) {
  size_t size_vec = vec_all.size();
  std::sort(vec_all.begin(), vec_all.end());
  //MIN
  res[0] = vec_all[0];
  //MAX
  res[3] = vec_all[size_vec-1];

  //MEDIAN
  if (size_vec  % 2 == 0)
    res[1] = (vec_all[size_vec / 2 - 1] + vec_all[size_vec / 2]) / 2;
  else
    res[1] = vec_all[size_vec / 2];

  //MEAN
  res[2] = std::accumulate(vec_all.begin(), vec_all.end(), 0.0)/size_vec;

}

// ///////////////////////////////////////////////////////
// read function
// ///////////////////////////////////////////////////////

int read_project(const val_param& v_param, map_project& project,
                 vec_str& headers, vec_str& samples, two_int& sum_kmers,
                 vec_int& tt_kmers) {
  int id_switch_group = 0;
  int cpt_switch = -1;
  string::size_type sz = 0;
  filebuf fb;
  if (fb.open(v_param[F_PROJECT], ios::in)) {
    istream is_file(&fb);
    string line = "";
    vec_str columns;
    string column = "";
    string sample = "";
    string tmp_group = "";
    int cpt_line = 0;

    getline(is_file, line);

    headers = split(line, '\t');
    //copy(istream_iterator<string>(iss),
    //     istream_iterator<string>(),
    //     back_inserter(headers));
    //stringstream line_ss(line);
    //split(headers, line, is_any_of("\t, "), token_compress_on);

    while (getline(is_file, line)) {
      columns = split(line, '\t');
      //iss(line);
      //copy(istream_iterator<string>(iss),
      //     istream_iterator<string>(),
      //     back_inserter(columns));
      //split(columns, line, is_any_of("\t, "), token_compress_on);
      for (uint col=0; col<columns.size(); ++col) {
        if (headers[col].compare("SAMPLE") == 0) {
          sample = columns[col];
          samples.push_back(columns[col]);
        } else {
          if (headers[col].compare("GROUP") == 0) {
            if (columns[col].compare(tmp_group) != 0) {
              ++cpt_switch;
              tmp_group = columns[col];
              id_switch_group = cpt_line;
            }
          } else {
            if (headers[col].compare("KMER") == 0) {
              tt_kmers.push_back(stoul(columns[col], &sz, 0));
              sum_kmers[cpt_switch] += tt_kmers[cpt_line];
            }
          }
          project[make_pair(sample, headers[col])] = columns[col];
        }
      }
      ++cpt_line;
    }
  } else {
      cerr << "We can't open the file:" << v_param[F_PROJECT] << endl;
      exit(EXIT_FAILURE);
  }

  if (cpt_switch > 1) {
      cerr << "Need to sort your project file on GROUP column"
           << " or there are more than two groups." << endl;
      exit(EXIT_FAILURE);
  }

  cerr << "samples: " << samples.size()  << endl;
  return(id_switch_group);
}

void read_kmer(const val_param& v_param,
               const vec_str& samples, map_kmer* kmers, 
               map_pkmer* kmers_not_used, vec_str* k_header) {
  string::size_type sz = 0;
  filebuf fb;
  int cpt_kmer = 0;
  if (fb.open(v_param[F_KMER], ios::in)) {
    istream is_file(&fb);
    string line;
    getline(is_file, line);

    int id_score = samples.size() + 1;

    *k_header = split(line, '\t');

    vec_str columns;
    while (getline(is_file, line) && line.size() > 2) {
      columns = split(line, '\t');
      kmers->insert(pair<string, Kmer>(columns[0], Kmer((*k_header), columns, id_score)));
      kmers_not_used->insert(pair<string, Kmer*>(columns[0], NULL));
      ++cpt_kmer;
    }
  } else {
      cerr << "We can't open the file:" << v_param[F_KMER] << endl;
      exit(EXIT_FAILURE);
  }
  cerr << cpt_kmer << " kmers read" << endl;
}

void read_kmer_with_dbsnp(const val_param& v_param,
                          const vec_str& samples, map_kmer* kmers, 
                          map_pkmer* kmers_not_used, vec_str* k_header) {
  string::size_type sz = 0;
  filebuf fb;
  int cpt_kmer = 0;
  int num_dbsnp = -1;

  FmIndex fm_index = FmIndex ("Nuc4");
  fm_index.Load(v_param[DBSNP].c_str());

  if (fb.open(v_param[F_KMER], ios::in)) {
    istream is_file(&fb);
    string line;
    getline(is_file, line);

    int id_score = samples.size() + 1;

    *k_header = split(line, '\t');

    vec_str columns;
    while (getline(is_file, line) && line.size() > 2) {
      columns = split(line, '\t');
      num_dbsnp = atoi(fm_index.CheckSeq(columns[0].c_str()));
      kmers->insert(pair<string, Kmer>(columns[0], Kmer((*k_header), columns, id_score, num_dbsnp)));
      kmers_not_used->insert(pair<string, Kmer*>(columns[0], NULL));
      ++cpt_kmer;
    }
  } else {
      cerr << "We can't open the file:" << v_param[F_KMER] << endl;
      exit(EXIT_FAILURE);
  }

  cerr << cpt_kmer << " kmers read" << endl;
}

// ///////////////////////////////////////////////////////
// graph assembly function
// ///////////////////////////////////////////////////////

bool check_kmer(const map_kmer& kmers, const string& cur_kmer,
               const string& next_kmer, const bool_param& b_param, 
               const val_param& v_param) {
  //cerr << "check_kmer" << endl;
  bool check = false;
  string rev_kmer = "";
  if (cur_kmer.compare(next_kmer) == 0) return(false);
  if (kmers.find(next_kmer) != kmers.end()) {
    if (kmers.at(next_kmer).group() == kmers.at(cur_kmer).group())
      check = true;
  }
  if (!check && b_param[REV]) {
    rev_kmer = reverse_complement(next_kmer);
    if (cur_kmer.compare(rev_kmer) == 0) return(false);
    if (kmers.find(rev_kmer) != kmers.end()) {
      if (kmers.at(rev_kmer).group() == kmers.at(cur_kmer).group())
        check = true;
    }
  }
  if (check && b_param[PVALUE]) {
    vector<double> cur_pvalue = kmers.at(cur_kmer).get_score("PVALUE");
    vector<double> next_pvalue;

    if (rev_kmer.compare("") == 0)  next_pvalue = kmers.at(next_kmer).get_score("PVALUE");
    else                            next_pvalue = kmers.at(rev_kmer).get_score("PVALUE");

    double diff = abs( abs(log10(cur_pvalue[0])) - abs(log10(next_pvalue[0])) );

    if (diff > atoi(v_param[PVALUE].c_str())) {
      // cerr << diff << endl;
      check = false;
    }
  }
  return check;
}

int get_kmers_left(const size_t& kmer_len, const map_kmer& kmers,
                   const string& seq_assembly, const string& seq_cur_kmer,
                   const bool& left, const bool_param& b_param, 
                   const val_param& v_param, four_string* next_kmers) {
  //cerr << "get_kmers_left: " << seq_assembly << endl;
  string next_kmer;
  int seq_len = seq_assembly.size();
  int nb_kmer = 0;
  for (int i=0; i < 4; ++i) {
    if (left) {
      next_kmer = NUC[i] + seq_assembly.substr(0, (kmer_len - 1));
    } else {
      next_kmer = NUC[i] + seq_assembly.substr((seq_len - kmer_len), (kmer_len - 1));
    }

//    if (seq_assembly.compare("TCAAGATCTCTGTCTGGCAGTGGAGGAAGTCT") == 0) {
//      cerr << "L: " << next_kmer << endl;
//    }

    if (check_kmer(kmers, seq_cur_kmer, next_kmer, b_param, v_param)) {
//      if (seq_assembly.compare("TCAAGATCTCTGTCTGGCAGTGGAGGAAGTCT") == 0) {
//        cerr << "OK" << endl;
//      }
      (*next_kmers)[i] = next_kmer;
      ++nb_kmer;
    } else { (*next_kmers)[i] = ""; }
  }
  return(nb_kmer);
}

int get_kmers_right(const size_t& kmer_len, const map_kmer& kmers,
                    const string& seq_assembly, const string& seq_cur_kmer,
                    const bool& left, const bool_param& b_param, 
                    const val_param& v_param, four_string* next_kmers) {
  //cerr << "get_kmers_right" << endl;
  string next_kmer;
  int seq_len = seq_assembly.size();
  int nb_kmer = 0;
  for (int i=0; i < 4; ++i) {
    if (left) {
      next_kmer = seq_assembly.substr(1, (kmer_len - 1)) + NUC[i];
    } else {
      next_kmer = seq_assembly.substr((seq_len - (kmer_len - 1)), (kmer_len - 1)) + NUC[i];
    }

//    if (seq_assembly.compare("TCAAGATCTCTGTCTGGCAGTGGAGGAAGTCT") == 0) {
//      cerr << "R: " << next_kmer << endl;
//    }

    if (check_kmer(kmers, seq_cur_kmer, next_kmer, b_param, v_param)) {
//      if (seq_assembly.compare("TCAAGATCTCTGTCTGGCAGTGGAGGAAGTCT") == 0) {
//        cerr << "OK" << endl;
//      }
      (*next_kmers)[i] = next_kmer;
      ++nb_kmer;
    } else { (*next_kmers)[i] = ""; }
  }
  return(nb_kmer);
}

bool check_edge(uint first_v, uint second_v, vec_str* seq_assembly,
                const size_t& kmer_len) {
  size_t seq_len = (*seq_assembly)[first_v].size();
  string first_s = (*seq_assembly)[first_v].substr(seq_len - (kmer_len - 1), kmer_len - 1);
  string second_s = (*seq_assembly)[second_v].substr(0, kmer_len - 1);

  if (first_s.compare(second_s) == 0) return(true);
  return(false);
}

bool extend(const Kmer& kmer, const string& nuc_backtrack,
            const size_t& kmer_len, map_pkmer* kmers_not_used,
            const map_kmer& kmers, map_uint* k_vertex,
            const bool_param& b_param, const val_param& v_param,
            const bool* left, vec_str* seq_assembly,
            vec_score* score_assembly, vec_int* dbsnp, vec_edge* edges);

void merge(const Kmer& last_kmer, const Kmer& new_kmer,
           vec_str* seq_assembly, vec_score* score_assembly,
           vec_int* dbsnp, vec_edge* edges, map_pkmer* kmers_not_used,
           map_uint* k_vertex, const size_t& i_nuc,
           const map_kmer& kmers, const bool& left,
           const size_t& kmer_len, const bool_param& b_param,
           const val_param& v_param) {
  //cerr << "merge: " << last_kmer.seq() << " " << new_kmer.seq() << endl;
  map_pkmer::iterator it_use = kmers_not_used->find(new_kmer.seq());
  uint l_vertex = (*k_vertex)[last_kmer.seq()];
//  if ((*seq_assembly)[l_vertex].substr(0,21).compare("TCAAGATCTCTGTCTGGCAGT") == 0) {
//    cerr << "merge: " << last_kmer.seq() << " " << new_kmer.seq() << endl;
//    cerr << "merge: " << (*seq_assembly)[l_vertex] << endl;
//  }

  if (it_use != kmers_not_used->end()) {
    kmers_not_used->erase(it_use);

    // new_kmer.set_vertex(vertex);
    (*k_vertex)[new_kmer.seq()] = l_vertex;


    map<string, vector<double> >::const_iterator it;
    for (it=new_kmer.get_score_begin(); it != new_kmer.get_score_end(); ++it) {
      (*score_assembly)[l_vertex].at(it->first).push_back(it->second[0]);
    }

    if (new_kmer.get_dbsnp() > 0) (*dbsnp)[l_vertex] += 1;

    string nuc = NUC[i_nuc];
    string nuc_backtrack = "Z";
    if (left) {
      (*seq_assembly)[l_vertex] = nuc + (*seq_assembly)[l_vertex];
      nuc_backtrack[0] = (*seq_assembly)[l_vertex][kmer_len];
    } else {
      (*seq_assembly)[l_vertex] = (*seq_assembly)[l_vertex] + nuc;
      uint seq_len = (*seq_assembly)[l_vertex].size();
      nuc_backtrack[0] = (*seq_assembly)[l_vertex][seq_len - (kmer_len + 1)];
    }

    extend(new_kmer, nuc_backtrack, kmer_len, kmers_not_used, kmers, k_vertex,
           b_param, v_param, &left, seq_assembly, score_assembly, dbsnp, edges);
  } else {
    if (left) {
      // if (check_edge((*k_vertex)[new_kmer.seq()], l_vertex, seq_assembly, kmer_len))
      edges->push_back(pair<uint, uint>((*k_vertex)[new_kmer.seq()], l_vertex));
    } else {
      // if (check_edge(l_vertex, (*k_vertex)[new_kmer.seq()], seq_assembly, kmer_len))
      edges->push_back(pair<uint, uint>(l_vertex, (*k_vertex)[new_kmer.seq()]));
    }
  }
}

void new_vertex(const Kmer& last_kmer, const Kmer& new_kmer,
                vec_str* seq_assembly, vec_score* score_assembly,
                vec_int* dbsnp, vec_edge* edges, map_pkmer* kmers_not_used,
                const map_kmer& kmers, map_uint* k_vertex,
                const uint& i_nuc, const bool& left,
                const size_t& kmer_len, const bool_param& b_param,
                const val_param& v_param) {
  //cerr << "new_vertex: " << last_kmer.seq() << " " << new_kmer.seq() << endl;
  map_pkmer::iterator it_use = kmers_not_used->find(new_kmer.seq());
  map_uint::iterator it_vertex = k_vertex->find(last_kmer.seq());
  uint l_vertex = it_vertex->second;

  if (it_use != kmers_not_used->end()) {
    kmers_not_used->erase(it_use);

    string tmp_seq;
    uint n_vertex = seq_assembly->size();
    (*k_vertex)[new_kmer.seq()] = n_vertex;

    score_assembly->push_back(new_kmer.get_map_score());

    if (new_kmer.get_dbsnp() > 0) dbsnp->push_back(1);
    else                          dbsnp->push_back(0);

    string nuc = NUC[i_nuc];
    string nuc_backtrack = "Z";
    if (left) {
      tmp_seq = (*seq_assembly)[l_vertex];
      tmp_seq = tmp_seq.substr(0, (kmer_len - 1));
      seq_assembly->push_back(nuc + tmp_seq);
      nuc_backtrack[0] = (*seq_assembly)[l_vertex][kmer_len - 1];

      edges->push_back(pair<uint, uint>(n_vertex, l_vertex));
    } else {
      tmp_seq = (*seq_assembly)[l_vertex];
      int seq_len = tmp_seq.size();
      tmp_seq = tmp_seq.substr((seq_len - (kmer_len - 1)), (kmer_len - 1));
      seq_assembly->push_back(tmp_seq + nuc);
      nuc_backtrack[0] = (*seq_assembly)[l_vertex][seq_len - kmer_len];

      edges->push_back(pair<uint, uint>(l_vertex, n_vertex));
    }

    extend(new_kmer, nuc_backtrack, kmer_len, kmers_not_used, kmers, k_vertex,
           b_param, v_param, &left, seq_assembly, score_assembly, dbsnp, edges);
  } else {
    // cerr << "New Edge!!!" << endl;
    if (left) {
      // if (check_edge((*k_vertex)[new_kmer.seq()], l_vertex, seq_assembly, kmer_len))
      edges->push_back(pair<uint, uint>((*k_vertex)[new_kmer.seq()], l_vertex));
    } else {
      // if (check_edge(l_vertex, (*k_vertex)[new_kmer.seq()], seq_assembly, kmer_len))
      edges->push_back(pair<uint, uint>(l_vertex, (*k_vertex)[new_kmer.seq()]));
    }
  }
}

void backtrack_new_vertex(const Kmer& kmer, vec_str* seq_assembly, 
                          vec_score* score_assembly, vec_int* dbsnp, map_uint* k_vertex,
                          vec_edge* edges, const bool& left,
                          const size_t& kmer_len, const val_param& v_param) {
  //cerr << "backtrack_new_vertex: " << kmer.seq() << endl;
  uint vertex = (*k_vertex)[kmer.seq()];
  int seq_len = (*seq_assembly)[vertex].size();

  uint v_tmp = seq_assembly->size();

  string new_seq;
  if (left) {
      new_seq = (*seq_assembly)[vertex].substr(0, kmer_len);
      (*seq_assembly)[vertex] = (*seq_assembly)[vertex].substr(1, seq_len - 1);
      edges->push_back(pair<uint, uint>(v_tmp, vertex));
  } else {
      new_seq = (*seq_assembly)[vertex].substr((seq_len - kmer_len), kmer_len);
      (*seq_assembly)[vertex] = (*seq_assembly)[vertex].substr(0, (seq_len - 1));
      edges->push_back(pair<uint, uint>(vertex, v_tmp));
  }

  map<string, vector<double> >::const_iterator it;
  for (it=kmer.get_score_begin(); it != kmer.get_score_end(); ++it) {
    (*score_assembly)[vertex].at(it->first).pop_back();// -= it->second;
  }

  

  //kmer.set_vertex(v_tmp);
  (*k_vertex)[kmer.seq()] = v_tmp;

  seq_assembly->push_back(new_seq);
  score_assembly->push_back(kmer.get_map_score());

  if (kmer.get_dbsnp() > 0) {
    dbsnp->push_back(1);
    (*dbsnp)[vertex] -= 1;
  }
  else {
    dbsnp->push_back(0);
  }
}

void extend_many(const Kmer& kmer, vec_str* seq_assembly, 
            vec_score* score_assembly, vec_int* dbsnp, vec_edge* edges,
            map_pkmer* kmers_not_used, const map_kmer& kmers,
            map_uint* k_vertex, const bool& left, const size_t& kmer_len,
            const bool_param& b_param, const val_param& v_param,
            const four_string& kmers_extend) {
  //cerr << "extend_many: " << kmer.seq() << endl;
  map_kmer::const_iterator it_new_kmer;
  for (uint i=0; i < 4; ++i) {
    if (kmers_extend[i].size() != 0) {
      it_new_kmer = kmers.find(kmers_extend[i]);
      if (it_new_kmer == kmers.end() && b_param[REV]) {
        string rev_kmer = reverse_complement(kmers_extend[i]);
        it_new_kmer = kmers.find(rev_kmer);
      }

      new_vertex(kmer, it_new_kmer->second, seq_assembly,
                 score_assembly, dbsnp, edges, kmers_not_used,
                 kmers, k_vertex, i, left, kmer_len, b_param, v_param);
    }
  }
}

void extend_one(const Kmer& kmer, vec_str* seq_assembly, 
            vec_score* score_assembly, vec_int* dbsnp, vec_edge* edges,
            map_pkmer* kmers_not_used, const map_kmer& kmers,
            map_uint* k_vertex, const bool& left, const size_t& kmer_len,
            const bool_param& b_param, const val_param& v_param,
            const four_string& kmers_extend) {
  //cerr << "extend_one: " << kmer.seq() << endl;
  map_kmer::const_iterator it_new_kmer;
  for (uint i=0; i < 4; ++i) {
    if (kmers_extend[i].size() != 0) {
      it_new_kmer = kmers.find(kmers_extend[i]);
      if (it_new_kmer == kmers.end() && b_param[REV]) {
        string rev_kmer = reverse_complement(kmers_extend[i]);
        it_new_kmer = kmers.find(rev_kmer);
      }
      merge(kmer, it_new_kmer->second, seq_assembly,
            score_assembly, dbsnp, edges, kmers_not_used, k_vertex, i,
            kmers, left, kmer_len, b_param, v_param);
      break;
    }
  }
}

bool extend(const Kmer& kmer, const string& nuc_backtrack,
            const size_t& kmer_len, map_pkmer* kmers_not_used,
            const map_kmer& kmers, map_uint* k_vertex,
            const bool_param& b_param, const val_param& v_param,
            const bool* left, vec_str* seq_assembly,
            vec_score* score_assembly, vec_int* dbsnp, vec_edge* edges) {
  //cerr << "extend: " << kmer.seq() << " bt: " << nuc_backtrack << endl;
  four_string kmers_left;
  four_string kmers_right;
  int nb_left = 0;
  int nb_right = 0;
  uint vertex = (*k_vertex)[kmer.seq()];

  if (left == NULL) {
    nb_left = get_kmers_left(kmer_len, kmers, (*seq_assembly)[vertex],
                             kmer.seq(), true, b_param, v_param, &kmers_left);
    nb_right = get_kmers_right(kmer_len, kmers, (*seq_assembly)[vertex],
                               kmer.seq(), false, b_param, v_param, &kmers_right);

//    if ((*seq_assembly)[vertex].substr(0,21).compare("TCAAGATCTCTGTCTGGCAGT") == 0) {
//      cerr << "I: " << (*seq_assembly)[vertex] << endl;
//      for (int i=0; i<4; i++) {
//        cerr << "L" << i << ": " << kmers_left[i] << endl;
//        cerr << "R" << i << ": " << kmers_right[i] << endl;
//      }
//    }

    if (nb_left == 1)
      extend_one(kmer, seq_assembly, score_assembly, dbsnp, edges,
                 kmers_not_used, kmers, k_vertex, true, kmer_len,
                 b_param, v_param, kmers_left);

    if (nb_left > 1)
      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, true, kmer_len,
                  b_param, v_param, kmers_left);

    if (nb_right == 1)
      extend_one(kmer, seq_assembly, score_assembly, dbsnp, edges,
                 kmers_not_used, kmers, k_vertex, false, kmer_len,
                 b_param, v_param, kmers_right);

    if (nb_right > 1)
      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, false, kmer_len,
                  b_param, v_param, kmers_right);

    return(true);
  }

  nb_left = get_kmers_left(kmer_len, kmers, (*seq_assembly)[vertex],
                           kmer.seq(), (*left), b_param, v_param, &kmers_left);
  nb_right = get_kmers_right(kmer_len, kmers, (*seq_assembly)[vertex],
                             kmer.seq(), (*left), b_param, v_param, &kmers_right);

  if (nuc_backtrack.size() != 0) {
    if (*left) {
      if (kmers_right[get_inuc(nuc_backtrack)].size() != 0) {
        kmers_right[get_inuc(nuc_backtrack)] = "";
        --nb_right;
      }
    } else {
      if (kmers_left[get_inuc(nuc_backtrack)].size() != 0) {
        kmers_left[get_inuc(nuc_backtrack)] = "";
        --nb_left;
      }
    }
  }

//  if ((*seq_assembly)[vertex].substr(0,21).compare("TCAAGATCTCTGTCTGGCAGT") == 0) {
//    cerr << (*seq_assembly)[vertex] << " - "<< *left << endl;
//    for (int i=0; i<4; i++) {
//      cerr << "L" << i << ": " << kmers_left[i] << endl;
//      cerr << "R" << i << ": " << kmers_right[i] << endl;
//    }
//  }

  if (*left) {
    if (nb_right >= 1) {
      if ((*seq_assembly)[vertex].size() != kmer_len)
        backtrack_new_vertex(kmer, seq_assembly, score_assembly, dbsnp, 
                             k_vertex, edges, (*left), kmer_len, v_param);

      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, false, kmer_len, b_param,
                  v_param, kmers_right);
    }

    if (nb_left == 1)
      extend_one(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, true, kmer_len, b_param,
                  v_param, kmers_left);

    if (nb_left > 1)
      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, true, kmer_len, b_param,
                  v_param, kmers_left);
  } else {
    if (nb_left >= 1) {
      if ((*seq_assembly)[vertex].size() != kmer_len)
        backtrack_new_vertex(kmer, seq_assembly, score_assembly, dbsnp, 
                             k_vertex, edges, (*left), kmer_len, v_param);

      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, true, kmer_len, b_param,
                  v_param, kmers_left);
    }

    if (nb_right == 1)
      extend_one(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, false, kmer_len, b_param,
                  v_param, kmers_right);

    if (nb_right > 1)
      extend_many(kmer, seq_assembly, score_assembly, dbsnp, edges,
                  kmers_not_used, kmers, k_vertex, false, kmer_len, b_param,
                  v_param, kmers_right);
  }
  return(true);
}

void linear_assembly(const map_kmer& kmers,
                     const vec_str& seq_assembly,
                     const vec_score& score_assembly,
                     const vec_int& dbsnp,
                     const GraphAssembly& G,
                     const vector<int>& component,
                     const size_t& kmer_len,
                     const vec_str& v_samples,
                     const bool_param& b_param,
                     const val_param& v_param) {
  //cerr << "linear_assembly" << endl;
  std::ofstream os_fasta;
  std::ofstream os_tab;
  os_fasta.open(v_param[F_OUT] + "/assembly.fasta");
  os_tab.open(v_param[F_OUT] + "/assembly.tab");
  vec_str score_headers = kmers.begin()->second.get_header_score();
  if (!os_tab.is_open()) {
    std::cerr << "Can't open/write in:"
              << v_param[F_OUT] << std::endl;
    exit(EXIT_FAILURE);
  } else {
    os_tab << "ID_SEQ\tID_GRAPH\tGROUP";
    for (vec_str::const_iterator it = score_headers.begin(); it != score_headers.end(); ++it) {
      os_tab << "\tMIN_"; os_tab << it->c_str();
      os_tab << "\tMEDIAN_" << it->c_str();
      os_tab << "\tMEAN_" << it->c_str();
      os_tab << "\tMAX_" << it->c_str();
    }

    if (b_param[DBSNP]) os_tab << "\tDB_SNP";

    os_tab << "\tLEN\tIN_DEGREE\tOUT_DEGREE\n";
  }

  string seq;
  string seq_kmer;
  string group;
  int cpt_id = 1;
  map<string, vector<double> >::const_iterator it;
  array<double,4> scores;

  for (size_t i=0; i != seq_assembly.size(); ++i) {
    seq = seq_assembly[i];
    seq_kmer = seq.substr(0, kmer_len);
    if (kmers.find(seq_kmer) != kmers.end())
      group = GROUP[kmers.at(seq_kmer).group()];
    else
      group = GROUP[kmers.at(reverse_complement(seq_kmer)).group()];

    os_fasta << ">" << cpt_id << "\n";
    os_fasta << seq << "\n";
    os_tab << cpt_id
           << "\t" << component[i]
           << "\t" << group;
    for (it=score_assembly[i].begin();
         it != score_assembly[i].end(); ++it) {
      calculate_scores(it->second, scores);
      os_tab << "\t" << scores[0];
      os_tab << "\t" << scores[1];
      os_tab << "\t" << scores[2];
      os_tab << "\t" << scores[3];
    }

    if (b_param[DBSNP]) os_tab << "\t" << dbsnp[i];

    os_tab << "\t" << seq.size()
           << "\t" << in_degree(i, G)
           << "\t" << out_degree(i, G)
           << "\n";

    ++cpt_id;
  }
}


void create_graph_assembly(const map_kmer& kmers,
                           map_pkmer* kmers_not_used,
                           const vec_str& v_samples,
                           const bool_param& b_param,
                           const val_param& v_param) {
  map_pkmer::iterator it_use;
  map_kmer::const_iterator it_kmer;
  vec_str seq_assembly;
  // vec_double score_assembly;
  vec_score score_assembly;
  vec_edge edges;
  vec_int dbsnp;
  map_uint k_vertex;
  int cpt_graph = 0;
  size_t kmer_len = kmers.begin()->second.len();
  string kmer_seq;
  while (!kmers_not_used->empty()) {
    //cerr << "NEW graph!!!!" << endl;
    //cerr << kmers.size() << " - " << kmers_not_used->size() << endl;
    //cerr << seq_assembly.size() << endl;

    it_use = kmers_not_used->begin();
    it_kmer = kmers.find(it_use->first);
    kmers_not_used->erase(it_use);

    kmer_seq = it_kmer->first;
    k_vertex[it_kmer->second.seq()] = score_assembly.size();
    // it_kmer.second.set_vertex(score_assembly.size());
    seq_assembly.push_back(kmer_seq);
    score_assembly.push_back(it_kmer->second.get_map_score());
    if (it_kmer->second.get_dbsnp() > 0) dbsnp.push_back(1);
    else                                 dbsnp.push_back(0);

    extend(it_kmer->second, "", kmer_len, kmers_not_used, kmers,
           &k_vertex, b_param, v_param, NULL, &seq_assembly,
           &score_assembly, &dbsnp, &edges);

    ++cpt_graph;
  }

  cerr << "End create graph: "<< cpt_graph << endl;
  cerr << "vertex: "<< score_assembly.size() << endl;

  GraphAssembly G(score_assembly.size());
  UndirectedGraph UG(score_assembly.size());
  for (size_t i=0; i < edges.size(); ++i) {
    add_edge(edges[i].first, edges[i].second, G);
    add_edge(edges[i].first, edges[i].second, UG);
  }

  vector<int> component(num_vertices(UG));
  int num = connected_components(UG, &component[0]);
  cerr << "number of connected components: " << num << endl;

  linear_assembly(kmers, seq_assembly, score_assembly, dbsnp, G, 
                  component, kmer_len, v_samples, b_param, v_param);
}

// ///////////////////////////////////////////////////////
// main function
// ///////////////////////////////////////////////////////
void print_help(string str_error = "") {
  if (str_error.compare("") != 0) {
    cerr << "Error: "<< str_error << endl;
  }
  cerr << "Usage: -p file -k file -o file -d [file] [-l int] [-r]" << endl;
  cerr << "-p file: input project file" << endl;
  cerr << "-k file: input kmer file" << endl;
  cerr << "-o file: output assembly file" << endl;
  cerr << "[-v int]: diff max between log10(pvalue) of 2 overlapped kmer" << endl;
  cerr << "[-d file]: input dbSNP file" << endl;
  cerr << "[-l int]: length of k-mer (default: 31)" << endl;
  cerr << "[-r ]: check the reverse complement during k-mer "
       << "assembly step" << endl;
  exit(EXIT_FAILURE);
}

void get_args(int argc, char *argv[], bool_param* b_param,
              val_param* v_param) {
  int opt;
  (*v_param)[K_LEN] = "31";
  (*v_param)[ASSEMBLY] = "1";
  (*v_param)[PVALUE] = "0";
  while ((opt = getopt(argc, argv, "p:k:o:v:d:l:r")) != -1) {
    switch (opt) {
    case 'p':
      (*v_param)[F_PROJECT] = optarg;
      break;
    case 'k':
      (*v_param)[F_KMER] = optarg;
      break;
    case 'o':
      (*v_param)[F_OUT] = optarg;
      break;
    case 'v':
      (*v_param)[PVALUE] = optarg;
      (*b_param)[PVALUE] = true;
      if (atoi((*v_param)[PVALUE].c_str()) < 0) {
        print_help("Value after -v need to be an int >= 0");
      }
      break;
    case 'd':
      (*v_param)[DBSNP] = optarg;
      (*b_param)[DBSNP] = true;
      break;
    case 'l':
      (*v_param)[K_LEN] = optarg;
      if (atoi((*v_param)[K_LEN].c_str()) <= 0) {
        print_help("Value after -l need to be an int > 0");
      }
      break;
    case 'r':
      (*b_param)[REV] = true;
      break;
    default:
      print_help();
    }
  }
}


int main(int argc, char *argv[]) {
  // chose_bidon (0);
  clock_t start, end;
  start = clock();
  // Check number of input files
  string file_out;
  string file_in;
  bool_param b_param = {false};
  val_param v_param;

  cerr << "Read args" << endl;
  get_args(argc, argv, &b_param, &v_param);

  vec_str p_headers;
  vec_str k_header;
  vec_str v_samples;
  map_project m_project;
  map_kmer kmers;
  map_pkmer kmers_not_used;

  two_int sum_kmers = {};
  vec_int tt_kmers;
  cerr << "Read project" << endl;
  //int id_switch_group = read_project(v_param, m_project, p_headers, v_samples,
  //                                   sum_kmers, tt_kmers);
  read_project(v_param, m_project, p_headers, v_samples,
                                     sum_kmers, tt_kmers);

  cerr << "Read kmer" << endl;
  if (b_param[DBSNP]) {
    read_kmer_with_dbsnp(v_param, v_samples, &kmers, &kmers_not_used, &k_header);
  } else {
    read_kmer(v_param, v_samples, &kmers, &kmers_not_used, &k_header);
  }

  end = clock();
  cerr  <<  "Time : " << static_cast<float>(end-start)/CLOCKS_PER_SEC << "\n";

  cerr << "Create graph assembly" << endl;
  create_graph_assembly(kmers, &kmers_not_used, v_samples, b_param, v_param);


  end = clock();
  cerr  <<  "Time : " << static_cast<float>(end-start)/CLOCKS_PER_SEC << "\n";
  return 0;
}
