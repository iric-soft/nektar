//                              -*- Mode: C++ -*-
// StrDict.h
//

#ifndef __StrDict_h_
#define __StrDict_h_

#include <stdlib.h>
#include <iostream>
#include <map>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <string>

using namespace std;

class StrDict
{

  struct ltstr {
    bool operator() (const string& s1, const string& s2) const {
      return strcmp(s1.c_str(), s2.c_str()) < 0;
    }
  };

  typedef map< string, unsigned int, ltstr > map_t;
  map_t mIndex;
  vector<string> mStr;

public:

  typedef map_t::iterator iterator;
  typedef map_t::const_iterator const_iterator;

  // LIFECYCLE ------------------------------------------------------------

  StrDict();

  // OPERATORS ------------------------------------------------------------

  // ACCESS ---------------------------------------------------------------

  unsigned int getId(const char *s);
  unsigned int newId(const char *s);
  bool hasId (const char *s) const ;
  const char* getStr(unsigned int i) const { return mStr[i].c_str(); }

  iterator begin() { return mIndex.begin(); }
  iterator end() { return mIndex.end(); }
  const_iterator begin() const { return mIndex.begin(); }
  const_iterator end() const { return mIndex.end(); }
  unsigned int size() const { return mStr.size(); }

  // METHODS --------------------------------------------------------------

  void clear();

  // I/O  -----------------------------------------------------------------

  char* readBinary(char *base);
  void writeBinary(ostream &out) const ;

  ostream &output(ostream &out) const;

};

ostream &operator<< (ostream &out, const StrDict &o);

#endif
