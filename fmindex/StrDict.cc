//                              -*- Mode: C++ -*-
// StrDict.cc
//

#include "StrDict.h"

// LIFECYCLE ------------------------------------------------------------

StrDict::StrDict ()
{ }

// OPERATORS ------------------------------------------------------------

// ACCESS ---------------------------------------------------------------

unsigned int StrDict::newId (const char *s) {
  string string_s(s);
  map_t::iterator i = mIndex.find(string_s);
  if (i == mIndex.end()) {
    unsigned int pos = mStr.size();
    mStr.push_back(string_s);
    mIndex.insert(make_pair(string_s, pos));
    return pos;
  } else {
    cerr << "StrDict: an entry for [" << s << "] already exists." << endl;
    abort ();
  }
}

unsigned int StrDict::getId (const char *s) {
  string string_s(s);
  map_t::iterator i = mIndex.find(string_s);
  if (i == mIndex.end ()) {
    unsigned int pos = mStr.size();
    mStr.push_back(string_s);
    mIndex.insert(make_pair(string_s, pos));
    return pos;
  }
  return i->second;
}

bool StrDict::hasId (const char *s) const {
  string string_s(s);
  return mIndex.find(string_s) != mIndex.end();
}

// METHODS --------------------------------------------------------------

void StrDict::clear () {
  mIndex.clear();
  mStr.clear();
}

// I/O  -----------------------------------------------------------------

char *StrDict::readBinary (char *base) {
  unsigned int n = *(unsigned int *)base;
  base += sizeof (unsigned int);
  for (unsigned int i = 0; i < n; ++i) {
    unsigned int l = *(unsigned int *)base;
    base += sizeof (unsigned int);
    string string_base(base);
    mStr.push_back(string_base);
    mIndex.insert(make_pair(string_base, i));
    base += l;
  }

  return base;
}

void StrDict::writeBinary (ostream &out) const {
  unsigned int n = mStr.size ();
  out.write ((char *)&n, sizeof (unsigned int));
  for (unsigned int i = 0; i < n; ++i) {
    unsigned int l = mStr[i].length() + 1;
    out.write ((char *)&l, sizeof (unsigned int));
    out.write (mStr[i].c_str(), l);
  }
}

ostream &StrDict::output (ostream &out) const
{
  for (map_t::const_iterator i = mIndex.begin (); i != mIndex.end (); ++i) {
    cout << '[' << i->first << "] = " << i->second << endl;
  }
  return out;
}

ostream &operator<< (ostream &out, const StrDict &o)
{
  return o.output (out);
}

