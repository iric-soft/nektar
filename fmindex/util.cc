//                              -*- Mode: C++ -*-
// util.cc
//

#include "util.h"

void long_write (ostream &out, const char *c, size_t n) {
  size_t max_chunk = (1L << 30) - 1;
  size_t written = 0;
  while (written < n) {
    size_t chunk = min (max_chunk, n - written);
    out.write (c + written, chunk);
    if (out.fail ()) {
      cerr << "long_write: Write failed" << endl;
      return;
    }
    //cout << "Wrote " << chunk << endl;
    written += chunk;
  }
}
