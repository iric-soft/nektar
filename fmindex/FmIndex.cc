//                              -*- Mode: C++ -*-
// FmIndex.cc
//

#include "FmIndex.h"

// LIFECYCLE ------------------------------------------------------------

FmIndex::FmIndex(): mSeq(NULL), mL(NULL), mOcc(NULL), mC(NULL), mPos(NULL), mNew(false)
{
  mBaseRep = new Nuc4_BaseRep();

  mSeq = SeqStore(mBaseRep);
}

FmIndex::FmIndex(const char *base): mSeq(NULL), mL(NULL), mOcc(NULL), mC(NULL), mPos(NULL), mNew(false)
{
  if (strcmp(base, "Nuc4") == 0)
    mBaseRep = new Nuc4_BaseRep();
  if (strcmp(base, "AA8") == 0)
    mBaseRep = new AA8_BaseRep();

  mSeq = SeqStore(mBaseRep);
}

FmIndex::FmIndex(const char *genome_fn, const char *base, const bool rand, const bool reverse): mSeq(NULL), mL(NULL), mOcc(NULL), mC(NULL), mPos(NULL), mNew(false)
{
  if (strcmp(base, "Nuc4") == 0)
    mBaseRep = new Nuc4_BaseRep();
  if (strcmp(base, "AA8") == 0)
    mBaseRep = new AA8_BaseRep();

  mSeq = SeqStore(mBaseRep);

  mSeq.readBinaryOrFasta(genome_fn);
  mSeq.detach();
  if (strcmp(base, "AA8") != 0)
    if (rand || strcmp(base, "Nuc2") == 0)
      mSeq.randomize_n();

  if (reverse)
    mSeq.unfold();

  mSeqLength = mSeq.size();
  cout << "Store size after unfold: " << mSeq.size() << endl;
}

FmIndex::~FmIndex() {
  delete mBaseRep;
  DeleteArray();
}

// OPERATORS ------------------------------------------------------------

// ACCESS ---------------------------------------------------------------

// METHODS --------------------------------------------------------------
void FmIndex::CreateIndex(unsigned int nBucketSize, unsigned int nMarkFreq) {
  DeleteArray();

  //check and create tmp directory
  struct stat st;
  if(stat("./tmp_fmindex",&st) != 0)
    system("mkdir ./tmp_fmindex");

  CreateSortedArray();

  // Allocate the FM index
  cout << "Allocation" << endl;

  mBucketSize = nBucketSize;
  mNbBucket = mSeqLength / mBucketSize + 1;

  mMarkFreq = nMarkFreq;
  mNbMark = mSeqLength / mMarkFreq + 1;

  mNew = true;
  mL = new seq_t[mSeq.wsize()];
  //cout << "mL:" << (void*)mL << endl;
  //cout << "mL_size: " << mSeq.wsize () << endl;

  mOcc = new count_t[mNbBucket * mBaseRep->alphabetSize];
  mC = new position_t[mBaseRep->alphabetSize + 1];
  mPos = new position_t[mNbMark];
  for (size_t i = 0; i < mBaseRep->alphabetSize + 1; ++i)
    mC[i] = mSeqLength;


  // Construct L
  cout << "Building the index" << endl;

  seq_t tmp_l = 0;
  count_t count_l = 0;
  position_t pos_l = 0;

  count_t symbol_count[mBaseRep->alphabetSize];  // Used for Occ
  for (unsigned int i = 0; i < mBaseRep->alphabetSize; ++i) symbol_count[i] = 0;

  seq_t last_base = 1 << mBaseRep->bitPerBase;  // a symbol that can't occur

  ifstream final_in("./tmp_fmindex/o_all.n");
  for (position_t i = 0; i < mSeqLength; ++i) {

    if (i % 10000000UL == 0) {
      cout << "i = " << i << endl;
    }

    position_t o_i = 42;
    final_in.read((char*)&o_i, sizeof(position_t));

    // L
    position_t prev_pos = o_i + mSeqLength - 1;
    if (prev_pos >= mSeqLength) prev_pos -= mSeqLength;
    count_t base = mSeq.getBit(SeqIter (prev_pos));
    count_t base_index = mBaseRep->bitToIndex(base);

    tmp_l <<= mBaseRep->bitPerBase;
    tmp_l += base;
    count_l++;
    if (count_l == mBaseRep->basePerWord) {
      mL[pos_l] = tmp_l;
      pos_l++;
      tmp_l = 0;
      count_l = 0;
    }

    // Occ
    position_t bucket_id = i / mBucketSize;
    if (bucket_id * mBucketSize == i) {
      count_t *ptr = &(mOcc[bucket_id * mBaseRep->alphabetSize]);
      for (unsigned int j = 0; j < mBaseRep->alphabetSize; ++j) {
        ptr[j] = symbol_count[j];
      }
    }

    // Pos
    position_t mark_id = i / mMarkFreq;
    if (mark_id * mMarkFreq == i) {
      mPos[mark_id] = o_i;  // Mark this position for retrieval
    }

    symbol_count[base_index]++;

    // C
    base = mSeq.getBit(SeqIter(o_i));
    base_index = mBaseRep->bitToIndex(base);
    if (base != last_base) {
      mC[base_index] = i;
      last_base = base;
    }
  }
  mC[mBaseRep->alphabetSize] = mSeqLength;

  final_in.close();
  system("rm -fr ./tmp_fmindex");
}

void FmIndex::CreateSortedArray() {
  cout << "Creating the sorted array..." << endl;

  size_t prefix_size = 2; // Increase to reduce memory footprint...
  comp_seq comp(&mSeq);

  size_t prefix_n = size_t(pow(mBaseRep->alphabetSize, prefix_size));
  cout << "prefix_n = " << prefix_n << endl;

  vector< ofstream > prefix_f(prefix_n);
  omp_lock_t prefix_lock[prefix_n];
  size_t prefix_count[prefix_n];

  for (size_t i = 0; i < prefix_n; ++i) {
    char buffer[4096];
    sprintf(buffer, "./tmp_fmindex/o_%lu.n", i);
    prefix_f[i].open(buffer);
    omp_init_lock(&(prefix_lock[i]));
    prefix_count[i] = 0;
  }

  #pragma omp parallel for
  for (position_t i = 0; i < mSeq.size(); ++i) {
    size_t prefix = 0;
    for (size_t j = 0; j < prefix_size; ++j) {
      size_t circ_ij = i + j;
      if (circ_ij >= mSeq.size()) circ_ij -= mSeq.size();
      prefix = (prefix * mBaseRep->alphabetSize) + mBaseRep->bitToIndex(mSeq.getBit(SeqIter(circ_ij)));
    }
    omp_set_lock(&(prefix_lock[prefix]));
    prefix_f[prefix].write((char*)&i, sizeof(position_t));
    // if (prefix == 1) {
    //   cout << "Writing [" << i << "]" << endl;
    // }
    prefix_count[prefix]++;
    omp_unset_lock(&(prefix_lock[prefix]));
  }

  for (size_t i = 0; i < prefix_n; ++i) {
    prefix_f[i].close();
  }

  cout << "Sorting..." << endl;

#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < prefix_n; ++i) {
    if (prefix_count[i] > 0) {
      position_t *tmp = new position_t[prefix_count[i]];
      char buffer[4096];
      sprintf(buffer, "./tmp_fmindex/o_%lu.n", i);
      ifstream tmp_in(buffer);
      tmp_in.read((char*)tmp, sizeof(position_t) * prefix_count[i]);
      tmp_in.close();
      //cout << "sorting " << i << " (" << prefix_count[i] << ")" << endl;
      // for (size_t j = 0; j < prefix_count[i]; ++j)
      //   cout << "  " << *(tmp + j) << endl;
      sort(tmp, tmp + prefix_count[i], comp);
      //cout << "done " << i << endl;
      ofstream tmp_out(buffer);
      tmp_out.write((char*)tmp, sizeof(position_t) * prefix_count[i]);
      tmp_out.close();
      delete[] tmp;
    }
  }

  cout << "Merging..." << endl;

  ofstream all_f("./tmp_fmindex/o_all.n");
  for (size_t i = 0; i < prefix_n; ++i) {
    if (prefix_count[i] > 0) {
      position_t *tmp = new position_t[prefix_count[i]];
      prefix_f[i].close();
      char buffer[4096];
      sprintf(buffer, "./tmp_fmindex/o_%lu.n", i);
      ifstream tmp_in(buffer);
      tmp_in.read((char*)tmp, sizeof(position_t) * prefix_count[i]);
      all_f.write((char*)tmp, sizeof(position_t) * prefix_count[i]);
      delete[] tmp;
    }
  }
  all_f.close();

  //cout << "done with the merge." << endl;
}


void FmIndex::RebuildOcc(count_t bucketSize) {
  if (mNew) delete[] mOcc;

  mBucketSize = bucketSize;
  mNbBucket = mSeqLength / mBucketSize + 1;

  mOcc = new count_t[mNbBucket * mBaseRep->alphabetSize];

  count_t symbol_count[mBaseRep->alphabetSize];
  for (size_t i = 0; i < mBaseRep->alphabetSize; ++i) symbol_count[i] = 0;

  for (position_t i = 0; i < mSeqLength; ++i) {

    if (i % 10000000UL == 0) {
      cout << (i / 1e6) << " / " << (mSeqLength / 1e6) << endl;
    }

    // Occ
    position_t bucket_id = i / mBucketSize;
    if (bucket_id * mBucketSize == i) {
      count_t *ptr = &(mOcc[bucket_id * mBaseRep->alphabetSize]);
      for (unsigned int j = 0; j < mBaseRep->alphabetSize; ++j) {
        ptr[j] = symbol_count[j];
      }
    }

    count_t base = mBaseRep->getBit(mL, i, false);
    count_t base_index = mBaseRep->bitToIndex(base);
    symbol_count[base_index]++;
  }

}


position_t FmIndex::Occ(count_t c_bit, count_t c_index, position_t pos) const {
  position_t bucket_id = pos / mBucketSize;
  position_t count = (position_t)mOcc[bucket_id * mBaseRep->alphabetSize + c_index];
  // count = 0;

  seq_t *base = &(mL[bucket_id * (mBucketSize / mBaseRep->basePerWord)]);
  return count + mBaseRep->countBit(base, pos - (bucket_id * mBucketSize), c_bit);

  // for (position_t i = (bucket_id * mBucketSize); i <= pos; ++i) {
  //   if (mBaseRep->getBit (mL, i, false) == c_bit) {
  //     count++;
  //   }
  // }

  // cout << count << "]" << endl;

  // return count;

}

position_t FmIndex::Locate(position_t i) const {
  position_t mark_id = i / mMarkFreq;
  position_t v = 0;
  count_t c_bit = mBaseRep->getBit(mL, i, false);
  count_t c_index = mBaseRep->bitToIndex(c_bit);

  while ((mark_id * mMarkFreq) != i) {
    i = mC[c_index] + Occ(c_bit, c_index, i - 1);
    mark_id = i / mMarkFreq;
    c_bit = mBaseRep->getBit(mL, i, false);
    c_index = mBaseRep->bitToIndex(c_bit);
    v++;
  }

  position_t pos = mPos[mark_id] + v;
  if (pos >= mSeqLength) {
    //cerr << "*********************** " << mSeqLength << endl;
    pos -= mSeqLength; // Should not happen since there is a mark at 0
  }
  return pos;
}

void FmIndex::InitLookup(count_t c_bit, position_t *out_sp, position_t *out_ep) const {
  count_t c_index = mBaseRep->bitToIndex(c_bit);
  *out_sp = mC[c_index];
  *out_ep = mC[c_index + 1] - 1;
}

void FmIndex::Lookup(position_t cur_sp, position_t cur_ep, count_t c_bit,
                      position_t *out_sp, position_t *out_ep) const
{
  if (cur_ep < cur_sp) {
    *out_sp = cur_sp;
    *out_ep = cur_ep;
    return;
  }
  count_t c_index = mBaseRep->bitToIndex(c_bit);

  position_t tmp_sp = mC[c_index] + ((cur_sp == 0) ? 0 : Occ(c_bit, c_index, cur_sp - 1));
  position_t tmp_ep = mC[c_index] + Occ(c_bit, c_index, cur_ep) - 1;

  *out_sp = tmp_sp;
  *out_ep = tmp_ep;
}

const char* FmIndex::CheckSeq(const char *seq) const {
  position_t begin = 0;
  position_t end = 0;
  position_t cpt = strlen(seq) - 1;
  unsigned int bit_char = mBaseRep->baseToBit(seq[cpt]);
  InitLookup(bit_char, &begin, &end);
  while (cpt > 0) {
    cpt--;
    if (end >= begin) {
      bit_char = mBaseRep->baseToBit(seq[cpt]);
      Lookup(begin, end, bit_char, &begin, &end);
    }
  }

  string str_out = "";
  if (end < begin) {
    str_out = "0";
  } else {
    ostringstream convert;
    convert << int(end - begin + 1);
    str_out = convert.str();
  }
  // Trick need a copy : str_out.c_str() only create a bug in python
  char *cstr = new char[str_out.length() + 1];
  strcpy(cstr, str_out.c_str());
  return(cstr);
}

const char* FmIndex::SearchSeq(const char *seq, int max_match) const {
  position_t begin = 0;
  position_t end = 0;
  position_t cpt = strlen(seq) - 1;
  position_t loc = 0;
  unsigned int bit_char = mBaseRep->baseToBit(seq[cpt]);
  InitLookup(bit_char, &begin, &end);
  while (cpt > 0) {
    cpt--;
    // cout << "b: " << begin << " e: " << end << endl;
    if (end >= begin) {
      // cout << seq[cpt] << endl;
      bit_char = mBaseRep->baseToBit(seq[cpt]);
      Lookup(begin, end, bit_char, &begin, &end);
    }
  }
  //cout << "b: " << begin << " e: " << end << endl;
  //cout << "end search seq" << endl;

  string str_out = "";
  if (end < begin) {
    str_out = "0";
  } else {
    int num_match = static_cast<int>(end - begin + 1);
    ostringstream convert;
    convert << num_match;
    str_out = convert.str();

    if (num_match <= max_match) {
      for (cpt=begin; cpt <= end; cpt++) {
        loc = Locate(cpt);
        str_out += ",";
        str_out += mSeq.output(SeqRange(loc, loc+strlen(seq)));
      }
    }
  }
  // Trick need a copy : str_out.c_str() only create a bug in python
  char *cstr = new char[str_out.length() + 1];
  strcpy(cstr, str_out.c_str());
  return(cstr);
}

// // I/O  -----------------------------------------------------------------

void FmIndex::Save(const char *fn) {
  cout << "Saving the index..." << endl;
  //cout << "A: " << mC[1] << endl;
  //cout << "C: " << mC[2] << endl;
  //cout << "G: " << mC[4] << endl;
  //cout << "T: " << mC[8] << endl;

  ofstream out(fn);
  mSeq.writeBinary(out);
  out.write((char*)&mSeqLength, sizeof(position_t));

  out.write((char*)&mBucketSize, sizeof(count_t));
  out.write((char*)&mNbBucket, sizeof(count_t));

  out.write((char*)&mMarkFreq, sizeof(count_t));
  out.write((char*)&mNbMark, sizeof(count_t));

  long_write(out, (char *)mL, mSeq.wsize() * sizeof(seq_t));
  out.write((char*)mOcc, (position_t)mNbBucket * mBaseRep->alphabetSize * sizeof(count_t));
  out.write((char*)mC, (mBaseRep->alphabetSize + 1) * sizeof(position_t));
  long_write(out, (char*)mPos, (size_t)mNbMark * sizeof(position_t));

  //cout << "  Store size = " << mSeq.size() << endl;
  //cout << "  BucketSize = " << mBucketSize << endl;
  //cout << "  MarkFreq = " << mMarkFreq << endl;

}


void FmIndex::Load(const char *fn) {
  DeleteArray();

  int fd = open(fn, O_RDONLY);
  struct stat sb;
  fstat(fd, &sb);  // File size

  // POPULATE removed since it is unsupported on OSX. Need to assess performance loss.
  // char *addr = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED | MAP_POPULATE, fd, 0);
  char *addr = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  char *ptr = addr;

  ptr = mSeq.readBinaryMM(ptr);

  mSeqLength = *(position_t*)ptr;  ptr += sizeof(position_t);
  mBucketSize = *(count_t*)ptr;    ptr += sizeof(count_t);
  mNbBucket = *(count_t*)ptr;      ptr += sizeof(count_t);
  mMarkFreq = *(count_t*)ptr;      ptr += sizeof(count_t);
  mNbMark = *(count_t*)ptr;        ptr += sizeof(count_t);

  //cout << "  Store size = " << mSeq.size() << endl;
  //cout << "  BucketSize = " << mBucketSize << endl;
  //cout << "  MarkFreq = " << mMarkFreq << endl;

  mL =   (seq_t*)ptr;       ptr += (position_t)(mSeq.wsize() * sizeof(seq_t));
  mOcc = (count_t*)ptr;     ptr += (position_t)mNbBucket * mBaseRep->alphabetSize * sizeof (count_t);
  mC =   (position_t*)ptr;  ptr += (mBaseRep->alphabetSize + 1) * sizeof(position_t);
  mPos = (position_t*)ptr;  ptr += (position_t)mNbMark * sizeof(position_t);

  //cout << "A: " << mC[1] << endl;
  //cout << "C: " << mC[2] << endl;
  //cout << "G: " << mC[4] << endl;
  //cout << "T: " << mC[8] << endl;

}


string FmIndex::GetLocateChr(position_t i) const{
  return mSeq.output(SeqIter(Locate(i)));
}


bool FmIndex::GetLocateStrand(position_t i) const{
  return mSeq.forward(SeqIter(Locate(i)));
}


// Private --------------------------------------------------------------

FmIndex::comp_seq::comp_seq(SeqStore *store_ptr) {
  begin = store_ptr->begin().getOffset();
  end = store_ptr->end().getOffset();
  base = store_ptr->getStore();
  rep = store_ptr->getBaseRep();
}


void FmIndex::DeleteArray(){
  if(mNew){
    if (mL) delete[] mL;
    if (mOcc) delete[] mOcc;
    if (mC) delete[] mC;
    if (mPos) delete[] mPos;
    mNew = false;
  }
}

bool FmIndex::comp_seq::operator()(position_t a, position_t b) {
  if (a == b) return false;

  count_t bit_a = rep->getBit(base, a, false);
  count_t bit_b = rep->getBit(base, b, false);
  while (bit_a == bit_b) {
    ++a;
    ++b;
    if (a == end) a = begin;
    if (b == end) b = begin;
    bit_a = rep->getBit (base, a, false);
    bit_b = rep->getBit (base, b, false);
  }

  return bit_a < bit_b;
  // if (a == b) return false;
  // SeqIter ai(a);
  // SeqIter bi(b);
  // SeqIter end = store->end();

  // while (store->getBit(ai) == store->getBit(bi)) {
  //   ++ai;
  //   ++bi;
  //   if (ai == end) ai = store->begin();
  //   if (bi == end) bi = store->begin();
  // }

  // return store->getBit(ai) < store->getBit(bi);
}
