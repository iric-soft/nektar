//                              -*- Mode: C++ -*-
// BaseRep.cc
//

#include <stdio.h>

#include "BaseRep.h"

using namespace std;
// GENERIC ----------------------------------------------------------------

count_t BaseRep::getBit (seq_t *base, position_t offset, bool comp) const {
  // cout << "offset = " << offset << endl;
  // cout << "  i = " << (offset >> offsetShift) << endl;
  // cout << "  shift_nuc = " << (~offset & offsetMask) << endl;
  // cout << "  shift_bit = " << ((~offset & offsetMask) << logBitPerBase) << endl;

  // printf("==\nAvt: base, %d, offset %d, shift: %d\n", base, offset, offsetShift);
  base += (offset >> offsetShift);
  //cout <<"base app: " << base << endl;
  // printf("App: base, %d, offset %d, shift: %d\n", base, offset, offsetShift);
  count_t b = (*base >> ((~offset & offsetMask) << logBitPerBase)) & bitMask;
  // cout << "  b = " << b << endl;
  // cout << endl;
  if (comp) b = complement (b);
  // count_t i = offset / basePerWord;
  // base += i;
  // seq_t tmp = (*base) << ((offset % basePerWord) * bitPerBase);
  // count_t b = tmp >> (sizeof (seq_t) * 8 - bitPerBase);
  // if (comp) b = complement (b);
  return b;
}

char BaseRep::getBase (seq_t *base, position_t offset, bool comp) const {
  return bitToBase (getBit (base, offset, comp));
}

void BaseRep::setBit (seq_t *base, position_t offset, bool comp, seq_t b) {
  if (comp) b = complement (b);

  count_t i = offset / basePerWord;
  base += i;
  count_t ii = (offset % basePerWord) * bitPerBase;

  count_t shift = (sizeof (seq_t) * 8 - ii - bitPerBase);
  seq_t mask = bitMask << shift;

  *base = ((*base) & (~mask)) | (b << shift);
}

void BaseRep::setBase (seq_t *base, position_t offset, bool comp, char c) {
  seq_t b = baseToBit (c);
  setBit (base, offset, comp, b);
}

seq_t BaseRep::word_complement (seq_t w) const {
  seq_t res = 0;
  for (size_t i = 0; i < basePerWord; ++i) {
    res = (res << bitPerBase) + complement ((w >> (i * bitPerBase)) & bitMask);
  }
  return res;
}

count_t BaseRep::countBit (seq_t *base, position_t offset, seq_t b) const {
  // Counts the number of b between the beginning of base to the offset
  position_t cur_offset = 0;
  seq_t mask_zero = (((seq_t)1 << bitPerBase) - 1) << ((basePerWord - 1) * bitPerBase);  // 1111000000000000...
  count_t count = 0;
  while (cur_offset <= offset) {
    seq_t mask_cur = mask_zero;
    seq_t b_cur = b << ((basePerWord - 1) * bitPerBase);
    for (size_t i = 0; i < basePerWord; ++i) {
      if ((*base & mask_cur) == b_cur) count++;
      mask_cur >>= bitPerBase;
      b_cur >>= bitPerBase;
      cur_offset++;
      if (cur_offset > offset) break;
    }
    base++;
  }
  return count;
}


// Nuc4 - SPECIFICS -------------------------------------------------------

unsigned int Nuc4_BaseRep::baseToBit (char c) const {
  switch (c) {
  case 'A':
  case 'a': return 1;
  case 'C':
  case 'c': return 2;
  case 'G':
  case 'g': return 4;
  case 'T':
  case 't': return 8;
  case 'N':
  case 'n': return 15;
  }
  return 0;
}

char Nuc4_BaseRep::bitToBase (seq_t b) const {
  switch (b) {
  case 0: return '.';
  case 1: return 'A';
  case 2: return 'C';
  case 4: return 'G';
  case 8: return 'T';
  case 15: return 'N';
  }
  return '?';
}

count_t Nuc4_BaseRep::complement (count_t b) const {
  switch (b) {
  case 0: return 0;
  case 1: return 8;
  case 2: return 4;
  case 4: return 2;
  case 8: return 1;
  case 15: return 15;
  }
  return 0;
}

count_t Nuc4_BaseRep::bitToIndex (count_t b) const {
  switch (b) {
  case 0: return 0;
  case 1: return 1;
  case 2: return 2;
  case 4: return 3;
  case 8: return 4;
  case 15: return 5;
  }
  return 0;
}

//== Tariq AAs --------------------------------
//returns 0 (empty '-') if a wrong chr is entered
//Defined in .h
//virtual count_t complement (count_t b) const { return b; }
//virtual count_t bitToIndex (count_t b) const { return b; }


//===== 5 bit representation
unsigned int AA8_BaseRep::baseToBit (char c) const {
	switch (c) {
		case '-':
			return 0;
		case 'A':
		case 'a':
			return 1;
		case 'L':
		case 'l':
			return 2;
		case 'R':
		case 'r':
			return 3;
		case 'K':
		case 'k':
			return 4;
		case 'N':
		case 'n':
			return 5;
		case 'M':
		case 'm':
			return 6;
		case 'D':
		case 'd':
			return 7;
		case 'F':
		case 'f':
			return 8;
		case 'C':
		case 'c':
			return 9;
		case 'P':
		case 'p':
			return 10;
		case 'Q':
		case 'q':
			return 11;
		case 'S':
		case 's':
			return 12;
		case 'E':
		case 'e':
			return 13;
		case 'T':
		case 't':
			return 14;
		case 'G':
		case 'g':
			return 15;
		case 'W':
		case 'w':
			return 16;
		case 'H':
		case 'h':
			return 17;
		case 'Y':
		case 'y':
			return 18;
		case 'I':
		case 'i':
			return 19;
		case 'V':
		case 'v':
			return 20;
    case 'U':
    case 'u':
      return 21;
    case 'O':
    case 'o':
      return 22;
		case '*':
			return 23;
	}
  return 0;
}

char AA8_BaseRep::bitToBase (seq_t b) const {
	switch (b) {
		case 0: return '-';
		case 1: return 'A';
		case 2: return 'L';
		case 3: return 'R';
		case 4: return 'K';
		case 5: return 'N';
		case 6: return 'M';
		case 7: return 'D';
		case 8: return 'F';
		case 9: return 'C';
		case 10: return 'P';
		case 11: return 'Q';
		case 12: return 'S';
		case 13: return 'E';
		case 14: return 'T';
		case 15: return 'G';
		case 16: return 'W';
		case 17: return 'H';
		case 18: return 'Y';
		case 19: return 'I';
		case 20: return 'V';
    case 21: return 'U';
    case 22: return 'O';
		case 23: return '*';
	}

  return '?';
}


