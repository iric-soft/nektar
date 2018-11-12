//                              -*- Mode: C++ -*-
// BaseRep.h
//

// Base: character representation of a base
// Bit:  bit representation of a base (least significant bits of the word)
// Index: bases are represented by contiguous numbers: 0, 1, 2, ... , n  Should be same order as bit.

// bitPerBase must be a power of 2

#ifndef __BaseRep_h_
#define __BaseRep_h_

#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;

#include "types.h"

class BaseRep {
public:

  unsigned int bitPerBase;
  count_t alphabetSize;

  unsigned int logBitPerBase;
  unsigned int basePerWord;
  unsigned int offsetShift;  // log_2 (basePerWord) must be integer
  seq_t bitMask;
  position_t offsetMask;

  // SPECIFICS ------------------------------------------------------------

  BaseRep () {}
  virtual ~BaseRep () {}

  virtual unsigned int baseToBit (char c) const = 0;
  virtual char bitToBase (seq_t b) const = 0;
  virtual count_t complement (count_t b) const = 0;
  virtual count_t bitToIndex (count_t b) const = 0;

  // GENERIC --------------------------------------------------------------

  void finalize () {
    logBitPerBase = (unsigned int)log2 (bitPerBase);
    basePerWord = ((sizeof (seq_t) * 8) / bitPerBase);
    offsetShift = 6 - logBitPerBase;  // log2 (64 / bitPerBase)
    offsetMask = ((position_t)1 << offsetShift) - 1;
    bitMask = (((seq_t)1 << bitPerBase) - 1);
  }

  count_t getBit (seq_t *base, position_t offset, bool comp) const ;
  char getBase (seq_t *base, position_t offset, bool comp) const ;

  void setBit (seq_t *base, position_t offset, bool comp, seq_t b);
  void setBase (seq_t *base, position_t offset, bool comp, char c);

  seq_t word_complement (seq_t w) const ;

  count_t countBit (seq_t *base, position_t offset, seq_t b) const ;

};

class Nuc4_BaseRep : public BaseRep {
  // .: 0000 : 0
  // A: 0001 : 1
  // C: 0010 : 2
  // G: 0100 : 4
  // T: 1000 : 8
  // N: 1111 : 15
public:

  // SPECIFICS ------------------------------------------------------------

  Nuc4_BaseRep () {
    bitPerBase = 4;
    alphabetSize = 6;
    finalize ();
  }

  virtual unsigned int baseToBit (char c) const ;
  virtual char bitToBase (seq_t b) const ;
  virtual count_t complement (count_t b) const ;
  virtual count_t bitToIndex (count_t b) const ;

};

class Qual_BaseRep : public BaseRep {
public:

  // SPECIFICS ------------------------------------------------------------

  Qual_BaseRep () {
    bitPerBase = 8;
    alphabetSize = 94;
    finalize ();
  }

  virtual unsigned int baseToBit (char c) const { return (unsigned int)c - 33; }
  virtual char bitToBase (seq_t b) const { return (char)(b + 33); }
  virtual count_t complement (count_t b) const { return b; }
  virtual count_t bitToIndex (count_t b) const { return b; }

};

//===== Tariq
class AA8_BaseRep : public BaseRep {//le nombre de bit doit un multiple de 2 (cf. sebastien)
//-: 000000, 
//A: 000001,
//L: 000010,
//R: 000011,
//K: 000100,
//N: 000101,
//M: 000110,
//D: 000111,
//F: 001000,
//C: 001001,
//P: 001010,
//Q: 001011,
//S: 001100,
//E: 001101,
//T: 001110,
//G: 001111,
//W: 010000,
//H: 010001,
//Y: 010010,
//I: 010011,
//V: 010100, 
//U: 010101, 
//O: 010110, 
//*: 010111
public:

  // SPECIFICS ------------------------------------------------------------

  AA8_BaseRep () {
    bitPerBase = 8;
    alphabetSize = 24; //+empty +stop
    finalize ();
  }

  virtual unsigned int baseToBit (char c) const ;
  virtual char bitToBase (seq_t b) const ;
  virtual count_t complement (count_t b) const { return b; }
  virtual count_t bitToIndex (count_t b) const { return b; }

};
#endif
