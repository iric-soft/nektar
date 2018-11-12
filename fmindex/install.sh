rm -r *.o

g++ -g -Wall -fPIC -c util.cc BaseRep.cc SeqIter.cc StrDict.cc SeqStore.cc SeqRange.cc MultiRange.cc FmIndex.cc FmIndex_wrap.cxx -I/soft/bioinfo/linux_RH6/python-2.7.6/include/python2.7/ -lz -std=c++0x -fopenmp -lm
g++ -shared util.o BaseRep.o SeqIter.o SeqStore.o SeqRange.o StrDict.o FmIndex.o FmIndex_wrap.o -o _FmIndex.so -lz -std=c++0x -fopenmp -lm
