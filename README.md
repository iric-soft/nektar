# NEKTAR : NEw Kmer Tool for Analysis of Reads

This repository holds several tools to extract the NEKTAR of kmers created by [jellyfish](http://www.genome.umd.edu/jellyfish.html).

# Requirements #
* gcc/4.9.3 jellyfish/2.2.3 python/2.7.6
* for kmer_assembly: ulimit -s unlimited
* [boost](https://www.boost.org/)

# Install #
* Go in fmindex folder and execute the script ./install.sh
* Go in nektar_cpp folder and choose your makefile (in function of your OS) and execute it (with *make*)

# Help #
* python *path_to*/nektar.py -h
* python *path_to*/nektar.py [cmd] -h
* *path_to*/kmer_assembly -h

# Usage #
kmer_assembly -k kmer.tsv -p design.tsv -o test.res 

* project.tsv: tabulated file with 3 columns
```
SAMPLE  JF      GROUP
14H124  path_to/14H124.jf  Ref
```
* kmer.tsv : tabulated file with 2 columns, which can be generate by `jellyfish dump -c -t path_to/14H124.jf`
```
kmer    14H124
GTTTCCTTCTACAGCATGTCAGCATCTCAAGTT       101
CTGATTCTCCAAGATTCCCTCATAGAGGAATTT       4
ATTGCAGAATACCAGCGTGTATTGCAGGAGAAC       6
```
