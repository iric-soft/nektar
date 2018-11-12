# NEKTAR : NEw Kmer Tool for Analysis of Reads

This repository holds several tools to extract the NEKTAR of kmers created by [jellyfish](http://www.genome.umd.edu/jellyfish.html).

# Requirements #
* nektar.py:
* * gcc/4.9.3 jellyfish/2.2.3 python/2.7.6
* nektar_cpp:
* * boost/1.60.0 jellyfish/2.2.3
* * for kmer_assembly: ulimit -s unlimited

# Installation #
* Go in fmindex folder and execute the script ./install.sh
* Go in nektar_cpp folder and choose your makefile (in function of your OS) and execute it (with *make*)

# Help #
* python *path_to*/nektar.py -h
* python *path_to*/nektar.py [cmd] -h

# File format #
* project.csv:
* * Tabulated file with a header formatted as follow: SAMPLE,JF,GROUP

|SAMPLE |JF             |GROUP  |
| ------|:-------------:| -----:|
|X      |/path/kmers.jf |Query  |
|Y      |/path/kmers.jf |Ref    |

# Pipeline to run nektar #
* Create jellyfish file for each sample
* Create the project file
