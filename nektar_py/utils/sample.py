class Sample(object):

    def __init__(self, name, group, num_kmer, files):
        self._name = name
        self._group = group
        self._num_kmer = num_kmer

        self._jf_file = files["JF"]
        self._fq1_file = files["FASTQ1"]
        self._fq2_file = files["FASTQ2"]

    @property
    def name(self):
        return self._name

    @property
    def group(self):
        return self._group

    @property
    def num_kmer(self):
        return self._num_kmer

    @property
    def jf_file(self):
        return self._jf_file

    @property
    def fastq_files(self):
        return (self._fq1_file, self._fq2_file)
