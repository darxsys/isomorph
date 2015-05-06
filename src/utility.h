#ifndef UTILITY_H
#define UTILITY_H

#include<string>
#include<seqan/seq_io.h>

class Reader {
    public:
    int read_sam(std::string filename); 
    int read_fasta(std::string filename,
                   StringSet<CharString>* ids,
                   StringSet<Dna5String>* seqs);
    int read_fastq(std::string filename,
                   StringSet<CharString>* ids,
                   StringSet<Dna5String>* seqs,
                   StringSet<CharString>* phred);
};

#endif // UTILITY_H
