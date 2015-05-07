#ifndef UTILITY_H
#define UTILITY_H

#include<string>
#include<seqan/seq_io.h>

namespace isomorph {
    struct FastQData {
        StringSet<CharString> ids;
        StringSet<Dna5String> seqs;
        StringSet<CharString> phred;
    };

    struct FastAData {
        StringSet<CharString> ids;
        StringSet<Dna5String> seqs;
    };

    std::string execute_command(const char* cmd);

    class Reader {
    public:
        int read_sam(std::string filename); 
        int read_fasta(CharString filename,
                       FastAData* data);
        int read_fastq(CharString filename,
                       FastQData* data);
    };
}

#endif // UTILITY_H
