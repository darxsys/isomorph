#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

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

    struct SamData {
        BamHeader header;
        std::vector<BamAlignmentRecord> records;
    };

    std::string execute_command(const char* cmd);

    class Reader {
    public:
        int read_sam(CharString filename,
                     SamData* data); 
        int read_fasta(CharString filename,
                       FastAData* data);
        int read_fastq(CharString filename,
                       FastQData* data);
    };
}

#endif // UTILITY_H
