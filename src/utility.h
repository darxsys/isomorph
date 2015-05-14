#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

namespace isomorph {
    struct FastQData {
        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;
        seqan::StringSet<seqan::CharString> phred;
    };

    struct FastAData {
        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;
    };

    struct SamData {
        seqan::BamHeader header;
        std::vector<seqan::BamAlignmentRecord> records;
    };

    std::string execute_command(const char* cmd);
    /*
        Prints all the important attributes of sam alignment records.
    */
    void print_sam_alignment_records(
            const std::vector<seqan::BamAlignmentRecord>& records);
    
    class Reader {
    public:
        int read_sam(seqan::CharString filename,
                     SamData* data); 
        int read_fasta(seqan::CharString filename,
                       FastAData* data);
        int read_fastq(seqan::CharString filename,
                       FastQData* data);
    };
}

#endif // UTILITY_H
