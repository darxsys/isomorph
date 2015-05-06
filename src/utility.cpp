#include<iostream>
#include <string>

#include<seqan/bam_io.h>

using namespace seqan;

class Reader {
    public:
    void read_sam_file(std::string filename) {
        BamFileIn samFile(filename.c_str());
        BamHeader header;
        
        BamFileOut bamFileOut;
        open(bamFileOut, "example.bam");

        readHeader(header, samFile);
        writeHeader(bamFileOut, header);
        BamAlignmentRecord record;

        while (!atEnd(samFile)) {
            readRecord(record, samFile);
            writeRecord(bamFileOut, record);
        }

        return;
    }
};
