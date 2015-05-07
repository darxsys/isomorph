#include <iostream>
#include <string>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

#include "utility.h"

/*
    http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
*/
std::string isomorph::execute_command(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;    
}

int isomorph::Reader::read_sam(std::string filename) {
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

    return 1;
}

/* 
    Reads and returns fasta file sequences.
*/
int isomorph::Reader::read_fasta(CharString filename, 
               FastAData* data) {
    
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(filename))) {
        std::cerr << "ERROR: Could not open the file " << filename << std::endl;
        return 1;
    }

    try {
        readRecords(data->ids, data->seqs, seqFileIn);
    } catch (Exception const & e) {
        std::cout << "ERROR:" << e.what() << std::endl;
        return 1;
    }

    return 0;

}

/*
    Reads and returns fastq file sequences.
*/
int isomorph::Reader::read_fastq(CharString filename,
               FastQData* data) {

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(filename))) {
        std::cerr << "ERROR: Could not open the file " << filename << std::endl;
    }

    try {
        readRecords(data->ids, data->seqs, data->phred, seqFileIn);
    } catch (Exception const & e) {
        std::cout << "ERROR:" << e.what() << std::endl;
        return 1;
    }

    return 0;
}
