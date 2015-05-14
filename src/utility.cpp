#include <iostream>
#include <string>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "utility.h"

using namespace seqan;

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

void isomorph::print_sam_alignment_records(
        const std::vector<BamAlignmentRecord>& records) {

    for (auto record : records) {
        std::cout << "###############" << std::endl;
        std::cout << "Name: " << record.qName << std::endl;
        std::cout << "Ref. ID: " << record.rID << std::endl;
        std::cout << "flag: " << record.flag << std::endl;
        std::cout << "rNextId: " << record.rNextId << std::endl;
        std::cout << "Sequence: " << record.seq << std::endl;
        std::cout << "###############" << std::endl;
    }

    return;
}

int isomorph::Reader::read_sam(CharString filename,
                               SamData* data) {
    BamFileIn samFileIn;
    if (!open(samFileIn, toCString(filename))) {
        std::cerr << "ERROR: Could not open the file " << filename << std::endl;
        return 1;
    }

    try {
        BamHeader header;
        
        readHeader(data->header, samFileIn);
        BamAlignmentRecord record;

        while (!atEnd(samFileIn)) {
            readRecord(record, samFileIn);
            data->records.emplace_back(record);
        }
    } catch (Exception const & e) {
        std::cout << "ERROR:" << e.what() << std::endl;
        return 1;
    }

    return 0;
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
