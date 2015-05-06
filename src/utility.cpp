#include<iostream>
#include <string>

#include<seqan/bam_io.h>
#include<seqan/seq_io.h>

using namespace seqan;

#include "utility.h"


int Reader::read_sam(std::string filename) {
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
int Reader::read_fasta(std::string filename, 
               StringSet<CharString>* ids, 
               StringSet<Dna5String>* seqs) {
    
    try {
        SeqFileIn seqFileIn(filename.c_str());
        readRecords(*ids, *seqs, seqFileIn);
    } catch (Exception const & e) {
        std::cout << "ERROR:" << e.what() << std::endl;
        return 1;
    }

   return 0;

}

/*
    Reads and returns fastq file sequences.
*/
int Reader::read_fastq(std::string filename,
               StringSet<CharString>* ids,
               StringSet<Dna5String>* seqs,
               StringSet<CharString>* phred) {
    try {
        SeqFileIn seqFileIn(filename.c_str());
        readRecords(*ids, *seqs, *phred, seqFileIn);
    } catch (Exception const & e) {
        std::cout << "ERROR:" << e.what() << std::endl;
        return 1;
    }

    return 0;
}
