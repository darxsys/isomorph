/** @file utility.cpp
@author Pavlovic:Dario
@version Revision 0.2
@brief Utility data and functions used by other classes are implemented here.
@date Tuesday, June 16, 2015
*/

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include <ctype.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "./utility.h"

using namespace std;
using namespace seqan;

/*
    http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
*/
std::string isomorph::execute_command(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}

void isomorph::estimate_insert_size(const SamData& alignments,
                          std::pair<double, double>& params) {
    double mean = 0;
    int num_aligns = 0;

    for (auto record : alignments.records) {
        if (hasFlagAllProper(record) && hasFlagFirst(record)) {
            mean += record.tLen;
            num_aligns++;
        }
    }

    mean /= num_aligns;
    double var = 0;
    for (auto record : alignments.records) {
        if (hasFlagAllProper(record) && hasFlagFirst(record)) {
            int a = record.tLen - mean;
            var += a * a;
        }
    }

    var /= (num_aligns-1);
    params.first = mean;
    params.second = sqrt(var);
}

void isomorph::run_alignment(const string& reads,
                             const string& pairs,
                             const string& transcripts,
                             const string& aligner_path,
                             const string& output_dir,
                             const bool paired_end) {
    string aligner;
    if (aligner_path == "") {
        aligner = "bowtie2";
    } else {
        aligner = aligner_path;
    }

    // builds bowtie index
    string command = "bowtie2-build -f " + transcripts + " " +
                      output_dir + "/isomorph-bowtie-index";
    execute_command(command.c_str());

    // runs the alignment
    // parameters are set to be the same as in RSEM with bowtie2
    if (paired_end) {
        command = aligner + " -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
                  --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant \
                  -p 1 -k 200 -x " + output_dir + "/isomorph-bowtie-index -1 " + reads + " -2 " +
                  pairs + " -S " + output_dir + "/isomorph.sam";
    } else {
        command = aligner + " -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
                   --score-min L,0,-0.1 -p 1 -k 200 -x " + output_dir + "/isomorph-bowtie-index -U " + reads +
                   " -S " + output_dir + "/isomorph.sam";
    }

    execute_command(command.c_str());
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
