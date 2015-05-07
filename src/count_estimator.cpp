#include <iostream>
#include <vector>
#include <unordered_set>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace std;

#include "utility.h"
#include "count_estimator.h"
#include "utility.h"

void isomorph::CountEstimator::estimate_abundances(CharString left_pairs,
                                                   CharString right_pairs,
                                                   CharString transcripts) {

    isomorph::Reader reader;
    isomorph::FastQData left_data;
    isomorph::FastQData right_data;
    isomorph::FastAData transcript_data;

    reader.read_fastq(left_pairs, &left_data);
    reader.read_fastq(right_pairs, &right_data);
    reader.read_fasta(transcripts, &transcript_data);

    string transcripts_str(toCString(transcripts)); 
    string left_pairs_str(toCString(left_pairs));
    string right_pairs_str(toCString(right_pairs));

    string dir = "bowtie-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + dir + "/isomorph-index";
    execute_command(command.c_str());

    // runs the alignment
    command = "bowtie2 -x " + dir + "/isomorph-index -1 " + left_pairs_str + " -2 " +
              toCString(right_pairs_str) + " -S " + dir + "/isomorph.sam";
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam(dir + "/isomorph.sam");
    reader.read_sam(sam, &sam_data);

    // calculates the counts
    auto ids = transcript_data.ids;
    vector<long long> counts(length(ids), 0);
    auto records = sam_data.records;
    unordered_set<string> used_reads;

    for (int j = 0; j < length(records); ++j) {
        if (records[j].rID >= 0 and 
            used_reads.find(toCString(records[j].qName)) == used_reads.end()) {
            counts[records[j].rID]++;
            used_reads.insert(toCString(records[j].qName));
        }
    }

    for (int i = 0; i < counts.size(); ++i) {
        cout << ids[i] << '\n' << counts[i] << endl;
    }

    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
}