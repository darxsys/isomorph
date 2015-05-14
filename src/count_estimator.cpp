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

void isomorph::CountEstimator::estimate_abundances(CharString reads,
                                                   CharString transcripts,
                                                   CharString pairs) {
    bool paired_end = pairs == "" ? false : true;
    isomorph::Reader reader;
    isomorph::FastQData reads_data;
    isomorph::FastQData pairs_data;
    isomorph::FastAData transcripts_data;

    reader.read_fastq(reads, &reads_data);
    if (paired_end) {
        reader.read_fastq(pairs, &pairs_data);
    }

    reader.read_fasta(transcripts, &transcripts_data);

    string transcripts_str(toCString(transcripts)); 
    string left_pairs_str(toCString(reads));
    string right_pairs_str(toCString(pairs));

    string dir = "bowtie-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + dir + "/isomorph-index";
    execute_command(command.c_str());

    // runs the alignment
    if (paired_end) {
        command = "bowtie2 -x " + dir + "/isomorph-index -1 " + left_pairs_str + " -2 " +
                toCString(right_pairs_str) + " -S " + dir + "/isomorph.sam";
    } else {
        command = "bowtie";
    }
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam(dir + "/isomorph.sam");
    reader.read_sam(sam, &sam_data);

    // calculates the counts
    auto ids = transcripts_data.ids;
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
