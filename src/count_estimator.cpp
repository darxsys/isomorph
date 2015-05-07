#include <iostream>

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

    string dir = "mkdir -p bowtie-tmp";
    execute_command(dir.c_str());

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + dir + "/isomorph-index";
    execute_command(command.c_str());

    // runs the alignment
    command = "bowtie2 -x " + dir + "/isomorph-index -1 " + left_pairs_str + " -2 " +
              toCString(right_pairs_str) + " -S " + dir + "/isomorph.sam";
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam("isomorph.sam");
    reader.read_sam(sam, &sam_data);

    // calculates the counts
    auto ids = transcript_data.ids;
    for (int i = 0; i < length(ids); ++i) {
        auto id = ids[i];
        long long count = 0;

        auto records = sam_data.records;
        for (int j = 0; j < length(records); ++j) {
            if (records[j].rID == i) {
                count++;
            }
        }

        cout << id << '\n' << count << endl;
    }

    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
}