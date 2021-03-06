/** @file em_estimator.cpp
@author Pavlovic:Dario
@version Revision 0.2
@brief CountEstimator class methods are implemented here. \n
Detailed docs are available in the corresponding .h file.
@date Tuesday, June 16, 2015
*/

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "./count_estimator.h"
#include "./utility.h"

using namespace seqan;
using namespace std;

void isomorph::CountEstimator::estimate_abundances(CharString reads,
                                                   CharString transcripts,
                                                   CharString aligner_path,
                                                   CharString pairs) {
    string dir = "isomorph-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());
    AlgoParams params;
    CountResult result;

    preprocess_data(transcripts,
                    reads,
                    pairs,
                    aligner_path,
                    dir,
                    params);

    calculate_read_count(params, result);

    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
    return;
}

void isomorph::CountEstimator::preprocess_data(const CharString& transcripts,
                                               const CharString& reads,
                                               const CharString& pairs,
                                               const CharString& aligner_path,
                                               const string& output_dir,
                                               AlgoParams& params) {
    cerr << "Preprocessing data." << endl;
    bool paired_end = pairs == "" ? false : true;
    isomorph::Reader reader;

    reader.read_fastq(reads, &params.reads);
    if (paired_end) {
        reader.read_fastq(pairs, &params.pairs);
    }

    reader.read_fasta(transcripts, &params.transcripts);

    string transcripts_str(toCString(transcripts));
    string left_pairs_str(toCString(reads));
    string right_pairs_str(toCString(pairs));
    string aligner_str(toCString(aligner_path));

    run_alignment(left_pairs_str,
                  right_pairs_str,
                  transcripts_str,
                  aligner_str,
                  output_dir,
                  paired_end);

    // reads sam data
    CharString sam(output_dir + "/isomorph.sam");
    reader.read_sam(sam, &params.alignments);
    params.paired_end = paired_end;
    cerr << "Done preprocessing data." << endl;
    return;
}

void isomorph::CountEstimator::calculate_read_count(const AlgoParams& params,
                                                    CountResult& result) {
    cerr << "Calculating read counts." << endl;
    // calculates the counts
    auto ids = params.transcripts.ids;
    vector<int64_t> counts(length(ids), 0);
    auto records = params.alignments.records;
    unordered_set<string> used_reads;
    double counts_sum = 0;

    for (int j = 0; j < length(records); ++j) {
        if (records[j].rID != records[j].INVALID_REFID &&
                used_reads.find(toCString(records[j].qName)) ==
                used_reads.end()) {
            counts[records[j].rID]++;
            counts_sum += 1;
            used_reads.insert(toCString(records[j].qName));
        }
    }

    for (int i = 0; i < counts.size(); ++i) {
        result.counts.emplace_back(counts[i] / counts_sum);
    }

    cerr << "Done calculating." << endl;
    return;
}

void isomorph::CountEstimator::output_result(const FastAData& transcripts,
                                             const CountResult& result,
                                             const string& output_file) {
    cerr << "Outputing results." << endl;
    ofstream output;
    output.open(output_file.c_str(), ofstream::out | ofstream::trunc);

    int num_transcripts = result.counts.size();
    double sum = 0;
    // needed to transform to TPM according to RSEM paper
    double v_sum = 0;
    vector<double> tpm(num_transcripts, 0);

    for (int i = 0; i < num_transcripts; ++i) {
        sum += result.counts[i];
        // calculate TPM
        int t_len = seqan::length(transcripts.seqs[i]);
        tpm[i] = (result.counts[i] / t_len);
        v_sum += tpm[i];
    }

    double tpm_sum = 0;
    for (int i = 0; i < num_transcripts; ++i) {
        tpm[i] *= 1e6 / v_sum;
        tpm_sum += tpm[i];

        // first ni_i then tpm_i
        output << ">" << transcripts.ids[i] << endl;
        output << result.counts[i];
        output << '\t' << tpm[i] << endl;
    }

    cerr << "Expression sum: " << sum << endl;
    cerr << "TPM sum: " << tpm_sum << endl;
    output.close();
    cerr << "Done outputing results." << endl;
    return;
}
