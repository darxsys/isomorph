#include <iostream>
#include <vector>
#include <unordered_set>

#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "rsem_estimator.h"
#include "utility.h"
#include "count_estimator.h"

using namespace seqan;
using namespace std;

void isomorph::RsemEstimator::estimate_abundances(CharString reads, CharString transcripts,
                                                  CharString pairs) {
    bool paired_end = pairs == CharString("") ? false : true;

    isomorph::Reader reader;
    isomorph::FastQData reads_data;
    isomorph::FastQData pairs_data;
    isomorph::FastAData transcripts_data;

    cout << reads << " " << transcripts << " " << pairs << endl;

    reader.read_fastq(reads, &reads_data);
    if (paired_end) {
        reader.read_fastq(pairs, &pairs_data);
    }

    reader.read_fasta(transcripts, &transcripts_data);

    string transcripts_str(toCString(transcripts)); 
    string reads_str(toCString(reads));
    string pairs_str(toCString(pairs));

    string dir = "bowtie-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + dir + "/isomorph-index";
    execute_command(command.c_str());

    // runs the alignment
    if (paired_end) {
        command = "bowtie2 -x " + dir + "/isomorph-index -1 " + reads_str + " -2 " +
                toCString(pairs_str) + " -S " + dir + "/isomorph.sam";
    } else {
        command = "bowtie2 --all -N 1 -L 25 -q -x " + dir + "/isomorph-index -U " + reads_str +
                " -S " + dir + "/isomorph.sam";
    }
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam(dir + "/isomorph.sam");
    reader.read_sam(sam, &sam_data);
    EMParams params;

    // print_sam_alignment_records(sam_data.records);
    if (paired_end) {
        preprocess_data(sam_data, 
                        transcripts_data, 
                        reads_data, 
                        pairs_data, 
                        params);
    } else {
        preprocess_data(sam_data, 
                        transcripts_data, 
                        reads_data, 
                        params);
    }

    // cleaning up
    // string rm = "rm -rf " + dir;
    // execute_command(rm.c_str());
    return;
}

/*
    Single-end data version of data preprocessing.

*/
void isomorph::RsemEstimator::preprocess_data(const SamData& alignments,
                                              const FastAData& transcripts, const FastQData& reads, 
                                              EMParams& params) {

    int num_reads = length(reads.ids);
    int num_transcripts = length(transcripts.ids);

    // this vector will also be used to indicate if a read has any good alignments.
    vector<short> reads_info(num_reads, 0);
    params.pi_x_n.insert(params.pi_x_n.begin(), 
                         num_transcripts+1,
                         reads_info);
    
    for (int i = 0; i < num_reads; ++i) {
        params.qNameToID[toCString(reads.ids[i])] = i;
    }

    // arange the alignments in the "neighboring matrix"
    for (auto record : alignments.records) {
        string qName = toCString(record.qName);
        int read_id = params.qNameToID[qName];

        if (record.rID != record.INVALID_REFID) {
            reads_info[read_id] = true;
            params.pi_x_n[record.rID][read_id] = 1;
        } 
    }
    
    // all the reads that have no good alignment are assigned to dummy isoform
    for (int i = 0; i < num_reads; ++i) {
        if (!reads_info[i]) {
            params.pi_x_n[num_transcripts][i] = 1;
        }
    }

    return;
}

/*
    Paired-end version of data preprocessing.
*/
void isomorph::RsemEstimator::preprocess_data(const SamData& alignments,
                                              const FastAData& transcripts, const FastQData& reads, 
                                              const FastQData& pairs, EMParams& params) {

}


