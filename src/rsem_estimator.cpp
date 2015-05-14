#include <iostream>
#include <vector>
#include <unordered_set>

#include <seqan/arg_parse.h>
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
        command = "bowtie2 -x " + dir + "/isomorph-index -r " + reads_str +
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
        preprocess_data(sam_data, transcripts_data, reads_data, pairs_data, params);
    } else {
        preprocess_data(sam_data, transcripts_data, reads_data, params);
    }

    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
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
    vector<short> v(num_reads, 0);
    params.pi_x_n.insert(v, num_transcripts);
    
    // filter out alignments with more than 5 mismatches
    for (auto record : alignments.records) {
         
    }


}

/*
    Paired-end version of data preprocessing.
*/
void isomorph::RsemEstimator::preprocess_data(const SamData& alignments,
                                              const FastAData& transcripts, const FastQData& reads, 
                                              const FastQData& pairs, EMParams& params) {

}


