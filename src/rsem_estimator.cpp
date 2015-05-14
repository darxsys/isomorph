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
        command = "bowtie2 -x " + dir + "/isomorph-index -r " + left_pairs_str +
                " -S " + dir + "/isomorph.sam";
    }
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam(dir + "/isomorph.sam");
    reader.read_sam(sam, &sam_data);
    
    print_sam_alignment_records(sam_data.records);
    return;
}

void isomorph::RsemEstimator::preprocess_data() {

}


