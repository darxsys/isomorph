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

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " isomorph-index";
    execute_command(command.c_str());

    // run alignment
    command = "bowtie2 -x isomorph-index -1 " + left_pairs_str + " -2 " +
              toCString(right_pairs_str) + " -S isomorph.sam";
    execute_command(command.c_str());

    // read sam data

}