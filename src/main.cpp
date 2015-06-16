/** @file main.cpp
@author Pavlovic:Dario
@version Revision 0.2
@brief Parses input parameters and runs the approriate estimator.
@details This file starts the program, parses input flags by
using seqan's arg_parse module and runs the selected estimator
on the input.
@date Tuesday, June 16, 2015
*/

/**
@mainpage
Isomorph is a software for RNA transcript abundance estimation.
It works with de novo assembled transcripts and supports
different modes of operation. \n \n
At its core, it contains
an EM algorithm that tries to optimize and find
parameters of a robust statistical model
proposed in the papers \n
"RNA-Seq gene expression estimation with read mapping uncertainty"\n
http://bioinformatics.oxfordjournals.org/content/26/4/493.fullm\n and
"RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome"\n
http://www.biomedcentral.com/1471-2105/12/323.

Isomorph, besides the above described statistical model, \n
also supports simple read count models described in the associated master thesis.

It requires seqan (http://www.seqan.de/), an open source C++ library for more efficient handling \n
of input and output data.
*/

#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "count_estimator.h"
#include "utility.h"
#include "em_estimator.h"

using namespace std;
using namespace seqan;
using namespace isomorph;

int main(int argc, char** argv) {
    ArgumentParser parser("isomorph");

    addOption(parser, ArgParseOption(
        "S", "single-end", "Use single-end reads. DEFAULT=false",
        ArgParseArgument::INTEGER, "INTEGER"));
    addOption(parser, ArgParseOption(
        "R", "reads", "Reads file.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "P", "pairs", "Read pairs, if applicable.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "T", "transcripts", "Transcript file.",
        ArgParseArgument::STRING, "STRING"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK) {
        cout << "ERROR while parsing arguments." << endl;
        return res == ArgumentParser::PARSE_ERROR;
    }

    CharString read_file = "";
    CharString pairs_file = "";
    CharString transcripts = "";
    int use_single_end = 0;
    
    getOptionValue(use_single_end, parser, "single-end");
    getOptionValue(read_file, parser, "reads");
    if (use_single_end == 0) {
        getOptionValue(pairs_file, parser, "pairs");
    }
    getOptionValue(transcripts, parser, "transcripts");

    cout << "left_file \t" << read_file << '\n'
         << "right file \t" << pairs_file << '\n'
         << "transcripts \t" << transcripts << endl;

    EMEstimator estimator;
    estimator.estimate_abundances(read_file, transcripts, pairs_file);

    return 0;
}
