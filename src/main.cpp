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

/**
 * Struct to group and encapsulate parsed arguments.
 */
struct Args {
    bool paired_end;
    CharString reads;
    CharString pairs;
    CharString transcripts;
    CharString bowtie_path;
    bool use_count;
};

/**
 * Parses command line arguments using seqan parsing API.
 */
int parse_arguments(int argc, char** argv, Args& args) {
    ArgumentParser parser("isomorph");

    addOption(parser, ArgParseOption(
        "S", "single-end", "Use single-end reads."));
    addOption(parser, ArgParseOption(
        "C", "count", "Use simple count instead of EM algorithm."));        
    addOption(parser, ArgParseOption(
        "R", "reads", "Reads file.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "P", "pairs", "Read pairs, if applicable.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "T", "transcripts", "Transcript file.",
        ArgParseArgument::STRING, "STRING"));
    
    addOption(parser, ArgParseOption(
        "B", "bowtie_path", "Path to bowtie executable",
        ArgParseArgument::STRING, "STRING"));

    setRequired(parser, "R");
    setRequired(parser, "T");

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK) {
        if (!isSet(parser, "help")) {
            cerr << "Error while parsing arguments." << endl;
        }
        
        return -1;
    }

    getOptionValue(args.reads, parser, "reads");
    
    args.paired_end = !isSet(parser, "single-end");
    if (args.paired_end) {
        getOptionValue(args.pairs, parser, "pairs");
    }
        
    getOptionValue(args.transcripts, parser, "transcripts");
    args.use_count = isSet(parser, "count"); 

    setDefaultValue(parser, "bowtie_path", "");
    getOptionValue(args.bowtie_path, parser, "bowtie_path");
    return 0;   
}

int main(int argc, char** argv) {
    Args args;
    int ret = parse_arguments(argc, argv, args);
    
    if (ret != 0) {
        return -1;
    }
    
    if (args.use_count) {
        CountEstimator estimator;
        estimator.estimate_abundances(args.reads, 
                                      args.transcripts, 
                                      args.pairs,
                                      args.bowtie_path);
    } else {
        EMEstimator estimator;
        estimator.estimate_abundances(args.reads,
                                      args.transcripts, 
                                      args.pairs,
                                      args.bowtie_path);
    }

    return 0;
}
