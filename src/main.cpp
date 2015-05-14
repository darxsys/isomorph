#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "count_estimator.h"
#include "utility.h"
#include "rsem_estimator.h"

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

    RsemEstimator estimator;
    estimator.estimate_abundances(read_file, transcripts, pairs_file);

    return 0;
}

