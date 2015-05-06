#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

#include "utility.h"

int main(int argc, char** argv) {
    ArgumentParser parser("isomorph");

    addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, ArgParseOption(
        "L", "left", "Left read pairs.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "R", "right", "Right read pairs.",
        ArgParseArgument::STRING, "STRING"));
    addOption(parser, ArgParseOption(
        "T", "transcripts", "Transcript file."));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK) {
        cout << "ERROR while parsing arguments." << endl;
        return res == ArgumentParser::PARSE_ERROR;
    }

    CharString left_file = "";
    CharString right_file = "";
    CharString transcripts = "";
    getOptionValue(left_file, parser, "left");
    getOptionValue(right_file, parser, "right");
    getOptionValue(transcripts, parser, "transcripts");

    cout << "left_file \t" << left_file << '\n'
         << "right file \t" << right_file << '\n'
         << "transcripts \t" << transcripts << endl;

    return 0;
}

