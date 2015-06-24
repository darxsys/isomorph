/** @file utility.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Utility data and functions used by other classes are declared here. 
@date Tuesday, June 16, 2015
*/

#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>

#include <string>
#include <unordered_map>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

const double PI = 3.14159265358979;

namespace isomorph {
    /**
        FastQData struct.
        Encapsulates seqan fastq file data classes for easier handling.
    */
    struct FastQData {
        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;
        seqan::StringSet<seqan::CharString> phred;
    };

    /**
        FastAData struct.
        Encapsulates seqan fasta file data classes for easier handling.
    */
    struct FastAData {
        seqan::StringSet<seqan::CharString> ids;
        seqan::StringSet<seqan::Dna5String> seqs;
    };

    /**
        SamData struct.
        Encapsulates seqan fsam file data classes for easier handling.
    */
    struct SamData {
        seqan::BamHeader header;
        std::vector<seqan::BamAlignmentRecord> records;
    };
    
    /**
        Executes a command and returns the output.
        
        @param cmd string that represents the command.
    */
    std::string execute_command(const char* cmd);
    
    /**
        Calculates the probability of x according to a normal distribution.
        
        @params mean normal distribution mean
        @params stdev standard deviation
        @params x the argument of probability function
    */
    inline double prob_normal(double mean, double stdev, double x) {
        return 1 / stdev / sqrt(2 * PI) * exp(-((x-mean) * (x-mean)) / 2 / (stdev * stdev));
    }
    
    /**
        Estimates paired end insert sizes from valid alignments.
        
        @param alignments valid sam alignments for analysis.
        @param params normal distribution parameters which will be calculated.
    */
    void estimate_insert_size(const SamData& alignments,
                              std::pair<double, double>& params);

    /**
        Runs the alignment of reads to transcripst using a given mapping tool.\n
        For now, only bowtie2 is supported.
        
        @param reads path to the reads file
        @param pairs path to the paired mates file, if used
        @param transcripts path to the transcripts file
        @param output_dir a temporary directory to store bowtie index
        @param paired_end a flag denoting if pairs are used
    */                              
    void run_alignment(const std::string& reads,
                       const std::string& pairs,
                       const std::string& transcripts,
                       const std::string& output_dir,
                       const bool paired_end);

    /**
        Prints all the important attributes of sam alignment records. \n
        Mainly used for debugging.
        
        @param records list of BAM alignment record class from seqan.
    */
    void print_sam_alignment_records(
            const std::vector<seqan::BamAlignmentRecord>& records);
    
    /**
        Reader class that encapsulates methods for data input and output.
    */
    class Reader {
    public:
        /**
            Reads SAM alignment file.
            @param filename path to the file
            @param data Pointer to samdata struct that will be filled.
        */
        int read_sam(seqan::CharString filename,
                     SamData* data); 
                     
        /**
            Reads a fasta file.
            @param filename path to the file
            @param data pointer to fastadata struct that will be filled.
        */
        int read_fasta(seqan::CharString filename,
                       FastAData* data);
                       
        /**
            Reads a fastq file.
            @param filename path to the file
            @param data pointer to fastqadat struct that will be filled.
        */                       
        int read_fastq(seqan::CharString filename,
                       FastQData* data);
    };
    
    /**
        Calculates the reverse complement of a character.
        @param c character to be reversed.
    */
    inline char reverse_complement(const char& c) {
        switch(toupper(c)) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        }
        
        return '\0';
    }
}

#endif // UTILITY_H
