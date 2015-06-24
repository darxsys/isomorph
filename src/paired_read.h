/** @file paired_read.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares paired_end read header.
@details Class modeling paired-end reads is declared here. 
@date Tuesday, June 16, 2015
*/

#ifndef PAIRED_READ_H
#define PAIRED_READ_H

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include "read.h"

namespace isomorph {
   /**
    *   Paired read class. \n
        All the data necessary for paired read analysis is placed here.
    */ 	
	struct PairedRead : public Read {
		public:
		/**
			Read id given to it by isomorph for easier indexing.
		*/
		int id;
		/**
			The full read ID string from first fastq file.
		*/
		seqan::CharString left_fa_id;
		/**
			Sequence of mate from the first fastq file.
		*/
		seqan::Dna5String left_seq;
		/**
			Corresponding phred quality symbols..
		*/
		seqan::CharString left_phred;
		
		/**
			The full ID string of the second mate.
		*/
		seqan::CharString right_fa_id;
		/**
			Sequence of the second mate.
		*/
		seqan::Dna5String right_seq;
		/**
			Corresponding phred quality symbols.
		*/
		seqan::CharString right_phred;
		
		/**
			The reduced space of transcripts this read maps well to \n
			together with mapped read sequence from SAM file.
		*/
		std::vector<std::tuple<int, 
							   int, 
							   int, 
							   int, 
							   seqan::CharString, 
							   seqan::CharString> > pi_x_n;				
	};
}

#endif // PAIRED_READ_H