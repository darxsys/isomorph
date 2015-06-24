/** @file single_read.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares single_end read header.
@details Class modeling single-end reads is declared here. 
@date Tuesday, June 16, 2015
*/

#ifndef SINGLE_READ_H
#define SINGLE_READ_H

#include "read.h"

namespace isomorph {
   /**
    *   A single-end read class. \n
        All the data necessary for single-read analysis is placed here.
    */ 		
	struct SingleRead : public Read {
		public:
		/**
			Read id given to it by isomorph for easier indexing.
		*/		
		int id;
		/**
			Full read ID from the fastq file.
		*/
		seqan::CharString fa_id;
		/**
			Read sequence from the fastq file.
		*/
		seqan::Dna5String seq;
		/**
			Associated phred quality symbols.
		*/
		seqan::CharString phred;
		
		/**
			The reduced space of the transcripts this read maps well to.
		*/
		std::vector<std::tuple<int, int> > pi_x_n;				
	};
}

#endif // SINGLE_READ_H