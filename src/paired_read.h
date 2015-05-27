#ifndef PAIRED_READ_H
#define PAIRED_READ_H

#include "read.h"

namespace isomorph {
	struct PairedRead : public Read {
		public:
		int id;
		// full id from fastq file
		seqan::CharString left_fa_id;
		seqan::Dna5String left_seq;
		seqan::CharString left_phred;
		
		seqan::CharString right_fa_id;
		seqan::Dna5String right_seq;
		seqan::CharString right_phred;
		
		// all the transcript alignments that are good
		// transcript id, left_position, right_position
		std::vector<std::tuple<int, int, int> > pi_x_n;				
	}
}

#endif // PAIRED_READ_H