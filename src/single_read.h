#ifndef SINGLE_READ_H
#define SINGLE_READ_H

namespace isomorph {
	struct SingleRead : public Read {
		public:
		int id;
		// full id from fastq file
		seqan::CharString fa_id;
		seqan::Dna5String seq;
		seqan::CharString phred;
		
		// all the transcript alignments that are good
		// transcript id, alignment position
		std::vector<std::tuple<int, int> > pi_x_n;				
	};
}

#endif // SINGLE_READ_H