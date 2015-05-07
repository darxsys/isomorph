#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <seqan/stream.h>
#include <seqan/basic.h>

/* 
    Abstract base class
*/
namespace isomorph {
    class Estimator {
    public:
        /*
            Estimates how many reads map to a transcript.
            Only takes into account the best alignment for each read.
        */
        virtual void estimate_abundances(CharString left_pairs, CharString right_pairs,
                                         CharString transcripts)=0;
    };
}

#endif // ESTIMATOR_H