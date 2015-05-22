#ifndef COUNT_ESTIMATOR_H
#define COUNT_ESTIMATOR_H

#include "utility.h"
#include "estimator.h"

namespace isomorph {
    class CountEstimator : public Estimator {
    public:
        virtual void estimate_abundances(seqan::CharString reads, seqan::CharString transcripts,
                                         seqan::CharString pairs="");
    };
}

#endif // COUNT_ESTIMATOR_H
