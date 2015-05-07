#ifndef COUNT_ESTIMATOR_H
#define COUNT_ESTIMATOR_H

#include "utility.h"
#include "estimator.h"

namespace isomorph {
    class CountEstimator : public Estimator {
    public:
        virtual void estimate_abundances(CharString left_pairs, CharString right_pairs,
                                         CharString transcripts);
    };
}

#endif // COUNT_ESTIMATOR_H