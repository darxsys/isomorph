#ifndef RSEM_ESTIMATOR_H
#define RSEM_ESTIMATOR_H

#include "estimator.h"

namespace isomorph {
    class RsemEstimator : public Estimator {
    public:
       virtual void estimate_abundances(seqan::CharString left_pairs, seqan::CharString right_pairs,
                                        seqan::CharString transcripts);
    private:
        void preprocess_data();
    };
}

#endif // RSEM_ESTIMATOR_H
