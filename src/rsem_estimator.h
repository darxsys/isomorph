#ifndef RSEM_ESTIMATOR_H
#define RSEM_ESTIMATOR_H

#include <seqan/stream.h>
#include <seqan/basic.h>

#include "estimator.h"

namespace isomorph {
    class RsemEstimator : public Estimator {
    public:
       virtual void estimate_abundances(CharString left_pairs, CharString right_pairs,
                                        CharString transcripts);
    private:
        void preprocess_data();
    };
}

#endif //RSEM_ESTIMATOR_H
