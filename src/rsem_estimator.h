#ifndef RSEM_ESTIMATOR_H
#define RSEM_ESTIMATOR_H

#include <vector>
#include <unordered_map>

#include "utility.h"
#include "estimator.h"

namespace isomorph {
    struct EMParams {
        std::unordered_map<std::string, int> qNameToID;
        std::vector<std::vector<short> > pi_x_n;
        FastQData reads;
        FastAData transcripts;
    };

    struct EMResult {
        std::vector<double> relative_expressions;
    };

    class RsemEstimator : public Estimator {
    public:
       virtual void estimate_abundances(seqan::CharString left_pairs, 
                                        seqan::CharString right_pairs,
                                        seqan::CharString transcripts);
    private:
        /*
            Methods preprocess data, eliminate extra stuff and create
            necessary data structures for the EM to operate on. Most
            of the ideas are described in the original RSEM paper as
            well as in the follow-up.
        */
        void preprocess_data(const SamData& alignments,
                             const FastAData& transcripts, const FastQData& reads, 
                             EMParams& params);

        void preprocess_data(const SamData& alignments,
                             const FastAData& transcripts, const FastQData& reads, 
                             const FastQData& pairs, EMParams& params);
        
        void EMAlgorithm(EMParams& params, EMResult& result);
        void precalc_posteriors(const EMParams& params, 
                                std::vector<std::vector<double> >& posteriors);
        void output_result(EMResult& result, std::string filename);
    };
}

#endif // RSEM_ESTIMATOR_H
