#ifndef RSEM_ESTIMATOR_H
#define RSEM_ESTIMATOR_H

#include <vector>
#include <unordered_map>

#include "utility.h"
#include "estimator.h"

namespace isomorph {

    class RsemEstimator : public Estimator {
    public:
       virtual void estimate_abundances(seqan::CharString left_pairs, 
                                        seqan::CharString right_pairs,
                                        seqan::CharString transcripts);
    private:
        struct EMParams {
            std::unordered_map<std::string, int> qNameToID;
            std::vector<std::vector<std::pair<int, int> > > pi_x_n;
            FastQData reads;
            FastAData transcripts;
            int eff_num_reads;
            
            EMParams() {
                eff_num_reads = 0;
            }
        };
    
        struct EMResult {
            std::vector<double> relative_expressions;
        };
        
        /*
            Method does preprocesing of data, elimination of extra stuff and creates
            necessary data structures for the EM to operate on. Most
            of the ideas are described in the original RSEM paper as
            well as in the follow-up.
        */
        void preprocess_data(const seqan::CharString& transcripts, 
                             const seqan::CharString& reads,
                             const seqan::CharString& pairs,
                             const std::string& output_dir, 
                             EMParams& params);
        
        void EMAlgorithm(EMParams& params, EMResult& result);
        void precalc_posteriors(const EMParams& params, 
                                std::vector<std::vector<double> >& posteriors);
                                
        void output_result(const FastAData& transcripts, 
                           const EMResult& result, 
                           const std::string filename);
    };
}

#endif // RSEM_ESTIMATOR_H
