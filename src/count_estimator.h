/** @file count_estimator.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares the count estimator class.
@date Tuesday, June 16, 2015
*/

#ifndef COUNT_ESTIMATOR_H
#define COUNT_ESTIMATOR_H

#include <vector>
#include <unordered_map>

#include "utility.h"
#include "estimator.h"
#include "single_read.h"
#include "paired_read.h"

namespace isomorph {
    
    class CountEstimator : public Estimator {
    public:
        virtual void estimate_abundances(seqan::CharString reads, 
                                         seqan::CharString transcripts,
                                         seqan::CharString pairs="");
                                         
    private:
        struct AlgoParams {
            std::unordered_map<std::string, int> qNameToID;
            FastQData reads;
            FastQData pairs;
            FastAData transcripts;
            SamData alignments;
            std::vector<std::unique_ptr<SingleRead> > single_reads;
            std::vector<std::unique_ptr<PairedRead> > paired_reads;
            int eff_num_reads;
            bool paired_end;
        };
        
        struct CountResult {
            std::vector<double> counts;
        };
    
        void preprocess_data(const seqan::CharString& transcripts,
                             const seqan::CharString& reads,
                             const seqan::CharString& pairs,
                             const std::string& output_dir,
                             AlgoParams& params);
                             
        void calculate_read_count(const AlgoParams& params,
                                  CountResult& result);
                                  
        void calculate_base_count(const AlgoParams& params,
                                  CountResult& result);
                                  
        void output_result(const FastAData& transcripts,
                           const CountResult& result,
                           const std::string& output_file);                             
                                             
    };
}

#endif // COUNT_ESTIMATOR_H
