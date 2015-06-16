/** @file em_estimator.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares the em estimator class.
@details The em estimator header class
is declared here. This is the core header of isomorph. 
@date Tuesday, June 16, 2015
*/
#ifndef RSEM_ESTIMATOR_H
#define RSEM_ESTIMATOR_H

#include <vector>
#include <unordered_map>

#include "utility.h"
#include "estimator.h"
#include "read.h"
#include "single_read.h"
#include "paired_read.h"

namespace isomorph {
   /**
    *   EMEstimator class. \n
        This class estimates transcript abundances using a robust statistical model \n
        by applying the EM algorithm to obtain maximum likelihood estimates \n
        of the model parameters of the model.
    */
    class EMEstimator : public Estimator {
    public:
       /**
        *   Implementation of estimate_abundances virtual method. \n
            The method estimates transcript abundances according to the strategy of this class. \n
            @param reads path to the reads fastq file
            @param transcripts path to the reconstructed transcripts fasta file
            @param pairs path to the paired end file, if paired reads are used
        */      
       virtual void estimate_abundances(seqan::CharString left_pairs, 
                                        seqan::CharString right_pairs,
                                        seqan::CharString transcripts);
    private:
        struct EMParams {
            std::unordered_map<std::string, int> qNameToID;
            FastQData reads;
            FastQData pairs;
            FastAData transcripts;
            std::vector<std::unique_ptr<SingleRead> > single_reads;
            std::vector<std::unique_ptr<PairedRead> > paired_reads;
            int eff_num_reads;
            bool paired_end;

            // insert size normal dist params
            double insert_mean;
            double insert_stdev;
            
            EMParams() {
                eff_num_reads = 0;
                paired_end = false;
                insert_mean = -1;
                insert_stdev = 0;
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
                                
        void create_paired_end(const SamData& alignments,
                               EMParams& params);
                               
        void create_single_end(const SamData& alignments,
                               EMParams& params);                               
                                
        void output_result(const FastAData& transcripts, 
                           const EMResult& result, 
                           const std::string filename);
    };
}

#endif // RSEM_ESTIMATOR_H
