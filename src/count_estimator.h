/** @file count_estimator.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares the count estimator class.
@details The count estimator header class
is declared here. 
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
   /**
    *   CountEstimator class. \n
        This class estimates transcript abundances using very simple \n
        ideas of counting number of reads mapped to a transcript as well as \n
        the number of bases mapped to a certain isoform and reports the result.
    */
    class CountEstimator : public Estimator {
    public:
       /**
        *   Implementation of estimate_abundances virtual method. \n
            The method estimates transcript abundances according to the strategy of this class. \n
            @param reads path to the reads fastq file
            @param transcripts path to the reconstructed transcripts fasta file
            @param pairs path to the paired end file, if paired reads are used
        */        
        virtual void estimate_abundances(seqan::CharString reads, 
                                         seqan::CharString transcripts,
                                         seqan::CharString pairs="");
                                         
    private:
        /**
         * The struct encapsulates all the parameters used by count algorithm and passed to it from the outside.
         */    
        struct AlgoParams {
            /**
             * Serves to map read names collected from fastq files to corresponding integer indices of reads.
             */            
            std::unordered_map<std::string, int> qNameToID;
            /**
             * Reads from file. Encapsulates seqan data structures.
             */            
            FastQData reads;
            /**
             * If paired-end is used, reads coming from the second file are stored here.
             */            
            FastQData pairs;
            /**
             * Stores transcript information.
             */            
            FastAData transcripts;
            /**
             * The alignments of the reads to reconstructed transcripts.
             */            
            SamData alignments;
            /**
             * If single reads are used, this is where their information regarding \n
               mappings will be stored for EM processing.
             */            
            std::vector<std::unique_ptr<SingleRead> > single_reads;
            /**
             * If paired reads are used, this is where their information regarding \n
               mappings will be stored for EM processing.
             */                
            std::vector<std::unique_ptr<PairedRead> > paired_reads;
            /**
             * This is the number of reads that have a mapping and are not ignored. \n
               This will be gone in future versions when single_reads and paired_reads \n
               will be containing only information for valid reads.
             */            
            int eff_num_reads;
            /**
             * If paired-end is used, this is set to true.
             */            
            bool paired_end;
        };

        /**
         * The resulting &tau; values are stored here and forwarded to postprocessing (if done) and output.
         */        
        struct CountResult {
            std::vector<double> counts;
        };

        /**
         * Does rudimentary data preprocessing.
         */    
        void preprocess_data(const seqan::CharString& transcripts,
                             const seqan::CharString& reads,
                             const seqan::CharString& pairs,
                             const std::string& output_dir,
                             AlgoParams& params);
        
        /**
         * Calculates the count of reads mapped to each transcript. \n
         For multi mapping reads, only the first mapping is considered valid.
         */                             
        void calculate_read_count(const AlgoParams& params,
                                  CountResult& result);
        /**
         * Outputs the result to a file called isomorph.abundances.fasta.
         */                                  
        void output_result(const FastAData& transcripts,
                           const CountResult& result,
                           const std::string& output_file);                                       
    };
}

#endif // COUNT_ESTIMATOR_H
