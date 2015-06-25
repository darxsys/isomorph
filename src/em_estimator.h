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
       virtual void estimate_abundances(seqan::CharString reads,
                                        seqan::CharString transcripts,
                                        seqan::CharString aligner_path,
                                        seqan::CharString pairs="");
    private:
        /**
         * The struct encapsulates all the parameters used by EM algorithm and passed to it from the outside.
         */
        struct EMParams {
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
            /**
             * When paired-end is used, estimate of the mean of the insert sizes is stored here.
             */
            double insert_mean;
            /**
             * When paired-end is used, estimate of the standard deviation of the insert sizes is stored here.
             */
            double insert_stdev;

            /**
             * The constructor that sets some of the information to invalid values.
             */
            EMParams() {
                eff_num_reads = 0;
                paired_end = false;
                insert_mean = -1;
                insert_stdev = 0;
            }
        };

        /**
         * The resulting &tau; values are stored here and forwarded to postprocessing (if done) and output.
         */
        struct EMResult {
            std::vector<double> relative_expressions;
        };

        /**
            Method does preprocesing of data, elimination of extra stuff and creates
            necessary data structures for the EM to operate on. Most
            of the ideas are described in the original RSEM paper as
            well as in the follow-up.
        */
        void preprocess_data(const seqan::CharString& transcripts,
                             const seqan::CharString& reads,
                             const seqan::CharString& pairs,
                             const seqan::CharString& aligner_path,
                             const std::string& output_dir,
                             EMParams& params);
        /**
         * The core EM algorithm is done here.
         */
        void EMAlgorithm(EMParams& params, EMResult& result);

        /**
         * The method precalculates some probabilities that don't need to be recalculated anymore.
         */
        void precalc_posteriors(const EMParams& params,
                                std::vector<std::vector<double> >& posteriors);

        /**
         * Parses the SAM mapping file and creates information for paired reads and their mappings.
         */
        void create_paired_end(const SamData& alignments,
                               EMParams& params);

        /**
         * Parses the SAM mapping file and creates information for single reads and their mappings.
         */
        void create_single_end(const SamData& alignments,
                               EMParams& params);

        /**
         * Outputs the EM algorithm result. Output format is: \n
         >transcript_name \n
         &tau; \tab TPM_value
         */
        void output_result(const FastAData& transcripts,
                           const EMResult& result,
                           const std::string filename);
    };
}

#endif // RSEM_ESTIMATOR_H
