/** @file estimator.h
@author Pavlovic:Dario
@version Revision 0.2
@brief Declares a generic abstract estimator base class \n 
inherited by other estimator classes.
@details This is the abstract base class header that other estimators \n
must include and inherit to keep everything according to the strategy pattern. 
@date Tuesday, June 16, 2015
*/

#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include <seqan/stream.h>
#include <seqan/basic.h>

namespace isomorph {
   /**
    *   Estimator class. \n
        Abstract base class all estimator classes must inherit.
    */    
    class Estimator {
    public:
        /**
            This method must be implemented in the classes inheriting this one.
            @param reads path to the reads fastq file
            @param transcripts path to the transcripts fasta file
            @param pairs path to the paired end fastq file, if used
            @param aligner_path path to the aligner used (only bowtie supported for now)
        */
        virtual void estimate_abundances(seqan::CharString reads, 
                                         seqan::CharString transcripts,
                                         seqan::CharString aligner_path,
                                         seqan::CharString pairs="")=0;
    };
}

#endif // ESTIMATOR_H
