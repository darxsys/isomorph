#include <ctype.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>

#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "rsem_estimator.h"
#include "utility.h"
#include "count_estimator.h"
#include "read.h"
#include "single_read.h"
#include "paired_read.h"

using namespace seqan;
using namespace std;

void isomorph::RsemEstimator::estimate_abundances(CharString reads,
                                                  CharString transcripts,
                                                  CharString pairs) {
    
    string dir = "bowtie-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());    
    EMParams params;
    preprocess_data(transcripts, 
                    reads, 
                    pairs,
                    dir, 
                    params);

    EMResult result;
    EMAlgorithm(params, result);
    output_result(params.transcripts, result, "isomorph.abundances.fasta");
    
    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
    return;
}

void isomorph::RsemEstimator::preprocess_data(const CharString& transcripts, 
                                              const CharString& reads,
                                              const CharString& pairs,
                                              const string& output_dir, 
                                              EMParams& params) {
    
    cerr << "Preprocessing data." << endl;
    bool paired_end = pairs == CharString("") ? false : true;
    params.paired_end = paired_end;
    Reader reader;
    string transcripts_str(toCString(transcripts)); 
    string reads_str(toCString(reads));
    string pairs_str(toCString(pairs));

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + output_dir + "/isomorph-bowtie-index";
    execute_command(command.c_str());

    // runs the alignment
    // parameters are set to be the same as in RSEM with bowtie2
    if (paired_end) {
        command = "bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
                  --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant \
                  -p 1 -k 200 -x " + output_dir + "/isomorph-bowtie-index -1 " + reads_str + " -2 " +
                  toCString(pairs_str) + " -S " + output_dir + "/isomorph.sam";
    } else {
        command = "bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
                   --score-min L,0,-0.1 -p 1 -k 200 -x " + output_dir + "/isomorph-bowtie-index -U " + reads_str +
                   " -S " + output_dir + "/isomorph.sam";
    }
    
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData alignments;
    CharString sam(output_dir + "/isomorph.sam");
    reader.read_sam(sam, &alignments);                                                      
    reader.read_fastq(reads, &params.reads);
    reader.read_fasta(transcripts, &params.transcripts);

    int num_reads = length(params.reads.ids);
    int num_transcripts = length(params.transcripts.ids);

    if (paired_end) {
        FastQData pairs_data;
        reader.read_fastq(pairs, &params.pairs);
        create_paired_end(alignments, params);
    } else {
        create_single_end(alignments, params);
    }

//    vector<pair<int, int> > v;
//    params.pi_x_n.insert(params.pi_x_n.begin(), 
//                         num_reads,
//                         v);
//    
//    vector<bool> read_processed(num_reads, false);
//        
//    for (int i = 0; i < num_reads; ++i) {
//        string qName = toCString(params.reads.ids[i]);
//        int pos = qName.find(' ', 0);
//        
//        if (pos != string::npos) {
//            qName = qName.substr(0, pos);
//        }
//
//        params.qNameToID[qName] = i;
//    }
//    
//    // arange the alignments in the "neighboring matrix"
//    int ID = 0;
//    int count = 0;
//    for (auto record : alignments.records) {
//        string qName = toCString(record.qName);
//        int read_id = params.qNameToID[qName];
//
//        if (record.rID != record.INVALID_REFID) {
//            count++;
//            // reads that have at least one good alignment
//            // here, I rely on bowtie to provide no discordant alignments
//            if (!read_processed[read_id]) {
//                params.eff_num_reads++;
//            }
//            
//            read_processed[read_id] = true;
//            params.pi_x_n[read_id].emplace_back(record.rID, record.beginPos);
//        } 
//    }
    
//    cerr << "Number of reads: " << num_reads << endl;
//    cerr << "Effective number of reads: " << params.eff_num_reads << endl;
//    cerr << "Done preprocessing data." << endl;
    return;
}                   

void isomorph::RsemEstimator::create_single_end(const SamData& alignments,
                                                EMParams& params) {
    
    auto& reads = params.reads;
    auto& transcripts = params.transcripts;
    auto& processed_reads = params.single_reads;
                                                        
    // creating the read objects
    int num_reads = seqan::length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    vector<bool> has_alignment(num_reads, false);
    
    for (int i = 0; i < num_reads; ++i) {
        string qName = toCString(params.reads.ids[i]);
        int pos = qName.find(' ', 0);
        
        if (pos != string::npos) {
            qName = qName.substr(0, pos);
        }

        params.qNameToID[qName] = i;
        processed_reads.emplace_back(unique_ptr<SingleRead>(new SingleRead));
        unique_ptr<SingleRead>& read = processed_reads[i];
        read->id = i;
        
        read->fa_id = reads.ids[i];
        read->seq = reads.seqs[i];
        read->phred = reads.phred[i];
    }
    
    // filling up the pi_x_n vectors    
    int ID = 0;
    int count = 0;
    for (auto record : alignments.records) {
        string qName = toCString(record.qName);
        int read_id = params.qNameToID[qName];
        unique_ptr<SingleRead>& read = processed_reads[read_id];

        // only records that are completely mapped
        if (record.rID != record.INVALID_REFID) {
            if (!has_alignment[read_id]) {
                params.eff_num_reads++;
                has_alignment[read_id] = true;
            }
            
            read->pi_x_n.emplace_back(record.rID, record.beginPos);
        } 
    }
    
    return;  
}

void isomorph::RsemEstimator::create_paired_end(const SamData& alignments,
                                                EMParams& params) {

    auto& reads = params.reads;
    auto& pairs = params.pairs;
    auto& transcripts = params.transcripts;
    auto& processed_reads = params.paired_reads;

    // creating the read objects
    int num_reads = seqan::length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    vector<bool> has_alignment(num_reads, false);
    
    for (int i = 0; i < num_reads; ++i) {
        string qName = toCString(params.reads.ids[i]);
        int pos = qName.find(' ', 0);
        
        if (pos != string::npos) {
            qName = qName.substr(0, pos);
        }

        params.qNameToID[qName] = i;
        processed_reads.emplace_back(unique_ptr<PairedRead>(new PairedRead));
        unique_ptr<PairedRead>& read = processed_reads[i];
        read->id = i;
        
        read->left_fa_id = reads.ids[i];
        read->left_seq = reads.seqs[i];
        read->left_phred = reads.phred[i];
        
        read->right_fa_id = pairs.ids[i];
        read->right_seq = pairs.seqs[i];
        read->right_phred = pairs.phred[i];
    }
    
    // filling up the pi_x_n vectors    
    int ID = 0;
    int count = 0;
    for (auto record : alignments.records) {
        string qName = toCString(record.qName);
        int read_id = params.qNameToID[qName];
        unique_ptr<PairedRead>& read = processed_reads[read_id];

        // only records that are completely mapped
        if (record.rID != record.INVALID_REFID && 
            record.rNextId != record.INVALID_REFID) {
            
            if (!has_alignment[read_id]) {
                params.eff_num_reads++;
                has_alignment[read_id] = true;
            }
            
            // find where this alignment should be grouped
            int insert_i = -1;
            for (int i = 0; i < read->pi_x_n.size(); ++i) {
                auto elem = read->pi_x_n[i];
                if (get<0>(elem) == record.rID) {
                    if (hasFlagFirst(record) && get<2>(elem) == record.pNext) {
                        get<1>(elem) = record.beginPos;
                        insert_i = i;
                    } else if (hasFlagLast(record) && get<1>(elem) == record.pNext) {
                        get<2>(elem) = record.beginPos;
                        insert_i = i;
                    }
                }
                
                if (insert_i == -1) {
                    if (hasFlagFirst(record)) {
                        read->pi_x_n.emplace_back(record.rID, record.beginPos, -1);
                    } else {
                        read->pi_x_n.emplace_back(record.rID, -1, record.beginPos);
                    }
                }
            
                insert_i = -1;
            }
        } 
    }
    
    return;    
}

void isomorph::RsemEstimator::precalc_posteriors(const EMParams& params,
                                                vector<vector<double> >& posteriors) {
    
    cerr << "Precalulcating posterior probabilities." << endl;
    auto& transcripts = params.transcripts;
    auto& reads = params.reads;
    auto& qNameToID = params.qNameToID;
    int num_reads = length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    int read_len = length(reads.seqs[0]);    
    
    if (!params.paired_end) {
        auto& single_reads = params.single_reads;
        
        for (int n = 0; n < num_reads; ++n) {
            vector<double> v;
            const unique_ptr<SingleRead>& read = single_reads[n];
            // isoforms assigned to this read
            for (auto record : read->pi_x_n) {
                int i = get<0>(record);
                int transcript_len = length(transcripts.seqs[i]);
                
                // P(rn|znijk=1)
                double P_sum = 0.;
                int position = get<1>(record);
                double Prf = 1;
                double Prr = 1;
                for (int j = 0; j < read_len && j + position < transcript_len; ++j) {
                    // forward direction
                    if (reads.seqs[n][j] != transcripts.seqs[i][position + j]) {
                        Prf *= 0.5;
                    }
                    
                    // reverse direction
                    if (reverse_complement(reads.seqs[n][j]) != toupper(transcripts.seqs[i][position + j])) {
                        Prr *= 0.5;
                    }
                }
                
                P_sum += Prf * 0.5; // orientation probability
                P_sum += Prr * 0.5;
                v.emplace_back(P_sum);            
            }
            posteriors.emplace_back(v);
        }
    } else {
        auto& pairs = params.pairs;
        auto& paired_reads = params.paired_reads;
        
        for (int n = 0; n < num_reads; ++n) {
            vector<double> v;
            const unique_ptr<PairedRead>& read = paired_reads[n];
            // isoforms assigned to this read
            for (auto record : read->pi_x_n) {
                int i = get<0>(record);
                int transcript_len = length(transcripts.seqs[i]);
                
                // P(rn|znijk=1)
                double P_sum = 0.;
                double Prf = 1;
                double Prr = 1;
                int position = get<1>(record);
                
                // first mate
                for (int j = 0; j < read_len && j + position < transcript_len; ++j) {
                    // forward direction
                    if (reads.seqs[n][j] != transcripts.seqs[i][position + j]) {
                        Prf *= 0.5;
                    }
                    // reverse direction
                    if (reverse_complement(reads.seqs[n][j]) != toupper(transcripts.seqs[i][position + j])) {
                        Prr *= 0.5;
                    }
                }
                
                // second mate
                position = get<2>(record);
                for (int j = 0; j < read_len && j + position < transcript_len; ++j) {
                    // forward direction
                    if (pairs.seqs[n][j] != transcripts.seqs[i][position + j]) {
                        Prf *= 0.5;
                    }
                    // reverse direction
                    if (reverse_complement(pairs.seqs[n][j]) != toupper(transcripts.seqs[i][position + j])) {
                        Prr *= 0.5;
                    }
                }                
                
                P_sum += Prf * 0.5; // orientation probability
                P_sum += Prr * 0.5;
                v.emplace_back(P_sum);            
            }
            posteriors.emplace_back(v);
        }        
    }
    
    cerr << "Done precalculating posteriors." << endl;                                                    
    return;
}

void isomorph::RsemEstimator::EMAlgorithm(EMParams& params, EMResult& result) {
    auto& transcripts = params.transcripts;
    auto& reads = params.reads;
    auto& qNameToID = params.qNameToID;
    int num_reads = length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    int read_len = length(reads.seqs[0]);

    bool strand_specific = false;
    
    // precalculate posterior sums
    vector<vector<double> > read_posteriors;
    precalc_posteriors(params, read_posteriors);

    vector<double> expressions(num_transcripts, 1./(num_transcripts));
    vector<double> pre_m_expressions(num_transcripts, 0);
    int iter = 0;
    
    cerr << "Entering the EM iterations." << endl;
    if (!params.paired_end) {
        auto& single_reads = params.single_reads;
        do {
            // E-step
            pre_m_expressions.assign(num_transcripts, 0);
            for (int n = 0; n < num_reads; ++n) {
                double read_expect_sum = 0;
                const unique_ptr<SingleRead>& read = single_reads[n];
                
                // isoforms joined to this read
                for (int t = 0; t < read->pi_x_n.size(); ++t) {
                    auto record = read->pi_x_n[t];
                    int i = get<0>(record);
                    int transcript_len = seqan::length(transcripts.seqs[i]);
                    double coeff = expressions[i] / transcript_len;
                    // P(rn|znijk=1)
                    double P_sum = read_posteriors[n][t];
                    read_expect_sum += P_sum * coeff;
                }
                
                for (int t = 0; t < read->pi_x_n.size(); ++t) {
                    auto record = read->pi_x_n[t];
                    int i = get<0>(record);
                    int transcript_len = length(transcripts.seqs[i]);
                    pre_m_expressions[i] += read_posteriors[n][t] * expressions[i] / transcript_len / read_expect_sum;
                }
            }
            
            // M-Step
            for (int i = 0; i < num_transcripts; ++i) {
                expressions[i] = pre_m_expressions[i] / params.eff_num_reads;   
            }
        } while (++iter < 1000);
    } else {
        auto& paired_reads = params.paired_reads;
        do {
            // E-step
            pre_m_expressions.assign(num_transcripts, 0);
            for (int n = 0; n < num_reads; ++n) {
                double read_expect_sum = 0;
                const unique_ptr<PairedRead>& read = paired_reads[n];
                
                // isoforms joined to this read
                for (int t = 0; t < read->pi_x_n.size(); ++t) {
                    auto record = read->pi_x_n[t];
                    int i = get<0>(record);
                    int transcript_len = seqan::length(transcripts.seqs[i]);
                    double coeff = expressions[i] / transcript_len;
                    // P(rn|znijk=1)
                    double P_sum = read_posteriors[n][t];
                    read_expect_sum += P_sum * coeff;
                }
                
                for (int t = 0; t < read->pi_x_n.size(); ++t) {
                    auto record = read->pi_x_n[t];
                    int i = get<0>(record);
                    int transcript_len = length(transcripts.seqs[i]);
                    pre_m_expressions[i] += read_posteriors[n][t] * expressions[i] / transcript_len / read_expect_sum;
                }
            }
            
            // M-Step
            for (int i = 0; i < num_transcripts; ++i) {
                expressions[i] = pre_m_expressions[i] / params.eff_num_reads;   
            }
        } while (++iter < 1000);        
    }
    
    for (int i = 0; i < num_transcripts; ++i) {
        result.relative_expressions.emplace_back(expressions[i]);
    }
    
    cerr << "EM is done." << endl;
}

void isomorph::RsemEstimator::output_result(const FastAData& transcripts, 
                                            const EMResult& result, 
                                            const string filename) {
    
    cerr << "Outputing results." << endl;                                                
    ofstream output;
    output.open(filename.c_str(), ofstream::out | ofstream::trunc);
    
    int num_transcripts = result.relative_expressions.size();
    double sum = 0;
    // needed to transform to TPM according to RSEM paper
    double v_sum = 0;
    vector<double> tpm(num_transcripts, 0);
    
    for (int i = 0; i < num_transcripts; ++i) {
        sum += result.relative_expressions[i];
        // calculate TPM
        int t_len = seqan::length(transcripts.seqs[i]);
        tpm[i] = (result.relative_expressions[i] / t_len);
        v_sum += tpm[i];
    }
    
    double tpm_sum = 0;
    for (int i = 0; i < num_transcripts; ++i) {
        tpm[i] *= 1e6 / v_sum;
        tpm_sum += tpm[i];
        
        // first ni_i then tpm_i
        output << ">" << transcripts.ids[i] << endl;
        output << result.relative_expressions[i];
        output << '\t' << tpm[i] << endl;
    }
    
    cerr << "Expression sum: " << sum << endl;
    cerr << "TPM sum: " << tpm_sum << endl;
    output.close();
    cerr << "Done outputing results." << endl;
    return;
}
