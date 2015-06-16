#include <ctype.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>

#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>

#include "em_estimator.h"
#include "utility.h"
#include "count_estimator.h"
#include "read.h"
#include "single_read.h"
#include "paired_read.h"

using namespace seqan;
using namespace std;

void isomorph::EMEstimator::estimate_abundances(CharString reads,
                                                  CharString transcripts,
                                                  CharString pairs) {
    
    string dir = "isomorph-tmp";
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

void isomorph::EMEstimator::preprocess_data(const CharString& transcripts, 
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
    
    run_alignment(reads_str, 
                  pairs_str, 
                  transcripts_str, 
                  output_dir, 
                  paired_end);

    // reads sam data
    isomorph::SamData alignments;
    CharString sam(output_dir + "/isomorph.sam");
    reader.read_sam(sam, &alignments);                                                      
    reader.read_fastq(reads, &params.reads);
    reader.read_fasta(transcripts, &params.transcripts);

    int num_reads = length(params.reads.ids);
    int num_transcripts = length(params.transcripts.ids);

    if (paired_end) {
        reader.read_fastq(pairs, &params.pairs);
        create_paired_end(alignments, params);

        pair<double, double> insert_params;
        estimate_insert_size(alignments, insert_params);
        params.insert_mean = insert_params.first;
        params.insert_stdev = insert_params.second;

    } else {
        create_single_end(alignments, params);
    }

    cerr << "Effective number of reads: " << params.eff_num_reads << endl;
    cerr << "Done preprocessing data." << endl;
    return;
}                   

void isomorph::EMEstimator::create_single_end(const SamData& alignments,
                                                EMParams& params) {
  
    cerr << "Creating the reads." << endl;    
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
    
    cerr << "Done creating the reads." << endl;
    return;  
}

void isomorph::EMEstimator::create_paired_end(const SamData& alignments,
                                                EMParams& params) {

    cerr << "Creating the reads." << endl;
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
            record.rNextId != record.INVALID_REFID &&
            (hasFlagRC(record) && !hasFlagNextRC(record) || 
                !hasFlagRC(record) && hasFlagNextRC(record))) {
                
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
                        get<4>(elem) = record.seq;
                        insert_i = i;
                    } else if (hasFlagLast(record) && get<1>(elem) == record.pNext) {
                        get<2>(elem) = record.beginPos;
                        get<5>(elem) = record.seq;
                        insert_i = i;
                    }
                }
            }
            
            if (insert_i == -1) {
                if (hasFlagFirst(record)) {
                    read->pi_x_n.emplace_back(record.rID, 
                                              record.beginPos, 
                                              -1, 
                                              record.tLen,
                                              record.seq,
                                              "");
                } else {
                    read->pi_x_n.emplace_back(record.rID, 
                                              -1, 
                                              record.beginPos, 
                                              record.tLen,
                                              "",
                                              record.seq);
                }
            }
            insert_i = -1;
        } 
    }
    
    cerr << "Done creating the reads." << endl;
    return;    
}

void isomorph::EMEstimator::precalc_posteriors(const EMParams& params,
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
                        if (j < transcript_len / 3.) {
                            Prf *= 0.3;
                        } else if (j < 2 * transcript_len / 3.) {
                            Prf *= 0.6;
                        } else {
                           Prf *= 0.9;
                        }
                    }
                    
                    // reverse direction
                    if (reverse_complement(reads.seqs[n][j]) != toupper(transcripts.seqs[i][position + j])) {
                        if (j < transcript_len / 3.) {
                            Prf *= 0.3;
                        } else if (j < 2 * transcript_len / 3.) {
                            Prf *= 0.6;
                        } else {
                           Prf *= 0.9;
                        }
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

        // precalculate normalization factors for each isoform for insert sizes
        vector<double> norm_factors(num_transcripts, 0);
        for (int i = 0; i < num_transcripts; ++i) {
            double sum = 0;
            for (int l = 0; l < length(transcripts.seqs[i]); ++l) {
                sum += prob_normal(params.insert_mean, params.insert_stdev, l);
            }
            norm_factors[i] = sum;
        }
        
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
                
                // first mate
                int position = get<1>(record);
                CharString& seq = get<4>(record);
                for (int j = 0; j < read_len && j + position < transcript_len; ++j) {
                    // forward direction
                    if (seq[j] != transcripts.seqs[i][position + j]) {
                        Prf *= 0.5;
                    }
                    
                    // reverse direction
                    if (reverse_complement(seq[j]) != 
                            toupper(transcripts.seqs[i][position + j])) {
                        Prr *= 0.5;
                    }
                }
                
                // second mate
                position = get<2>(record);
                seq = get<5>(record);
                for (int j = 0; j < read_len && j + position < transcript_len; ++j) {
                    if (seq[j] != toupper(transcripts.seqs[i][position + j])) {
                        Prf *= 0.5;
                    }
                    
                    if (reverse_complement(seq[j]) != 
                            toupper(transcripts.seqs[i][position + j])) {
                        Prr *= 0.5;
                    }
                }                
                
                double insert_prob = prob_normal(params.insert_mean, params.insert_stdev, get<3>(record));
                Prf *= insert_prob / norm_factors[i]; // insert size normalized
                Prr *= insert_prob / norm_factors[i];

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

void isomorph::EMEstimator::EMAlgorithm(EMParams& params, EMResult& result) {
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
    
    cerr << "Entering EM iterations." << endl;
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
        } while (++iter < 2000);
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
                    // P(rn | Znijk)
                    auto record = read->pi_x_n[t];
                    int i = get<0>(record);
                    int transcript_len = seqan::length(transcripts.seqs[i]);
                    double coeff = expressions[i] / transcript_len;
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
        } while (++iter < 2000);   
    }
    
    for (int i = 0; i < num_transcripts; ++i) {
        result.relative_expressions.emplace_back(expressions[i]);
    }
    cerr << "EM is done." << endl;
}

void isomorph::EMEstimator::output_result(const FastAData& transcripts, 
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
