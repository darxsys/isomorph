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

using namespace seqan;
using namespace std;

void isomorph::RsemEstimator::estimate_abundances(CharString reads,
                                                  CharString transcripts,
                                                  CharString pairs) {
    bool paired_end = pairs == CharString("") ? false : true;
    Reader reader;
    string transcripts_str(toCString(transcripts)); 
    string reads_str(toCString(reads));
    string pairs_str(toCString(pairs));

    string dir = "bowtie-tmp";
    string mkdir = "mkdir -p " + dir;
    execute_command(mkdir.c_str());

    // builds bowtie index
    string command = "bowtie2-build " + transcripts_str + " " + dir + "/isomorph-index";
    execute_command(command.c_str());

    // runs the alignment
    if (paired_end) {
        command = "bowtie2 -x " + dir + "/isomorph-index -1 " + reads_str + " -2 " +
                toCString(pairs_str) + " -S " + dir + "/isomorph.sam";
    } else {
        command = "bowtie2 -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 \
                   --score-min L,0,-0.1 -p 1 -k 200 -x " + dir + "/isomorph-index -U " + reads_str +
                   " -S " + dir + "/isomorph.sam";
    }
    execute_command(command.c_str());

    // reads sam data
    isomorph::SamData sam_data;
    CharString sam(dir + "/isomorph.sam");
    reader.read_sam(sam, &sam_data);
    EMParams params;

    // print_sam_alignment_records(sam_data.records);
    if (paired_end) {
        preprocess_data(sam_data, 
                        transcripts, 
                        reads, 
                        pairs, 
                        params);
    } else {
        preprocess_data(sam_data, 
                        transcripts, 
                        reads, 
                        params);
    }
    
    EMResult result;
    EMAlgorithm(params, result);
    output_result(params.transcripts, result, "isomorph.abundances.fasta");
    
    // cleaning up
    string rm = "rm -rf " + dir;
    execute_command(rm.c_str());
    return;
}

/*
    Single-end version of data preprocessing.

*/
void isomorph::RsemEstimator::preprocess_data(const SamData& alignments,
                                              const CharString& transcripts, 
                                              const CharString& reads, 
                                              EMParams& params) {
    
    cerr << "Preprocessing data" << endl;                                                  
    Reader reader;
    reader.read_fastq(reads, &params.reads);
    reader.read_fasta(transcripts, &params.transcripts);

    int num_reads = length(params.reads.ids);
    int num_transcripts = length(params.transcripts.ids);

    // this vector will also be used to indicate if a read has any good alignments.
    vector<int> reads_info;
    params.pi_x_n.insert(params.pi_x_n.begin(), 
                         num_reads,
                         reads_info);
    
    vector<bool> read_processed(num_reads, false);
        
    for (int i = 0; i < num_reads; ++i) {
        string qName = toCString(params.reads.ids[i]);
        int pos = qName.find(' ', 0);
        if (pos != string::npos) {
            qName = qName.substr(0, pos);
        }

        params.qNameToID[qName] = i;
    }

    // arange the alignments in the "neighboring matrix"
    int ID = 0;
    int count = 0;
    for (auto record : alignments.records) {
        string qName = toCString(record.qName);
        int read_id = params.qNameToID[qName];

        if (record.rID != record.INVALID_REFID) {
            count++;
            if (read_processed[read_id] == false) {
                params.eff_num_reads++;
            }
            read_processed[read_id] = true;
            params.pi_x_n[read_id].emplace_back(record.rID);
        } 
    }
    
    cerr << "Number of reads: " << num_reads << endl;
    cerr << "Number of alignments: " << count << endl;
    cerr << "Effective number of reads: " << params.eff_num_reads << endl;
    cerr << "Done preprocessing data." << endl;
    return;
}

/*
    Paired-end version of data preprocessing.
*/
void isomorph::RsemEstimator::preprocess_data(const SamData& alignments,
                                              const CharString& transcripts, 
                                              const CharString& reads, 
                                              const CharString& pairs, 
                                              EMParams& params) {

}

void isomorph::RsemEstimator::EMAlgorithm(EMParams& params, EMResult& result) {
    auto& transcripts = params.transcripts;
    auto& reads = params.reads;
    auto& pi_x_n = params.pi_x_n;
    auto& qNameToID = params.qNameToID;
    int num_reads = length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    int read_len = length(reads.seqs[0]);

    bool strand_specific = false;
    
    vector<vector<double> > read_posteriors;
    precalc_posteriors(params, read_posteriors);
        
    vector<double> expressions(num_transcripts, 1./(num_transcripts));
    vector<double> pre_m_expressions(num_transcripts, 0);
    int iter = 0;
    
    cerr << "Entering the EM iterations." << endl;
    do {
        // E-step
        pre_m_expressions.assign(num_transcripts, 0);
        for (int n = 0; n < num_reads; ++n) {
            double read_expect_sum = 0;
            
            // isoforms joined to this read
            for (int t = 0; t < pi_x_n[n].size(); ++t) {
                int i = pi_x_n[n][t];
                int transcript_len = length(transcripts.seqs[i]);
                double coeff = expressions[i] / transcript_len;
                // P(rn|znijk=1)
                double P_sum = read_posteriors[n][t];
                read_expect_sum += P_sum * coeff;
            }
            
            for (int t = 0; t < pi_x_n[n].size(); ++t) {
                int i = pi_x_n[n][t];
                int transcript_len = length(transcripts.seqs[i]);
                pre_m_expressions[i] += read_posteriors[n][t] * expressions[i] / transcript_len / read_expect_sum;
            }
        }
        
        // M-Step
        for (int i = 0; i < num_transcripts; ++i) {
            expressions[i] = pre_m_expressions[i] / params.eff_num_reads;   
        }
        
    } while (++iter < 1000);
    
    for (int i = 0; i < num_transcripts; ++i) {
        result.relative_expressions.emplace_back(expressions[i]);
    }
}

void isomorph::RsemEstimator::precalc_posteriors(const EMParams& params,
                                                vector<vector<double> >& posteriors) {
    
    cout << "Precalulcating posterior probabilities" << endl;
    auto& transcripts = params.transcripts;
    auto& reads = params.reads;
    auto& pi_x_n = params.pi_x_n;
    auto& qNameToID = params.qNameToID;
    int num_reads = length(reads.ids);
    int num_transcripts = length(transcripts.ids);
    int read_len = length(reads.seqs[0]);    
    
    for (int n = 0; n < num_reads; ++n) {
        vector<double> v;
        // isoforms joined to this read
        for (int i : pi_x_n[n]) {
            int transcript_len = length(transcripts.seqs[i]);
            // P(rn|znijk=1)
            double P_sum = 0.;
            for (int j = 0; j < transcript_len - read_len; ++j) {
                double Prn = 1.;
                // forward direction
                for (int l = 0; l < read_len; ++l) {
                    if (reads.seqs[n][l] != transcripts.seqs[i][l+j]) {
                        Prn *= 0.5;
                    }
                }
                P_sum += Prn * 0.5; // orientation probability
                // reverse direction
                Prn = 1;
                for (int l = 0; l < read_len; ++l) {
                    if (reverse_complement(reads.seqs[n][l]) != toupper(transcripts.seqs[i][l+j])) {
                        Prn *= 0.5;
                    }
                }
                P_sum += Prn * 0.5;
                
            }
            v.emplace_back(P_sum);            
        }
        
        if (pi_x_n[n].size() == 0) {
            
        }
        
        posteriors.emplace_back(v);
    }
    
    cout << "Done precalculating posteriors." << endl;                                                    
}                   

void isomorph::RsemEstimator::output_result(const FastAData& transcripts, 
                                            const EMResult& result, 
                                            const string filename) {
                                                
    ofstream output;
    output.open(filename.c_str(), ofstream::out | ofstream::trunc);
    
    double sum = 0;
    for (int i = 0; i < result.relative_expressions.size(); ++i) {
        output << ">" << transcripts.ids[i] << endl;
        output << result.relative_expressions[i] << endl;
        sum += result.relative_expressions[i];
    }
    
    cerr << "Expression sum: " << sum << endl;
    
    output.close();
    return;
}