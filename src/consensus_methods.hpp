#ifndef MOLEPHASE_HDR
#define MOLEPHASE_HDR

#include "kmer_clustering.hpp"
#include "utils.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "htslib/faidx.h"
#include <spoa/alignment_engine.hpp>
#include <spoa/graph.hpp>
#include <seqan/align.h>
#include <seqan/align/dp_profile.h>
#include <seqan/score/score_edit.h>
#include <seqan/sequence.h>
#include "rapidfuzz/fuzz.hpp"

#include "../include/abpoa.h"


// Function that searches for sequence names in phased Molephase fq files
// and assigns sequences to clusters based on that phasing.
inline bool get_molephase_sequence_assignments(
    std::vector<int8_t> & sequence_assignments,
    std::string phasing_dir,
    std::string phasing_region,
    std::vector<std::string> seq_names
)
{
    sequence_assignments = std::vector<int8_t>(seq_names.size(), -1);

    // determine phasing region files and check their existence
    std::string hap1_file = phasing_dir + "/" + phasing_region + ".hap1.fq";
    std::string hap2_file = phasing_dir + "/" + phasing_region + ".hap2.fq";
    std::string ambig_file = phasing_dir + "/" + phasing_region + ".unphased.fq";

    if (!std::filesystem::exists(hap1_file) || !std::filesystem::exists(hap2_file) || !std::filesystem::exists(ambig_file))
    {
        std::cerr << "[ERROR] Molephase output appears to be missing" << std::endl;
        return false;
    }

    // process sequence names -> everything before second underscore
    for (int i = 0; i < seq_names.size(); ++i) {
	if (seq_names[i].find("ccs") == std::string::npos) // not PacBio: format S1_xxx_[start_end]
            seq_names[i] = seq_names[i].substr(0, seq_names[i].find("_", seq_names[i].find("_") + 1));
        else // PacBio: format S1/xxx/ccs[_start_end]
	    seq_names[i] = seq_names[i].substr(0, seq_names[i].find("_"));
    }

    bool indices_built = true;
    // Make this critical because the fasta-indexing can otherwise cause issues in multithreading
    #pragma omp critical
    {
        // open fastq files
        int status;

        status = fai_build(hap1_file.c_str());
        if (status < 0)
        {
            std::cerr << "Error building index for file " << hap1_file << std::endl;
            indices_built =  false;
        }
        status = fai_build(hap2_file.c_str());
        if (status < 0)
        {
            std::cerr << "Error building index for file " << hap2_file << std::endl;
            indices_built =  false;
        }
        // status = fai_build(ambig_file.c_str());
        // if (status < 0)
        // {
        //     std::cerr << "Error building index for file " << hap1_file << std::endl;
        //     indices_built = false;
        // }

        if (indices_built)
        {
            // open all three fastq files
            faidx_t * faidx_1 = fai_load(hap1_file.c_str());
            faidx_t * faidx_2 = fai_load(hap2_file.c_str());
            // faidx_t * faidx_u = fai_load(ambig_file.c_str());

            // search for read assignment using indexed fasta files
            for (int i = 0; i < seq_names.size(); ++i)
            {
                if (faidx_has_seq(faidx_1, seq_names[i].c_str()))
                    sequence_assignments[i] = 0;
                else if (faidx_has_seq(faidx_2, seq_names[i].c_str()))
                    sequence_assignments[i] = 1;
                // else if (faidx_has_seq(faidx_u, seq_names[i].c_str()))
                //     sequence_assignments[i] = -1;
            }

            // close files
            fai_destroy(faidx_1);
            fai_destroy(faidx_2);
            // fai_destroy(faidx_u);
        }
    }

    return indices_built;
}

// Function that determines haplotype for each sequence based on sequence name (S1_*: Hap1, S2_*: Hap2)
inline bool get_true_sequence_assignments(
    std::vector<int8_t> & sequence_assignments,
    std::vector<std::string> seq_names
)
{
    sequence_assignments = std::vector<int8_t>(seq_names.size(), -1);
    for (uint16_t i = 0; i < seq_names.size(); ++i)
    {
        if (seq_names[i].substr(0, 2) != "S1" && seq_names[i].substr(0,2) != "S2")
        {
            std::cerr << "[ERROR] In TRUE_HAP: unexpected sequence name; sequences cannot be assigned." << std::endl;
            return false;
        }
        sequence_assignments[i] = std::stoi(seq_names[i].substr(1,1)) - 1;
    }
    return true;
}


// Function to cluster a set of sequences based on their k-mer counts.
// Uses k-means clustering.
inline bool get_kmeans_sequence_assignments(
    std::vector<int8_t> & sequence_assignments,
    std::vector<std::string> sequences
)
{
    sequence_assignments = std::vector<int8_t>(sequences.size(), -1);
    std::vector<std::vector<uint8_t>> cnt = count_kmers(sequences);
    if (!reduce_kmers(cnt))
        return false;

    std::vector<int> cl_assignments = k_clust_kmers(cnt, 2);
    for (int i = 0; i < cl_assignments.size(); ++i)
        sequence_assignments[i] = cl_assignments[i];

    return true;
}

// Function that groups sequences based on their Levenshtein distance (calculated with RapidFuzz)
// NOTE: Use of Ratio (instead of PartialRatio) is not ideal because the extracted sequences may differ in length;
// therefore every sequence comparison is performed twice, cropping the longer sequence once on the left and once on the right
inline bool get_fuzzy_sequence_assignments(
    std::vector<int8_t> & sequence_assignments,
    std::vector<std::string> sequences
)
{
    sequence_assignments = std::vector<int8_t>(sequences.size(), -1);

    int n = sequences.size();

    std::vector<double> first_scores(n);
    std::vector<double> last_scores(n);
    first_scores[0] = 1.;

    // calculate similarity of every sequence to the longest sequence in the set
    for (int i = 1; i < n; ++i)
    {
        double s1 = rapidfuzz::fuzz::ratio(
            sequences[0].substr(0, sequences[i].length()),
            sequences[i]
        );
        double s2 = rapidfuzz::fuzz::ratio(
            sequences[0].substr(sequences[0].length() - sequences[i].length()),
            sequences[i]
        );
        first_scores[i] = std::max(s1, s2);
    }

    // choose the median similarity as the threshold for the first cluster
    std::vector<int> indices(sequences.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(
        indices.begin(), indices.end(),
        [&](const int & a, const int & b){return first_scores[a] < first_scores[b];}
    );
    double first_threshold = first_scores[indices[(int) (0.5*indices.size())]];

    // choose the first quartile sequence as the representative of the second cluster
    int second_rep = indices[(int) (0.25*indices.size())];
    last_scores[second_rep] = 1.;

    // calculate similarity of all sequences to the second cluster representative
    for (int i = 0; i < n; ++i)
    {
        if (i != second_rep)
        {
            double s1, s2;
            if (i > second_rep) // sequences earlier in the vector are longer, higher indices contain shorter sequences
            {
                s1 = rapidfuzz::fuzz::ratio(
                    sequences[second_rep].substr(0, sequences[i].length()),
                    sequences[i]
                );
                s2 = rapidfuzz::fuzz::ratio(
                    sequences[second_rep].substr(sequences[second_rep].length() - sequences[i].length()),
                    sequences[i]
                );
            }
            else
            {
                s1 = rapidfuzz::fuzz::ratio(
                    sequences[second_rep],
                    sequences[i].substr(0, sequences[second_rep].length())
                );
                s2 = rapidfuzz::fuzz::ratio(
                    sequences[second_rep],
                    sequences[i].substr(sequences[i].length() - sequences[second_rep].length())
                );
            }
            last_scores[i] = std::max(s1, s2);
        }
    }

    std::iota(indices.begin(), indices.end(), 0);
    std::sort(
        indices.begin(), indices.end(),
        [&](const int & a, const int & b){return last_scores[a] < last_scores[b];}
    );
    double second_threshold = last_scores[indices[(int) (0.5*indices.size())]];

    // assign each sequence to the cluster with the closer score; only assign scores above the median value
    for (int i = 0; i < sequences.size(); ++i)
    {
        double & s1 = first_scores[i];
        double & s2 = last_scores[i];
        if (s1 > s2 && s1 > first_threshold)
            sequence_assignments[i] = 0;
        else if (s2 > s1 && s2 > second_threshold)
            sequence_assignments[i] = 1;
    }

    return true;
}

inline bool merge_sequences(
    std::vector<std::string> & consensus_seqs,
    std::vector<int8_t> & sequence_assignments
)
{
    // align with seqan and calculate similarity
    seqan::StringSet<seqan::DnaString> stringSet;

    if (consensus_seqs[0].length() < consensus_seqs[1].length()) {
        seqan::appendValue(stringSet, consensus_seqs[0].c_str());
        seqan::appendValue(stringSet, consensus_seqs[1].c_str());
    } else {
        seqan::appendValue(stringSet, consensus_seqs[1].c_str());
        seqan::appendValue(stringSet, consensus_seqs[0].c_str());
    }

    // global alignment with free end gaps (first / horizontal / shorter) sequence;
    // score penalizes mismatches and gaps with -1 to calculate dissimilarity
    seqan::Align<seqan::DnaString> align(stringSet);
    int score = seqan::globalAlignment(
        align, seqan::Score<int, seqan::Simple>(0, -1, -1),
        seqan::AlignConfig<false, true, true, false>()
    );

    size_t seq_len = std::min(consensus_seqs[0].length(), consensus_seqs[1].length());
    float similarity = 1 + ((float) score / seq_len); // score is negative mismatch + gap count

    bool merge = similarity > (0.99 - (float) 30 / (float) seq_len);
    if (merge)
    {
        if (consensus_seqs[0].size() < consensus_seqs[1].size())
        {
            consensus_seqs.erase(consensus_seqs.begin());
            for (auto & x : sequence_assignments)
                if (x == 1)
                    x = 0;
        }
        else
        {
            consensus_seqs.erase(consensus_seqs.begin() + 1);
            for (auto & x : sequence_assignments)
                if (x == 1)
                    x = 0;
        }
    }
    return merge;
}

// K-mer based clustering, RapidFuzz, Molephase and True Haplotype methods only differ in
// the way the read clusters are defined. This function determines the
// sequence assignment flexibly based on the chosen method,
// then builds one POA graph per sequence cluster and extracts consensus sequences.
// If two sequences pass the quality filters, another alignment is performed to
// check whether they should be merged due to high similarity.
inline int8_t create_consensus_general(
    std::vector<std::string> & consensus_seqs,
    std::vector<int8_t> & sequence_assignments,
    std::string & sequence_counts,
    std::vector<std::string> sequences,
    std::vector<std::string> seq_names,
    std::string phasing_dir,
    std::string phasing_region,
    aln_params p,
    int margin,
    std::string method,
    bool semiglobal
)
{
    int read_threshold = 2; // number of sequences that must be exceeded for consensus creation

    // divide sequences into clusters with the given method
    bool status = false;
    if (method == "Molephase")
        status = get_molephase_sequence_assignments(sequence_assignments, phasing_dir, phasing_region, seq_names);
    else if (method == "TRUE_HAP")
        status = get_true_sequence_assignments(sequence_assignments, seq_names);
    else if (method == "KMEANS")
        status = get_kmeans_sequence_assignments(sequence_assignments, sequences);
    else if (method == "fuzzy")
        status = get_fuzzy_sequence_assignments(sequence_assignments, sequences);
    else
        std::cerr << "[ERROR] Invalid consensus method: " << method << std::endl;

    if (!status)
        return -1;


    // Once all reads are assigned, create two separate partial order graphs (spoa)
    std::vector<spoa::Graph> poa_graphs(2);
    std::vector<int> read_support {0,0};

    spoa::AlignmentType aln_type = p.global ? spoa::AlignmentType::kNW : spoa::AlignmentType::kSW;
    if (semiglobal)
        aln_type = spoa::AlignmentType::kOV;
    auto alignment_engine = spoa::AlignmentEngine::Create(
        aln_type,
        p.match,
        p.mismatch,
        p.open,
        p.ext
    );

    int assigned_reads = 0;
    int skipped = 0;
    std::vector<bool> sufficient_reads = {false, false}; // keep track of this for debugging / method improvements
    std::vector<bool> sufficient_size = {false, false};

    // Determine the number of reads assigned to each group
    for (int8_t i = 0; i < 2; ++i)
        for (int8_t j = 0; j < sequence_assignments.size(); ++j)
            if (sequence_assignments[j] == i)
                ++read_support[i];
    assigned_reads = read_support[0] + read_support[1];

    // Create POA graph for every sufficiently supported group
    bool any_unsupported = false;
    for (int8_t i = 0; i < 2; ++i)
    {
        sufficient_reads[i] = read_support[i] > read_threshold;
        any_unsupported |= read_support[i] == 0;
        if (sufficient_reads[i])
        {
            for (int j = 0; j < sequence_assignments.size(); ++j)
            {
                if (sequence_assignments[j] == i)
                {
                    auto alignment = alignment_engine->Align(
                        sequences[j],
                        poa_graphs[i]
                    );
                    poa_graphs[i].AddAlignment(alignment, sequences[j]);
                }
            }
        }
    }

    for (int8_t i = 0; i < 2; ++i)
    {
        std::string consensus = "";
        if (sufficient_reads[i])
            consensus = poa_graphs[i].GenerateConsensus(read_threshold);
        sufficient_size[i] = consensus.size() > margin;
        if (sufficient_size[i])
        {
           consensus_seqs.push_back(consensus);
            if (skipped > 0) // sequence-to-consensus assignment consistent even if consensus sequences are discarded
                for (auto & x : sequence_assignments)
                    if (x == i)
                        x -= skipped;
        }
        else
        {
            for (auto & x : sequence_assignments)
                if (x == i)
                    x = -1;
            ++skipped;
        }
    }

    int8_t consensus_status = 0;
    for (int i = 0; i < 2; ++i)
    {
        if (!sufficient_reads[i])
            consensus_status += (1 << i);
        else if (!sufficient_size[i])
            consensus_status += (1 << (i+2));
    }

    sequence_counts = "";
    for (int i = 0; i < 2; ++i)
        sequence_counts += (std::to_string(read_support[i]) + ",");
    sequence_counts += std::to_string(sequences.size());

    // check whether sequences need to be merged
    if (consensus_seqs.size() == 2)
        if (merge_sequences(consensus_seqs, sequence_assignments))
            consensus_status += 16;

    if (assigned_reads < 0.6 * sequences.size()) // Flag weak read assignment
        consensus_status += 32;

    return consensus_status;
}



inline int8_t create_consensus_iterative(
    std::vector<std::string> & consensus_seqs,
    std::vector<int8_t> & sequence_assignments,
    std::string & sequence_counts,
    std::vector<std::string> sequences,
    aln_params p,
    int margin,
    float score_threshold,
    bool semiglobal
)
{
    sequence_assignments = std::vector<int8_t>(sequences.size(), -1);

    int read_threshold = 2; // number of sequences that must be exceeded for consensus creation
    int max_graphs = 8;

    // Once all reads are assigned, create two separate partial order graphs (spoa)
    std::vector<spoa::Graph> poa_graphs(max_graphs);
    std::vector<int> read_support(max_graphs, 0);

    spoa::AlignmentType aln_type = p.global ? spoa::AlignmentType::kNW : spoa::AlignmentType::kSW;
    if (semiglobal)
        aln_type = spoa::AlignmentType::kOV;
    auto alignment_engine = spoa::AlignmentEngine::Create(
        aln_type,
        p.match,
        p.mismatch,
        p.open,
        p.ext
    );


    int assigned_reads = 0;
    int skipped = 0;
    int last_resort_alns = 0;
    std::vector<bool> sufficient_reads(max_graphs, false); // keep track of this for debugging / method improvements
    std::vector<bool> sufficient_size(max_graphs, false);


    // Traverse sequences and find the best alignment for each
    // initialize first graph from first sequence
    if (sequences.size() > 0)
    {
        auto alignment = alignment_engine->Align(sequences[0], poa_graphs[0]);
        poa_graphs[0].AddAlignment(alignment, sequences[0]);
        sequence_assignments[0] = 0;
        read_support[0] = 1;
    }

    int n_graphs = 1;
    for (int8_t i = 1; i < sequences.size(); ++i)
    {
        bool added = false;
        int32_t max_score = 0;
        int8_t max_idx = -1;
        spoa::Alignment best_aln;

        // try to find a matching graph among the existing
        for (int j = 0; j < n_graphs; ++j)
        {
            int32_t score;
            auto alignment = alignment_engine->Align(
                sequences[i],
                poa_graphs[j],
                &score
            );
            if ((max_idx < 0 || score > max_score) && alignment.size() > 0)
            {
                max_idx = j;
                max_score = score;
                best_aln = alignment;
            }
        }

        // found a possible match
        if (max_idx >= 0)
        {
            float aln_score = (float) max_score / sequences[i].length();
            int aln_bases = 0;
            for (auto it : best_aln)
                if (it.second >= 0 && it.second < sequences[i].length())
                    ++aln_bases;

            if (aln_score > score_threshold)
            {
                poa_graphs[max_idx].AddAlignment(best_aln, sequences[i]);
                ++read_support[max_idx];
                added = true;
            }
            // else if (n_graphs == max_graphs)
            // {
            //     ++last_resort_alns;
            //     poa_graphs[max_idx].AddAlignment(best_aln, sequences[i]);
            //     ++read_support[max_idx];
            //     added = true;
            // }

            if (added)
                sequence_assignments[i] = max_idx;
        }

        if (!added && n_graphs < max_graphs)
        {
            auto alignment = alignment_engine->Align(sequences[i], poa_graphs[n_graphs]);
            poa_graphs[n_graphs].AddAlignment(alignment, sequences[i]);
            sequence_assignments[i] = n_graphs;
            ++read_support[n_graphs];
            ++n_graphs;
        }
    }

    // extract consensus sequences for the two strongest graphs
    std::vector<int> indices(n_graphs);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(
        indices.begin(), indices.end(),
        [&](const int & a, const int & b) {return read_support[a] > read_support[b];}
    );

    sequence_counts = "";
    for (auto & idx : indices)
        sequence_counts += (std::to_string(read_support[idx]) + ",");
    sequence_counts += std::to_string(sequences.size());

    // traverse consensus sequences in order of read support and generate fina sequences (with min. read support)
    int new_idx = 0;
    int8_t consensus_status = 0;
    std::vector<int8_t> read_assignments(sequences.size(), -1);
    for (auto & idx : indices)
    {
        sufficient_reads[idx] = read_support[idx] > read_threshold;
        if (!sufficient_reads[idx])
        {
            consensus_status += (1 << new_idx);
            if (new_idx == 0)
                consensus_status += (1 << (new_idx+1));
            break;
        }

        std::string consensus = "";
        consensus = poa_graphs[idx].GenerateConsensus(read_threshold); // require at least two reads

        sufficient_size[idx] = consensus.length() > margin;
        if (sufficient_size[idx])
        {
            consensus_seqs.push_back(consensus);
            for (int j = 0; j < sequence_assignments.size(); ++j)
                if (sequence_assignments[j] == idx)
                    read_assignments[j] = new_idx;
            ++new_idx;
            assigned_reads += read_support[idx];
        } else {
            consensus_status += (1 << (new_idx + 2));
        }

        // create only two sequences max.
        if (new_idx == 2)
            break;
    }

    sequence_assignments = read_assignments;

    if (consensus_seqs.size() == 2)
        if (merge_sequences(consensus_seqs, sequence_assignments))
            consensus_status += 16;

    if (assigned_reads <= (int) (0.6 * sequences.size()) || last_resort_alns > 0.3 * sequences.size())
        consensus_status += 32;

    return consensus_status;
}


// Function to create up to two consensus sequences with abPOA
inline int8_t create_consensus_abpoa(
    std::vector<std::string> & consensus,
    std::vector<int8_t> & consensus_id,
    std::string & sequence_counts,
    const std::vector<std::string> & sequences,
    uint8_t n_cons,
    aln_params p,
    int margin,
    bool merge,
    bool progressive_poa
)
{
    int read_threshold = 2;
    std::unordered_map<char, uint8_t> base2int {
        {'A', 0},
        {'a', 0},
        {'C', 1},
        {'c', 1},
        {'G', 2},
        {'g', 2},
        {'T', 3},
        {'t', 3},
        {'N', 4},
        {'n', 4}
    };

    std::vector<char> int2base {'A', 'C', 'G', 'T', 'N'};

    consensus_id = std::vector<int8_t> (sequences.size(), -1);
    consensus = std::vector<std::string>();

    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    abpt->align_mode = p.global ? 0 : 1;

    abpt->progressive_poa = progressive_poa ? 1 : 0;
    abpt->gap_mode = ABPOA_AFFINE_GAP;
    abpt->gap_open1 = std::abs(p.open);
    abpt->gap_ext1 = std::abs(p.ext);
    abpt->mismatch = std::abs(p.mismatch);
    abpt->match = p.match;


    abpt->out_msa = 1;
    abpt->out_cons = 1;
    abpt->w = 6;
    abpt->k = 9;
    abpt->min_w = 10;

    abpt->max_n_cons = n_cons;

    abpoa_post_set_para(abpt);

    // store sequence lengths in array and transform vector of sequences to array of arrays of integers
    int *seq_lens = new int [sequences.size()];
    uint8_t ** seqs = new uint8_t* [sequences.size()];
    int **weights = new int* [sequences.size()];
    for (size_t i = 0; i < sequences.size(); ++i)
    {
        seq_lens[i] = sequences[i].length();
        seqs[i] = new uint8_t [seq_lens[i]];
        weights[i] = new int [seq_lens[i]];
        for (size_t j = 0; j < seq_lens[i]; ++j)
        {
            weights[i][j] = 1;
            seqs[i][j] = base2int[sequences[i][j]];
        }
    }

    // POA
    abpt->use_qv = 1;
    abpoa_msa(ab, abpt, sequences.size(), NULL, seq_lens, seqs, weights, NULL);
    abpoa_cons_t *abc = ab->abc;

    int8_t consensus_status = 0;
    int assigned_reads = 0;
    std::vector<bool> sufficient_reads(2, false);
    std::vector<bool> sufficient_size(2, false);
    sequence_counts = "";
    for (uint8_t i = 0; i < abc->n_cons; ++i)
    {
        // only include sequences with sufficient read support
        assigned_reads += abc->clu_n_seq[i];
        sequence_counts += (std::to_string(abc->clu_n_seq[i]) + ",");
        sufficient_reads[i] = abc->clu_n_seq[i] > read_threshold;
        if (!sufficient_reads[i])
            continue;

        std::string temp_seq = "";
        for (uint16_t j = 0; j < abc->cons_len[i]; ++j)
            if (abc->cons_cov[i][j] >= read_threshold) // the same method that is used in SPOA, extracting only bases covered by two reads
                temp_seq += int2base[abc->cons_base[i][j]];

        // do not include too short sequences
        sufficient_size[i] = temp_seq.length() >= margin;
        if (!sufficient_size[i])
            continue;

        for (uint16_t j = 0; j < abc->clu_n_seq[i]; ++j)
            consensus_id[abc->clu_read_ids[i][j]] = consensus.size();
        consensus.push_back(temp_seq);
    }
    sequence_counts += std::to_string(sequences.size());

    // free sequence variables
    for (size_t i = 0; i < sequences.size(); ++i)
    {
        delete[] seqs[i];
        delete[] weights[i];
    }
    delete[] seqs;
    delete[] weights;
    delete[] seq_lens;
    // free abpoa variables
    abpoa_free(ab);
    abpoa_free_para(abpt);

    // check whether sequences need to be merged
    if (consensus.size() == 2) // && merge?
        if (merge_sequences(consensus, consensus_id))
            consensus_status += 16;

    for (int i = 0; i < 2; ++i)
    {
        if (!sufficient_reads[i])
            consensus_status += (1 << i);
        else if (!sufficient_size[i])
            consensus_status += (1 << (i+2));
    }

    if (assigned_reads < 0.6 * sequences.size())
        consensus_status += 32;

    return consensus_status;
}



#endif
