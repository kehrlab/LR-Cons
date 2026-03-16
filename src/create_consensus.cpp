#include "sequence_snippets.hpp"
#include "consensus_methods.hpp"

#ifndef MAX_COV
#define MAX_COV 256
#endif

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>

#include <seqan/align.h>
#include <seqan/align/dp_profile.h>
#include <seqan/score/score_edit.h>
#include <seqan/sequence.h>
#include <omp.h>

int8_t create_consensus(
    std::vector<std::string> & consensus_seqs,
    std::vector<int8_t> & sequence_assignments,
    std::string & sequence_counts,
    std::vector<std::string> sequences,
    std::vector<std::string> seq_names,
    std::string method,
    std::string phasing_dir,
    std::string phasing_regions,
    aln_params p,
    bool semiglobal
) {
    int8_t status = 0;
    consensus_seqs.clear();
    sequence_assignments.clear();

    if (method == "ITER_POA")
    {
        status = create_consensus_iterative(consensus_seqs, sequence_assignments, sequence_counts, sequences, p, 500, 0.96, semiglobal);
    } else if (method == "ITER_POA_PERMISSIVE")
    {
        status = create_consensus_iterative(consensus_seqs, sequence_assignments, sequence_counts, sequences, p, 500, 0.90, semiglobal);
    } else if (method == "abPOA")
    {
        status = create_consensus_abpoa(consensus_seqs, sequence_assignments, sequence_counts, sequences, 2, p, 500, false, true);
    }
    else if(method == "Molephase" || method == "TRUE_HAP" || method == "KMEANS" || method == "fuzzy")
    {
        status = create_consensus_general(
            consensus_seqs, sequence_assignments, sequence_counts,
            sequences, seq_names,
            phasing_dir, phasing_regions,
            p, 500,
            method,
	    semiglobal
        );
    }
    return status;
}


int main(int argc, char** argv)
{
    if (argc < 12) {
        std::cerr << "Invalid number of positional arguments: <INFILE> <OUTDIR> <OUTFILE> <STATS_FILE> <PHASING_DIR> <MATCH> <MISMATCH> <OPEN> <EXTEND> <MODE> <METHODS...>." << std::endl;
        return -1;
    }

    std::string filepath = argv[1];     // file with <READ_FASTQ_FILE>\t<VARIANT_ID>\t<PHASING_REGIONS> in each line
    std::string outpath = argv[2];      // where to write the resulting consensus sequences to
    std::string outfile = argv[3];
    std::string stats_file = argv[4];
    std::string phasing_dir = argv[5];  // does not have to be valid unless one of the methods is Molephase

    int8_t match = std::stoi(argv[6]);
    int8_t mismatch = std::stoi(argv[7]);
    int8_t gap = std::stoi(argv[8]);
    int8_t ext = std::stoi(argv[9]);
    std::string global = argv[10];
    bool semiglobal = global == "semiglobal";
    if (global != "global" && global != "local" && global != "semiglobal")
    {
	    std::cerr << "[ERROR] Invalid alignment mode: " << global << ". Must be on of [global, local, semiglobal]" << std::endl;
	    return -1;
    }
    aln_params p {(global == "global"), match, mismatch, gap, ext};

    std::unordered_set<std::string> valid_methods {
        "ITER_POA", "ITER_POA_PERMISSIVE", "KMEANS", "abPOA", "Molephase", "TRUE_HAP", "fuzzy"
    };

    std::vector<std::string> methods;
    for (int i = 11; i < argc; ++i)
    {
        if (valid_methods.find(argv[i]) == valid_methods.end())
        {
            std::cerr << "[ERROR] Invalid method: " << argv[i] << std::endl;
            return -1;
        }
        methods.push_back(argv[i]);
    }
    if (methods.size() < 1)
    {
        std::cerr << "[ERROR] At least one valid consensus method (ITER_POA, ITER_POA_PERMISSIVE, abPOA, Molephase, TRUE_HAP, KMEANS, fuzzy) required." << std::endl;
        return -1;
    }

    std::string fastq_file, variant_id, phasing_regions;

    std::ifstream f_in(filepath);
    if (!f_in.is_open())
    {
        std::cerr << "[ERROR] Could not open " << filepath << " for reading." << std::endl;
        return -1;
    }

    std::ofstream f_out(outfile);
    if (!f_out.is_open())
    {
        std::cerr << "[ERROR] Could not open " << outfile << " for writing" << std::endl;
        return -1;
    }

    std::ofstream t_out(stats_file);
    if (!t_out.is_open())
    {
        std::cerr << "[ERROR] Could not open " << stats_file << " for writing" << std::endl;
        return -1;
    }

    if (!std::filesystem::exists(outpath))
        std::filesystem::create_directories(outpath);
    if (!std::filesystem::is_directory(outpath))
    {
        std::cerr << "[ERROR] Directory " << outpath << " does not exist and could not be created" << std::endl;
        return -1;
    }

    while (true)
    {
        // load the next 1000 entries
        int n = 0;
        std::vector<std::vector<std::string>> sequence_lists;
        std::vector<std::vector<std::string>> seq_name_lists;
        std::vector<std::string> variant_ids;
        std::vector<std::string> phasing_region_list;

        bool full_buffer = false;
        bool eof = true;
        while (f_in >> fastq_file >> variant_id >> phasing_regions)
        {
            std::vector<std::string> sequences;
            std::vector<std::string> seq_names;
            bool read_status = load_from_fastq(sequences, seq_names, fastq_file);
            if (sequences.size() < 1 || !read_status) {
                std::cerr << "[INFO] No sequences loaded for variant " << variant_id << std::endl;
                continue;
            }

            std::vector<std::string> sorted_seqs;
            std::vector<std::string> sorted_names;
            std::vector<int> indices (sequences.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&](const int & a, const int & b){return sequences[a].size() > sequences[b].size();});
            for (auto idx : indices)
            {
                sorted_seqs.push_back(sequences[idx]);
                sorted_names.push_back(seq_names[idx]);
            }
            sequences.clear();
            seq_names.clear();


            // remove sequences longer than 10kbp
            int idx = 0;
            while (sorted_seqs[idx].length() > 10000)
                ++idx;
            if (idx > 0) {
                sorted_seqs.erase(sorted_seqs.begin(), sorted_seqs.begin() + idx);
                sorted_names.erase(sorted_names.begin(), sorted_names.begin() + idx);
            }

            if (sorted_seqs.size() < 1) {
                std::cerr << "[INFO] No sequences smaller than 10kbp remaining for variant " << variant_id << std::endl;
                continue;
            }

            if (sorted_seqs.size() >= MAX_COV)
            {
                sorted_seqs.erase(sorted_seqs.begin() + MAX_COV, sorted_seqs.end());
                sorted_names.erase(sorted_names.begin() + MAX_COV, sorted_names.end());
                std::cerr << "[INFO] Too many sequences for variant " << variant_id << ". Reduced to " << MAX_COV << " largest." << std::endl;
            }

            sequence_lists.push_back(sorted_seqs);
            seq_name_lists.push_back(sorted_names);
            variant_ids.push_back(variant_id);
            phasing_region_list.push_back(phasing_regions);

            ++n;
            if (n == 2000)
            {
                full_buffer = true;
                eof = false;
                break;
            }
        }

        std::vector<std::string> fa_list_strings;
        std::vector<std::string> stats_strings;

        #pragma omp parallel for num_threads(14) // 16
        for (int i = 0; i < sequence_lists.size(); ++i)
        {
            std::vector<std::string> & sequences = sequence_lists[i];
            std::vector<std::string> & seq_names = seq_name_lists[i];
            std::string & v_name = variant_ids[i];
            // std::cout << v_name << std::endl;

            for (int j = 0; j < methods.size(); ++j)
            {
                int8_t status = 0;
                std::string fa_path = outpath + "/consensus_" + v_name + "_" + methods[j] + ".fa";
                std::vector<std::string> consensus_seqs;
                std::vector<std::string> consensus_ids;
                std::vector<int8_t> sequence_assignments;
                std::string sequence_counts;
                std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                status = create_consensus(
                    consensus_seqs, sequence_assignments, sequence_counts,
                    sequences, seq_names,
                    methods[j],
                    phasing_dir, phasing_region_list[i],
                    p,
		    semiglobal
                );
                std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

                for (int k = 0; k < consensus_seqs.size(); ++k)
                    consensus_ids.push_back(v_name + "_consensus_" + methods[j] + "_" + std::to_string(k+1));

                write_to_fasta(consensus_seqs, consensus_ids, fa_path);

                #pragma omp critical
                stats_strings.push_back(
                    v_name + "\t" + methods[j] + "\t" + std::to_string(ms) + "\t" + std::to_string((int16_t) status) + "\t" + sequence_counts
                );

                #pragma omp critical
                fa_list_strings.push_back(
                    v_name + "\t" + methods[j] + "\t" + fa_path
                );
            }
        }

        for (auto & s : stats_strings)
            t_out << s << std::endl;
        for (auto & s : fa_list_strings)
            f_out << s << std::endl;

        if (eof)
            break;
    }

    f_in.close();
    f_out.close();
    t_out.close();

    return 0;
}
