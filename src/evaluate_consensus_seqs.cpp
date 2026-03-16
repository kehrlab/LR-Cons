#include "sequence_snippets.hpp"

#include <filesystem>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>

#include <htslib/faidx.h>
#include <seqan/align.h>
#include <seqan/align/dp_profile.h>
#include <seqan/score/score_edit.h>
#include <seqan/sequence.h>


int main(int argc, char **argv)
{
    // take a file containing all variants + location of consensus sequences for a POA method
    // take a truth path

    // parse arguments
    if (argc != 4)
    {
        std::cout << "Usage: ./evaluate_consensus_seqs <CONSENSUS_LIST.tsv> <TRUTH_DIR> <OUT_FILE>" << std::endl;
        return -1;
    }

    std::string consensus_list = argv[1];
    std::string truth_dir = argv[2];
    std::string out_file = argv[3];

    // open consensus list
    std::ifstream f_in(consensus_list);
    if (!f_in.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << consensus_list << " for reading." << std::endl;
        return -1;
    }

    // open output file
    std::ofstream f_out(out_file);
    if (!f_out.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << out_file << " for writing." << std::endl;
        return -1;
    }
    f_out << "VARIANT\tMETHOD\tCONSENSUS_IDX\tCONSENSUS_LENGTH\tALLELE_MATCH\tSCORE\tALIGNED_BASES\tMATCHED_BASES" << std::endl;

    std::string variant_id, method, fa_path;
    while (true)
    {
        int n = 0;
        bool eof = true;

        std::vector<std::string> truth_files;
        std::vector<std::string> consensus_files;
        std::vector<std::string> variant_ids;
        std::vector<std::string> methods;

        while (f_in >> variant_id >> method >> fa_path)
        {
            variant_ids.push_back(variant_id);
            methods.push_back(method);
            consensus_files.push_back(fa_path);
            truth_files.push_back(truth_dir + "/" + variant_id + "_truth.fa");

            if (n == 8000)
            {
                eof = false;
                break;
            }
        }

        std::vector<std::string> output_strings;

        #pragma omp parallel for num_threads(16)
        for (int idx = 0; idx < variant_ids.size(); ++idx)
        {
            if (!std::filesystem::exists(truth_files[idx]))
                throw std::runtime_error("Allele sequence file " + truth_files[idx] + " does not exist.");

            if (!std::filesystem::exists(consensus_files[idx]))
                throw std::runtime_error("Consensus sequence file " + consensus_files[idx] + " does not exist.");

            std::vector<std::string> allele_seqs = load_from_fastq(truth_files[idx]);

            if (allele_seqs.size() != 2)
                throw std::runtime_error("Unexpected number of allele sequences for variant " + variant_ids[idx]);


            // load the consensus sequences
            std::vector<std::string> consensus_seqs = load_from_fastq(consensus_files[idx]);

            // align consensus sequences against allele sequences (semi-global, edit-distance)
            std::vector<int> consensus_assignments(consensus_seqs.size());
            std::vector<int> alignment_scores(consensus_seqs.size());
            std::vector<int> aln_lengths(consensus_seqs.size());
            std::vector<int> match_counts(consensus_seqs.size());
            for (int i = 0; i < consensus_seqs.size(); ++i)
            {
                int max_idx = -1;
                int max_score = -1;
                for (int j = 0; j < allele_seqs.size(); ++j)
                {
                    seqan::StringSet<seqan::DnaString> stringSet;

                    seqan::appendValue(stringSet, consensus_seqs[i].c_str());
                    seqan::appendValue(stringSet, allele_seqs[j].c_str());

                    seqan::Align<seqan::DnaString> align(stringSet);

                    int score = seqan::globalAlignment(
                        align,
                        seqan::Score<int, seqan::Simple>(0, -1, -1),
                        seqan::AlignConfig<false, true, true, false>()
                    );


                    auto &r0 = seqan::row(align, 0);
                    auto &r1 = seqan::row(align, 1);
                    unsigned aligned_bases = 0;
                    unsigned n_matches = 0;
                    for (unsigned k = 0; k < seqan::length(r0); ++k)
                    {
                        if (!seqan::isGap(r0, k) && !seqan::isGap(r1, k))
                        {
                            ++aligned_bases;
                            if (r0[k] == r1[k])
                                ++n_matches;
                        }
                    }

                    if (max_idx < 0 || score > max_score)
                    {
                        max_idx = j;
                        max_score = score;
                        aln_lengths[i] = aligned_bases;
                        match_counts[i] = n_matches;
                    }

                    consensus_assignments[i] = max_idx;
                    alignment_scores[i] = max_score;
                }
            }

            // write one line per consensus sequence
            #pragma omp critical
            for (int i= 0; i < consensus_seqs.size(); ++i)
                output_strings.push_back(
                    variant_ids[idx] + "\t" + methods[idx] + "\t" + std::to_string(i+1) + "\t" + std::to_string(consensus_seqs[i].length()) + "\t" + std::to_string(consensus_assignments[i]) + "\t" + std::to_string(alignment_scores[i]) + "\t" + std::to_string(aln_lengths[i]) + "\t" + std::to_string(match_counts[i])
                );
        }

        for (auto & s : output_strings)
            f_out << s << std::endl;

        if (eof)
            break;
    }

    f_out.close();
    f_in.close();
}
