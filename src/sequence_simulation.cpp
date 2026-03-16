#include "sequence_snippets.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>


std::vector<int> write_seqs(std::string s1, std::string s2, std::string dirname, std::string type, int k, int n_temps, bool random_af)
{
    std::string id_base = type + "_" + std::to_string(k) + "_";
    std::vector<std::string> seqs;
    std::vector<std::string> ids;
    // std::vector<std::string> ids {type + "_" + std::to_string(k) + "_S1", type + "_" + std::to_string(k) + "_S2"};
    std::vector<int> n_reads(2, 0);
    if (random_af)
    {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<> distrib(0., 1.);

        for (int i = 0; i < n_temps; ++i) {
            if (distrib(gen) < 0.5) {
                ++n_reads[0];
                seqs.push_back(s1);
                ids.push_back(id_base + "_S1_" + std::to_string(i));
            } else {
                ++n_reads[1];
                seqs.push_back(s2);
                ids.push_back(id_base + "_S2_" + std::to_string(i));
            }
        }
    }
    else
    {
        for (int i = 0; i < n_temps / 2; ++i) {
            ++n_reads[0];
            seqs.push_back(s1);
            ids.push_back(id_base + "_S1_" + std::to_string(i));
        }
        for (int i = n_temps / 2; i < n_temps; ++i) {
            ++n_reads[1];
            seqs.push_back(s2);
            ids.push_back(id_base + "_S2_" + std::to_string(i));
        }
    }

    write_to_fasta(seqs, ids, dirname + type + "_" + std::to_string(k) + ".fasta");
    return n_reads;
}

bool write_truth(std::vector<std::string> seqs, std::string filename)
{
    std::vector<std::string> seq_names;
    for (int i = 0; i < seqs.size(); ++i)
        seq_names.push_back("HAP_" + std::to_string(i + 1));

    write_to_fasta(seqs, seq_names, filename);

    return true;
}

int main(int argc, char **argv)
{
    // for different sequence lengths, variant sizes and snp rates,
    // create a homozygous reference, a heterozygous deletion,
    // a heterozygous inversion and a heterozygous duplication.
    // home_ref: 3000, heterozygous 1000 each
    std::vector<int> seq_lens = {1500};
    std::vector<int> seq_pos = {500};

    std::vector<int> variant_sizes = {200};
    std::vector<float> snp_rates = {0.005};

    int n_temps = 20;

    std::string fasta_file = "data/simulated_sequences/fasta_list.tsv";
    std::ofstream f_out(fasta_file);
    if (!f_out.is_open())
    {
        std::cerr << "[ERROR] Could not open " << fasta_file << " for writing" << std::endl;
	return -1;
    }

    for (int i = 0; i < seq_lens.size(); ++i)
    {
        int l = seq_lens[i];
        int p = seq_pos[i];
        for (int d : variant_sizes)
        {
            for (float r : snp_rates)
            {
                std::string param_string = std::to_string(seq_lens[i]) +
                    "_" + std::to_string(d) +
                    "_" + std::to_string(r).substr(0,5);
                std::string dirname = "data/simulated_sequences/" + param_string + "/";
                if (!std::filesystem::exists(dirname))
                    std::filesystem::create_directory(dirname);

                std::cout << "Create sequences for l=" << l << ", d=" << d << ", r=" << r << std::endl;

                for (int k = 0; k < 500; ++k) // 1000
                {
                    std::string s = create_random_sequence(5000);
                    std::vector<int> n_reads;

                    // Insert two inversions on different alleles and at different positions
          		    std::string s1 = create_inversion(s, p, d).substr(0,l);
          		    std::string s2 = create_inversion(s, (int) p + d + 100, d * 0.9).substr(0,l);
          		    // std::cout << s1 << std::endl << s2 << std::endl << std::endl;

                    std::vector<std::string> comp_het_seqs = {s1, s2};

                    mutate_seqs(comp_het_seqs, r);
                    n_reads = write_seqs(
                        comp_het_seqs[0], comp_het_seqs[1],
                        dirname, "COMPHET", k, n_temps, false
                    );
                    write_truth(
                        comp_het_seqs,
                        "data/simulated_sequences/allele_sequences/COMPHET_" + std::to_string(l) + "_" + std::to_string(d) + "_" + std::to_string(r).substr(0,5) + "_" + std::to_string(k) + "_truth.fa"
                    );
                    f_out << "data/simulated_sequences/reads/" << param_string << ("/COMPHET_" + std::to_string(k) + ".fq") << "\t" << "COMPHET_" << l << "_" << d << "_" << r << "_" << k << "\t-" << std::endl;
                }
            }
        }
    }
    f_out.close();
}
