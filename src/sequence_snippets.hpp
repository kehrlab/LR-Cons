#include <cassert>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <filesystem>
#include <iostream>

#include "htslib/faidx.h"

static std::vector<char> comp {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
};


inline std::string create_rev_comp (std::string s) {
    std::string s1 = s;
    for (int j = s.size() - 1; j >= 0; --j)
        s1[s.size() - 1 - j] = comp[s[j]];
    return s1;
}


inline std::string create_random_sequence(int l) {
    std::vector<char> nt = {'A', 'C', 'T', 'G'};

    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::vector<float> weights {0.3, 0.2, 0.3, 0.2};
    std::discrete_distribution<> distrib(weights.begin(), weights.end());

    std::string s;
    for (int i = 0; i < l; ++i)
        s += nt[distrib(gen)];

    return s;
}


// create function that introduces random small mutations with given probabilites (very simple), no indels yet
inline void mutate_seqs(std::vector<std::string> & seqs, float p)
{
    std::vector<char> nt = {'A', 'C', 'T', 'G'};

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> distrib(0., 1.);
    std::uniform_int_distribution<> distrib1(0, 3);

    for (auto & s : seqs)
    {
        for (int i = 0; i < s.size(); ++i)
        {
            float r = distrib(gen);

            if (r < p) // mutate the base
            {
                char x = nt[distrib1(gen)];
                while (x == s[i])
                    x = nt[distrib1(gen)];
                s[i] = x;
            }
        }
    }
}


inline std::string create_deletion(std::string s, int pos, int size)
{
    return (s.substr(0, pos) + s.substr(pos + size));
}


inline std::string create_inversion(std::string s, int pos, int size)
{
    return (s.substr(0, pos) + create_rev_comp(s.substr(pos, size)) + s.substr(pos + size));
}


inline std::string create_tandem_duplication(std::string s, int pos, int size)
{
    return (s.substr(0, pos) + s.substr(pos, size) + s.substr(pos, size) + s.substr(pos + size));
}


inline void write_to_fasta(const std::vector<std::string> & sequences, const std::vector<std::string> & ids, std::string filename)
{
    assert(sequences.size() == ids.size());

    std::ofstream f_out(filename);
    if (f_out.is_open())
    {
        for (size_t i = 0; i < sequences.size(); ++i)
        {
            if (sequences[i].length() == 0 || sequences[i] == "")
                continue;

            f_out << ">" << ids[i] << std::endl;

            int64_t j = 0;
            int n_bases = 70;
            while (j < (int64_t) sequences[i].size() - n_bases)
            {
                f_out << sequences[i].substr(j, n_bases) << std::endl;
                j += n_bases;
            }
            f_out << sequences[i].substr(j, sequences[i].size() - j) << std::endl;
        }
    }
}

inline std::vector<std::string> load_from_fastq(std::string filepath)
{
    if (!std::filesystem::exists(filepath))
    {
        std::cerr << "File " << filepath << " does not exists." << std::endl;
        return {};
    }

    int status = fai_build(filepath.c_str());
    if (status < 0)
    {
        std::cerr << "Error building index for file " << filepath << std::endl;
        return {};
    }

    std::vector<std::string> sequences;
    #pragma omp critical
    {
        faidx_t * faidx = fai_load(filepath.c_str());

        int n_seqs = faidx_nseq(faidx);
        for (int i = 0; i < n_seqs; ++i)
        {
            int l = 0;
            std::string s = fai_fetch(faidx, faidx_iseq(faidx, i), &l);
            if (l > 0)
                sequences.push_back(s);
        }
        fai_destroy(faidx);
    }
    return sequences;
}


inline bool load_from_fastq(std::vector<std::string> & sequences, std::vector<std::string> & seq_names, std::string filepath)
{
    sequences.clear();
    seq_names.clear();

    if (!std::filesystem::exists(filepath))
    {
        std::cerr << "File " << filepath << " does not exists." << std::endl;
        return false;
    }

    int status = fai_build(filepath.c_str());
    if (status < 0)
    {
        std::cerr << "Error building index for file " << filepath << std::endl;
        return false;
    }

    #pragma omp critical
    {
        faidx_t * faidx = fai_load(filepath.c_str());

        int n_seqs = faidx_nseq(faidx);
        for (int i = 0; i < n_seqs; ++i)
        {
            int l = 0;
            std::string n = faidx_iseq(faidx, i);
            std::string s = fai_fetch(faidx, faidx_iseq(faidx, i), &l);
            if (l > 0) {
                sequences.push_back(s);
                seq_names.push_back(n);
            }
        }
        fai_destroy(faidx);
    }
    return true;
}
