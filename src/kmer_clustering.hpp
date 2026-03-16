#ifndef CLUST_HDR
#define CLUST_HDR

#include "kmeans/kmeans.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};



inline uint16_t hash_kmer(std::string kmer) {
    uint16_t hash = 0;
    std::vector<uint16_t> base_vals = {1, 4, 16, 64, 256, 1024, 2048};
    for (uint8_t i = 0; i < 7; ++i)
        hash += base_vals[i] * seq_nt4_table[kmer[i]];
    return hash;
}


inline std::vector<std::vector<uint8_t>> count_kmers(std::vector<std::string> strings) {
    std::vector<std::vector<uint8_t>> cnt(16384, std::vector<uint8_t>(strings.size(), 0));
    for (uint8_t j = 0; j < strings.size(); ++j)
    {
        for (uint32_t i = 0; i < (int64_t) strings[j].length() - 7; ++i)
        {
            uint16_t hash = hash_kmer(strings[j].substr(i, 7));
            if (cnt[hash][j] <= 100) // anything that is too frequent can be excluded
                cnt[hash][j]++;
        }
    }
    return cnt;
}


inline float mean(const std::vector<uint8_t> & cnt) {
    uint16_t s = 0;
    for (auto & c : cnt)
    {
        if (c == 100)
            return -3;
        s += c;
    }

    if (s == 0)
        return -2;

    if (s < 0.15 * cnt.size()) // not worth further investigation, not present in enough reads
        return -1.;

    return (float) s / cnt.size();
}


inline float sd(float m, const std::vector<uint8_t> & cnt) {
    float s = 0.;
    for (auto & c : cnt)
        s += (m-c) * (m-c);
    return s / (cnt.size() - 1);
}


inline bool reduce_kmers(std::vector<std::vector<uint8_t>> & cnt) {
    std::vector<float> sds(cnt.size(), 0.);
    std::vector<int> indices(cnt.size(), 0);
    for (int i = 0; i < sds.size(); ++i)
        indices[i] = i;

    int insufficient = 0;
    int missing = 0;
    int overrepresented = 0;
    for (uint16_t i = 0; i < cnt.size(); ++i)
    {
        float m = mean(cnt[i]);
        float s = 0;

        if (m < 0) {
            if (m == -2)
                ++missing;
            else if (m == -1)
                ++insufficient;
            else
                ++overrepresented;
        } else {
            s = sd(m, cnt[i]);
        }
        sds[i] = s;
    }

    std::sort(indices.begin(), indices.end(), [&](int a, int b) -> bool {
        return sds[a] > sds[b];
    });

    #ifdef VERBOSE
    std::cout << "Number of k-mers: " << 16384 - missing << std::endl;
    std::cout << "K-Mers with insufficient occurrence: " << insufficient << std::endl;
    #endif

    size_t n = std::min(100, 16384 - missing - insufficient);

    if (n < 10)
    {
	    std::cerr << "[INFO] Insufficient number of remaining k-mers." << std::endl;
	    return false;
    }

    #ifdef VERBOSE
    std::cout << "Reduce to " << n << " k-mers with highest variance across sequences." << std::endl;
    #endif

    size_t n_reads = cnt[0].size();
    std::vector<std::vector<uint8_t>> reduced_cnt(n_reads, std::vector<uint8_t>(n));

    #ifdef VERBOSE
    std::cout << "Number of sequences: " << reduced_cnt.size() << std::endl << std::flush;
    #endif

    for (uint16_t i = 0; i < n; ++i)
        for (uint16_t j = 0; j < cnt[indices[i]].size(); ++j)
            reduced_cnt[j][i] = cnt[indices[i]][j];

    cnt = reduced_cnt;
    return true;
}


inline std::vector<int> k_clust_kmers(std::vector<std::vector<uint8_t>> cnt, uint8_t k)
{
    int n_obs = cnt.size();
    int n_dim = cnt[0].size();

    if (n_dim == 0)
        return std::vector<int>();

    std::vector<float> matrix(n_obs * n_dim);
    for (int i = 0; i < cnt.size(); ++i)
        for (int j = 0; j < cnt[i].size(); ++j)
            matrix[i * n_dim + j] = cnt[i][j];

    kmeans::SimpleMatrix<int, float> kmat(n_dim, n_obs, matrix.data());
    auto res = kmeans::compute(
        kmat,
        kmeans::InitializeKmeanspp<int, float, int, float>(),
        kmeans::RefineLloyd<int, float, int, float>(),
        (int) k
    );

    return res.clusters;
}

#endif
