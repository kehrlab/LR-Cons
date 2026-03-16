#ifndef UTILHEADER
#define UTILHEADER

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <ranges>


/**
 @brief Data structure holding alignment parameters for SPOA
*/
struct aln_params {
    bool global;
    int8_t match;
    int8_t mismatch;
    int8_t open;
    int8_t ext;
};


/**
 @brief Datatype for handling genomic regions.
*/
struct region {
    std::string r_name; // name of chromosome
    uint32_t start;     // start of region
    uint32_t end;       // end of region

    region(){}

    region(std::string r_name, uint32_t start, uint32_t end) : r_name(r_name), start(start), end(end)
    {}

    region(std::string region_string)
    {
        int pos1 = region_string.find(":");
        int pos2 = region_string.find("-", pos1);

        r_name = region_string.substr(0, pos1);
        ++pos1;

        start = std::stoi(region_string.substr(pos1, pos2 - pos1));
        end = std::stoi(region_string.substr(pos2 + 1));
    }

    std::string to_string() {
        return (r_name + ":" + std::to_string(start) + "-" + std::to_string(end));
    }
};





// Collection of small helper functions used throughout the algorithm or for debugging

/**
 @brief Split a string into substrings at a given delimiter. Intended for use on SA tag in bam files.
 NOTE: For simplicity, the string must NOT end on a delimiter.
 @param s string to be split
 @param delim character at which to split the string
 @return vector of substrings
*/
inline std::vector<std::string> split_string(std::string s, char delim)
{
    std::vector<std::string> substrings;
    for (const auto substring : std::views::split(s, delim))
        substrings.push_back(std::string(std::string_view(substring)));

    return substrings;
}

#endif
