#ifndef RECORD_HEADER
#define RECORD_HEADER

#include "utils.hpp"

#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>

#include <htslib/sam.h>

/*
 @brief Datatype for storing junction classifications
 */
enum junction_type {
  del,
  inv,
  ins,
  duplication,
  trans,
  clip,
  mismatch_region,
  unknown
};

/**
 @brief Datatype for handling secondary alignments.
*/
struct aln_info {
    std::string r_name;                            // Name of chromosome
    uint32_t start;                                // First aligned position on the chromosome
    uint32_t end;                                  // Last aligned position on the chromosome
    uint32_t seq_start;                            // First aligned position on the read
    uint32_t seq_len;                              // Length of the aligned sequence
    std::vector<std::tuple<char, uint32_t>> cigar; // CIGAR operations for the record, given by tuples of type and length

    bool reverse;                                  // Strand information
    uint8_t mapq;                                  // Mapping quality
    bool primary;                                  // Is primary alignment
};




/**
 @brief Class for storing / handling information extracted from bam1_t, in order
 to avoid those pesky pointers and potential memory allocation issues and make secondary alignments more accessible.
 Might include compressed sequence storage and functions for breakpoint extraction in the future.
*/
class AlnRecord {

    protected:
    std::string template_name;
    std::string seq;
    std::vector<uint8_t> qual;

    bool primary_is_reverse;

    std::vector<aln_info> alignments;

    /**
     @private
     @brief Private function for converting read sequence to string
     @param b pointer to alignment as read by htslib
     @return String of nucleotides (ACTG).
    */
    inline std::string extract_seq(bam1_t * b);

    /**
     @private
     @brief Private function for extracting CIGAR string in a more readable format
     @param b pointer to alignment as read by htslib
     @return vector of tuples consisting of CIGAR operation and length
    */
    inline std::vector<std::tuple<char, uint32_t>> extract_cigar_operations(bam1_t * b);

    /**
     @private
     @brief Private function for converting array of qualities into a vector
     @param b pointer to alignment as read by htslib
     @return vector of quality values
    */
    inline std::vector<uint8_t> extract_qual(bam1_t * b);

    /**
     @private
     @brief Read SA tag and extract relevant information (pos, orientation, CIGAR string)
     @param b pointer to alignment as read by htslib
     @return NULL. Extracted secondary alignments are stored in member variable'.
    */
    void extract_secondary_alignments(bam1_t * b);

    /**
     @private
     @brief Sort all alignments based on their starting position on the read sequence.
    */
    inline void sort_alignments();

    inline void clean_alignments();

    /**
     @private
     @brief The first and last position of a subsequence - determined via CIGAR string - are checked for validity and
     adjusted for the orientation of the alignment if necessary.
     @param first_pos Interval start (index of first base to be extracted from read sequence). Passed by reference and possibly changed by this function.
     @param end_pos Interval end (index of last base to be extracted from read sequence). Passed by reference and possibly changed by this function.
     @param estimated_start Indicates whether the first position was estimated from clipped sequence
     @param estimated_end Indicates whether end position was estimated from clipped sequence
     @param aln Alignment from which the positions were inferred.
     @return Boolean value indicating whether the interval is valid
    */
    inline bool adjust_sequence_interval(int32_t & first_pos, int32_t & last_pos,  bool & estimated_start, bool & estimated_end, const aln_info & aln);

    /**
     @private
     @brief Go through all alignments, adjust CIGAR strings based on orientation and determine read start position and alignment length on read
     @param temp_aln Alignment to adjust
    */
    inline void finalize_alignment(aln_info & temp_aln);

    /**
     @brief Function to determine whether an alignment is reversed with respect to the read's primary alignment (and thus the stored sequence).
     @param aln Alignment in question
     @return True if the alignment is reversed
    */
    inline bool is_reversed_wrt_primary(const aln_info & aln);


    public:
    /**
     @brief Default constructor.
    */
    AlnRecord();

    /**
     @brief Constructor for extracting information from bam1_t object (htslib).
     @param b pointer to alignment as read by htslib
     @param hdr pointer to htslib header object
     @param Use MD tag to distinguish matches and mismatches for inversion rescue
    */
    AlnRecord(bam1_t * b, sam_hdr_t * hdr);


    /**
     @brief Constructor for mock record containing multiple alignments, only used for testing
     @param sequence read sequence
     @param alignments vector of aln_info objects containing primary and secondary alignments
    */
    AlnRecord(std::string sequence, std::vector<aln_info> alignments);


    /**
     @brief Determine a read interval overlapping a region of variation
     @param interval The return tuple. Will contain the start and end position of the Interval.
     @param r Region defined by the Node representing the variant region
     @return integer indicating success (1), no valid region (0) or faulty position (-1)
    */
    // alternative definition used for the simplified simulated data in POA method benchmark
    int8_t determine_region_interval(std::tuple<int32_t, int32_t> & interval, region r);


    /**
     @brief Reverse the orientation of the read alignments and the stored sequence.
    */
    void flip();


    /**
     @brief Function for accessing the read sequence (read only)
     @return Const reference to nucleotide string
    */
    const std::string & get_seq();

    /**
     @brief Function for obtaining a subsequence [start, end] of the read sequence
     @param start Starting position of segment to be extracted
     @param end End position of segment to be extracted
     @return Nucleotide string
    */
    std::string get_subseq(uint32_t start, uint32_t end);

    /**
     @brief Function for obtaining a reverse complemented copy of the sequence
     @return Nucleotide string of reversed sequence
    */
    std::string get_reverse_seq();

    /**
     @brief Function for obtaining a reversed subsequence [start, end] of the read sequence.
     Note that this function first extracts the interval [start, end] from the original
     sequence and then creates the reverse complement of this fragment.
     @param start Starting position of segment to be extracted
     @param end End position of segment to be extracted
     @return Nucleotide string of reversed segment
    */
    std::string get_reverse_subseq(uint32_t start, uint32_t end);

    /**
     @brief Function for accessing alignments
     @return Vector of aln_info objects containing all alignments for the read.
    */
    const std::vector<aln_info> & get_alignments();

    /**
     @brief Function for accessing base qualities of read (read only)
     @return Vector of integer quality values of same length as read sequence
    */
    const std::vector<uint8_t> & get_base_qual();

    /**
     @brief Function for getting template name
     @return Name of template to which an alignment belongs
    */
    std::string get_template_name();

    /**
     @brief Getter function for reverse_status of primary alignment.
     */
    bool get_primary_is_reverse();

    /**
     @brief Process one CIGAR operation, moving either reference or read position forward
     @param operation Tuple representing the next CIGAR operation (see extract_cigar_operations)
     @param ref_pos Position on the reference as determined by CIGAR operations up to this point, initialized as alignment start position.
     Given as reference, will be changed in the function depending on CIGAR operation.
     @param read_pos Position on the read as determined by CIGAR operations up to this point.
     Given as reference, will be changed in the function depending on CIGAR operation.
     @return Integer representing whether the read (+1), the reference (-1) or neither (0) was consumed
    */
    static int consume_cigar_op(const std::tuple<char, uint32_t> & operation, uint32_t & ref_pos, uint32_t & read_pos);


    static std::unordered_map<char, char> base_complement; // Dictionary for fast complementing of bases
};

#endif
