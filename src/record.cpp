#include "record.hpp"
#include "utils.hpp"

#include <cassert>
#include <cstdlib>
#include <ranges>
#include <algorithm>
#include <regex>
#include <iostream>

#include "htslib/sam.h"

std::unordered_map<char, char> AlnRecord::base_complement = std::unordered_map<char, char>{
        {'A', 'T'},
        {'T', 'A'},
        {'C', 'G'},
        {'G', 'C'},
        {'N', 'N'}
    };


AlnRecord::AlnRecord()
{}


AlnRecord::AlnRecord(bam1_t * record, sam_hdr_t * hdr)
{
    this->template_name = bam_get_qname(record);
    this->seq = extract_seq(record);
    this->qual = extract_qual(record);
    this->primary_is_reverse = bam_is_rev(record);

    // extract information from the given record (primary alignment)
    aln_info record_alignment;
    record_alignment.r_name = sam_hdr_tid2name(hdr, record->core.tid);
    record_alignment.start = record->core.pos;
    record_alignment.end = bam_endpos(record);
    record_alignment.reverse = bam_is_rev(record);
    record_alignment.mapq = record->core.qual;
    record_alignment.cigar = extract_cigar_operations(record);

    //////////////////////////////////////////////////////////////////////////
    // extract secondary alignments of the read from tags
    extract_secondary_alignments(record);

    // if no other alignments are present, designate this one as primary by default,
    // even if the flag says something different; 2304: not primary, or supplementary
    record_alignment.primary = (!((record->core.flag & 2304) > 0)) || this->alignments.size() == 0;

    // assign the correct alignment to be primary and store the read sequence in the orientation indicated by the primary alignment
    if (!record_alignment.primary)
    {
        this->alignments[0].primary = true; // by convention, the first supplementary alignment is primary
        this->primary_is_reverse = this->alignments[0].reverse;

        if (record_alignment.reverse != this->primary_is_reverse)
            this->seq = get_reverse_seq();
    }
    this->alignments.push_back(record_alignment);

    for (auto & temp_aln : this->alignments)
        finalize_alignment(temp_aln);

    // sort alignments in the order they occur on the read
    if (this->alignments.size() > 1) {
        sort_alignments();
        clean_alignments();
    }
}


inline std::string AlnRecord::extract_seq(bam1_t* record){
    std::string seq = "";
    for (uint32_t i = 0; i < record->core.l_qseq; ++i)
        seq += seq_nt16_str[bam_seqi(bam_get_seq(record), i)];
    return seq;
}


inline std::vector<std::tuple<char, uint32_t>> AlnRecord::extract_cigar_operations(bam1_t * record)
{
    std::vector<std::tuple<char, uint32_t>> cigar;
    for (uint16_t i = 0; i < record->core.n_cigar; ++i)
        cigar.push_back(
            std::tuple<char, uint32_t>(
                bam_cigar_opchr(bam_get_cigar(record)[i]),
                bam_cigar_oplen(bam_get_cigar(record)[i])
            )
        );
    return cigar;
}


inline std::vector<uint8_t> AlnRecord::extract_qual(bam1_t * record)
{
    std::vector<uint8_t> qualities;
    for (uint32_t i = 0; i < record->core.l_qseq; ++i)
        qualities.push_back(bam_get_qual(record)[i]);
    return qualities;
}


void AlnRecord::extract_secondary_alignments(bam1_t * record)
{
    // get tag
    const uint8_t * sa_start = bam_aux_get(record, "SA");
    if (sa_start == NULL)
        return;

    std::string sa_string = bam_aux2Z(sa_start);
    // remove trailing delimiter if necessary
    if (sa_string[sa_string.length() - 1] == ';')
        sa_string = sa_string.substr(0, sa_string.length() - 1);

    // extract information
    std::vector<std::string> alignment_strings = split_string(sa_string, ';');
    for (auto & s : alignment_strings)
    {
        std::vector<std::string> aln_specs = split_string(s, ',');
        aln_info temp_aln;
        temp_aln.primary = false;

        // extract the low hanging fruit
        temp_aln.r_name = aln_specs[0];
        temp_aln.start = std::stoi(aln_specs[1]) - 1; // convert to 0-based
        temp_aln.reverse = aln_specs[2] == "-";
        temp_aln.mapq = std::stoi(aln_specs[4]);

        // convert CIGAR string to vector of tuples corresponding to operations
        std::regex op_regex("[0-9]+[MIDNSHP=XB]");
        auto op_begin = std::sregex_iterator(aln_specs[3].begin(), aln_specs[3].end(), op_regex);
        auto op_end = std::sregex_iterator();

        for (; op_begin != op_end; ++op_begin)
        {
            std::string op = (*op_begin).str();
            uint32_t op_len = std::stoi(op.substr(0, op.length() - 1)); // operation length
            char op_type = op[op.length() - 1]; // last position in string is the type of CIGAR operation
            temp_aln.cigar.push_back(std::tuple<char, uint32_t>{op_type, op_len});
        }

        this->alignments.push_back(temp_aln);
    }
}


inline void AlnRecord::finalize_alignment(aln_info & temp_aln)
{
    uint8_t min_clip_len = 100;
    uint8_t min_indel_len = 30;
    // infer position and length of alignment on the read sequence and end position of alignment on the reference from CIGAR string
    temp_aln.seq_start = (std::get<0>(temp_aln.cigar[0]) == 'S') ? std::get<1>(temp_aln.cigar[0]) : 0;
    temp_aln.seq_len = 0;
    temp_aln.end = temp_aln.start;
    for (const auto & op : temp_aln.cigar)
        consume_cigar_op(op, temp_aln.end, temp_aln.seq_len); // update temp_aln.end and temp_aln.seq_len with each CIGAR operation
    temp_aln.seq_len -= temp_aln.seq_start;
    if (std::get<0>(temp_aln.cigar[temp_aln.cigar.size() - 1]) == 'S')
        temp_aln.seq_len -= std::get<1>(temp_aln.cigar[temp_aln.cigar.size() - 1]);
    return;
}


inline void AlnRecord::clean_alignments()
{
    for (auto it = this->alignments.begin(); it != this->alignments.end(); )
    {
        bool removed_first = false;
        for (auto it1 = it + 1; it1 != this->alignments.end(); )
        {
            // it must start before it1; need only compare the ends
            if (it->seq_start + it->seq_len > it1->seq_start)
            {
                float overlap = (float) std::min(it1->seq_len, it->seq_start + it->seq_len - it1->seq_start) / (float) std::min(it->seq_len, it1->seq_len);
                if (overlap > 0.25)
                {
                    if (it->mapq < it1->mapq || (it->mapq == it1->mapq && it->seq_len < it1->seq_len))
                    {
                        removed_first = true;
                        it = this->alignments.erase(it);
                        break;
                    } else {
                        it1 = this->alignments.erase(it1);
                        continue;
                    }
                }
            }
            ++it1;
        }
        if (!removed_first)
            ++it;
    }
}


inline void AlnRecord::sort_alignments()
{
    /**
    Comparison function for sort() causes alignments to be switched if the second alignment starts earlier on the read.
    If two alignments have the same position, the primary alignment will come first.
    */
    std::sort(
        this->alignments.begin(),
        this->alignments.end(),
        [&](aln_info & a, aln_info & b){
            uint32_t a_begin = (is_reversed_wrt_primary(a) ? this->seq.length() - (a.seq_start + a.seq_len) : a.seq_start);
            uint32_t b_begin = (is_reversed_wrt_primary(b) ? this->seq.length() - (b.seq_start + b.seq_len) : b.seq_start);
            return a_begin < b_begin || (a_begin == b_begin && a.primary);
        });
};


int8_t AlnRecord::determine_region_interval(std::tuple<int32_t, int32_t> & interval, region r)
{
    int8_t status = -1;

    int32_t start_distance = this->seq.length();
    int32_t end_distance = this->seq.length();
    int32_t final_start = -1;
    int32_t final_end = -1;

    // attempt to find a matching position on both sides of the region
    for (auto & aln: this->alignments)
    {
        int32_t start_position = -1;
        int32_t end_position = -1;

        uint32_t ref_pos = aln.start;
        uint32_t read_pos = 0;

        // after every operation, check whether the reference position is closer to the start / end of the region
        // and whether it within the region of interest;
        for (auto & op : aln.cigar)
        {
            consume_cigar_op(op, ref_pos, read_pos);

            if (std::get<0>(op) != 'M' && std::get<0>(op) != '=')
                continue;

            int64_t d_start = (int64_t) ref_pos - r.start;
            int64_t start_correction = d_start > 0 ? std::min((int64_t) std::get<1>(op), (int64_t) d_start) : 1;

            int64_t d_end = (int64_t) ref_pos - r.end;
            int64_t end_correction = d_end <= 0 ? 1 : std::min((int64_t) std::get<1>(op), (int64_t) d_end);

            if (std::abs(d_start - start_correction) <= std::abs(start_distance) && d_start - start_correction >= -250)
            {
                start_position = read_pos - start_correction;
                start_distance = d_start - start_correction;
            }

            if (std::abs(d_end - end_correction) <= std::abs(end_distance) && d_end - end_correction <= 250)
            {
                end_position = read_pos - end_correction;
                end_distance = d_end - end_correction;
            }
        }

        if (is_reversed_wrt_primary(aln))
        {
            if (start_position >= 0)
                start_position = this->seq.length() - start_position - 1;
            if (end_position >= 0)
                end_position = this->seq.length() - end_position - 1;
            std::swap(start_position, end_position);
        }
        if (start_position >= 0 && start_distance >= -250 && start_distance <= 250)
            final_start = (final_start == -1) ? start_position : std::min(final_start, start_position);

        if (end_position >= 0 && end_distance <= 250 && end_distance >= -250)
            final_end = (final_end == -1) ? end_position : std::max(final_end, end_position);
    }

    interval = {final_start, final_end};

    bool valid_positions = final_start >= 0 && final_end < this->seq.length() && final_end > final_start;
    if (valid_positions)
        status = 1;
    else if (final_start >= 0 && final_end == -1)
        status = 2;
    else if (final_end >= 0 && final_start  == -1)
        status = 4;

    return status;
}


inline bool AlnRecord::adjust_sequence_interval(int32_t & first_pos, int32_t & last_pos, bool & estimated_start, bool & estimated_end, const aln_info & aln)
{
    // Invalid positions are detected here
    bool valid_positions = first_pos >= 0 && last_pos > first_pos && last_pos < this->seq.length();

    // Whether an interval needs to be flipped no longer depends on the primary alignments,
    // as the sequence is stored in the orientation indicated by aln.reverse
    bool flip_alignment = is_reversed_wrt_primary(aln);

    if (valid_positions && flip_alignment)
    {
        first_pos = this->seq.length() - first_pos - 1;
        last_pos = this->seq.length() - last_pos - 1;
        std::swap(first_pos, last_pos);
        std::swap(estimated_start, estimated_end);
    }
    return valid_positions;
}


void AlnRecord::flip()
{
    // alignment start/end position on the read and its indicated alignment orientation
    for (auto & aln : this->alignments)
        aln.reverse = !aln.reverse;
    this->primary_is_reverse = !this->primary_is_reverse;
}


inline bool AlnRecord::is_reversed_wrt_primary(const aln_info & aln)
{
    return (aln.reverse != this->primary_is_reverse);
}


const std::string & AlnRecord::get_seq()
{
    return this->seq;
}


std::string AlnRecord::get_subseq(uint32_t start, uint32_t end)
{
    if (end >= this->seq.length() || end <= start)
        std::cerr << "[ERROR] " << this->seq.length() << "; [" << start << "," << end << "]..." << std::endl;
    assert(end > start);
    assert(end < this->seq.length());

    auto v = this->seq | std::ranges::views::take(end+1) | std::ranges::views::drop(start);

    return std::string(v.begin(), v.end());
}


std::string AlnRecord::get_reverse_seq()
{
    std::string rev_seq;
    rev_seq.resize(this->seq.length());

    auto v = this->seq | std::ranges::views::reverse;
    std::ranges::transform(v, rev_seq.begin(), [](char & c){return AlnRecord::base_complement[c];});

    return rev_seq;
}


std::string AlnRecord::get_reverse_subseq(uint32_t start, uint32_t end)
{
    assert(end > start);
    assert(end < this->seq.size());

    std::string rev_seq = "";
    rev_seq.resize(end - start + 1);

    auto v = this->seq | std::ranges::views::take(end+1) | std::ranges::views::drop(start) | std::ranges::views::reverse;
    std::ranges::transform(v, rev_seq.begin(), [](char c){return AlnRecord::base_complement[c];});

    return rev_seq;
}


const std::vector<aln_info> & AlnRecord::get_alignments()
{
    return this->alignments;
}


const std::vector<uint8_t> & AlnRecord::get_base_qual()
{
    return this->qual;
}


std::string AlnRecord::get_template_name()
{
    return this->template_name;
}


int AlnRecord::consume_cigar_op(const std::tuple<char, uint32_t> &operation, uint32_t &ref_pos, uint32_t &read_pos)
{
    const int8_t & cigar_op = bam_cigar_table[std::get<0>(operation)];

    // depending on operation, move either ref_pos or read_pos by specified length
    int return_val = 0;
    if (bam_cigar_type(cigar_op) & 1) // consumes read (query)
    {
        read_pos += std::get<1>(operation);
        return_val -= 1;
    }
    if (bam_cigar_type(cigar_op) & 2){    // consumes reference, either forwards or backwards
        ref_pos += std::get<1>(operation);
        return_val += 1;
    }
    return return_val;
}


bool AlnRecord::get_primary_is_reverse()
{
    return this->primary_is_reverse;
}
