#include "utils.hpp"
#include "record.hpp"

#include <cstdio>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include <htslib/sam.h>


void parse_input_file(std::vector<std::string> & names, std::vector<region> & regions, std::ifstream & f_in)
{
    std::string id, chrom, type, zyg, region_string;
    while (f_in >> id >> chrom >> type >> zyg >> region_string)
    {
        std::string r;
        std::vector<std::string> temp_regions;
        std::stringstream ss(region_string);
        while(std::getline(ss, r, ','))
            temp_regions.push_back(r);

        if (temp_regions.size() != 1)
            continue;

        names.push_back(id);
        regions.emplace_back(temp_regions[0]);
    }
    std::cout << "[STATUS] Loaded " << names.size() << " variant regions from file." << std::endl;
}



int main(int argc, char** argv)
{
    if (argc != 6)
    {
        std::cerr << "Usage: ./extract_region_sequences <region_file> <output_dir> <output_list> <BAM> <MODE>" << std::endl;
        return -1;
    }

    std::string region_file = argv[1];
    std::string out_dir = argv[2];
    std::string out_file = argv[3];
    std::string bam_path = argv[4];
    std::string mode = argv[5];
    if (mode != "COMPLETE" && mode != "PARTIAL")
    {
	    std::cerr << "[ERROR] Invalid mode: " << argv[5] << ". Must be one of [COMPLETE, PARTIAL]." << std::endl;
	    return -1;
    }
    bool rescue_reads = (mode == "PARTIAL");

    // make sure target directory exists
    if (! std::filesystem::exists(out_dir))
        std::filesystem::create_directory(out_dir);

    if (!std::filesystem::is_directory(out_dir))
    {
        std::cerr << "[ERROR] Directory " << out_dir << " does not exist and could not be created" << std::endl;
        return -1;
    }



    // open output file
    std::ofstream f_out(out_file);
    if (!f_out.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << out_file << " for writing" << std::endl;
        return -1;
    }

    // Load and parse variant regions from input
    std::vector<std::string> variant_names;
    std::vector<region> variant_regions;

    std::ifstream f_in(region_file);
    if (!f_in.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << region_file << " for reading" << std::endl;
        return -1;
    }
    parse_input_file(variant_names, variant_regions, f_in);
    f_in.close();

    if (variant_names.size() == 0)
    {
        std::cerr << "[ERROR] No variant regions loaded." << std::endl;
        return -1;
    }

    std::cout << "[STATUS] Start iteration" << std::endl;
    // Iterate over variant regions
    for (int i = 0; i < variant_names.size(); ++i)
    {
        // std::cout << "\t" << i << std::endl;
        std::unordered_map<std::string, AlnRecord> records;

        std::string v_name = variant_names[i];
        region & r = variant_regions[i];
        std::string region_string = r.to_string();

        // open BAM file - within the loop in order to enable parallelization later
        samFile *bam_in = hts_open(bam_path.c_str(), "r");
        hts_idx_t *idx = sam_index_load(bam_in, (bam_path + ".bai").c_str());
        bam_hdr_t *bam_hdr = sam_hdr_read(bam_in);

        // iterate over records overlapping the region; -1 means no more data
        hts_itr_t *itr = sam_itr_querys(idx, bam_hdr, region_string.c_str());
        bam1_t *record = bam_init1();

        while (sam_itr_next(bam_in, itr, record) >= 0) {
            if (record->core.tid < 0) continue; // invalid chromosome id
            if (sam_hdr_tid2name(bam_hdr, record->core.tid) != r.r_name)
                break; // switched to wrong chromosome

            // basic qc: unmapped, duplicated or qc-failed (1540)
            if (record->core.flag & 1540)
                continue;

            if (record->core.qual < 1) // maybe this threshold should be set flexibly by the user?
                continue;

            // no need to process the read again if it is already present
            if (records.find(bam_get_qname(record)) == records.end())
                records[bam_get_qname(record)] = AlnRecord(record, bam_hdr); // if it is not present, create new AlnRecord() from bam1_t
        }

        // close and destroy everything
        sam_itr_destroy(itr);
        bam_destroy1(record);

        hts_idx_destroy(idx);
        bam_hdr_destroy(bam_hdr);
        sam_close(bam_in);

        // Unify the read orientations
        for (auto it = records.begin(); it != records.end(); ++it)
        {
            region start_region = {r.r_name, r.start, r.start + 20};
            region end_region = {r.r_name, r.end - 20, r.end};

            std::unordered_set<int8_t> start_orientations;
            std::unordered_set<int8_t> end_orientations;
            for (const aln_info & aln : it->second.get_alignments())
            {
                if (start_region.start <= aln.end && start_region.end >= aln.start && aln.r_name == start_region.r_name)
                    start_orientations.insert((aln.reverse ? -1 : 1));
                if (end_region.start <= aln.end && end_region.end >= aln.start && aln.r_name == start_region.r_name)
                    end_orientations.insert((aln.reverse ? -1 : 1));
            }
            if (start_orientations.size() == 1)
            {
                if (*start_orientations.begin() < 0)
                    it->second.flip();
            }
            else if (end_orientations.size() == 1)
            {
                if ( *end_orientations.begin() < 0)
                    it->second.flip();
            }
        }

        // extract sequence for the given region from all reads
        std::vector<std::string> sequences;
        std::vector<std::string> seq_names;
        // for every read
        // determine positions spanning the region
        int64_t min_len = r.end - r.start;
        std::unordered_map<std::string, std::tuple<int32_t, int32_t>> incomplete_intervals;
        for (auto it = records.begin(); it != records.end(); ++it)
        {
            std::tuple<int32_t, int32_t> interval{0,0};
            int8_t status = it->second.determine_region_interval(interval, r);

            if (status > 1)
                incomplete_intervals[it->first] = interval;
            if (status != 1)
                continue;


            // If a valid interval was found, extract the sequence
            if (it->second.get_primary_is_reverse())
                sequences.push_back(it->second.get_reverse_subseq(std::get<0>(interval), std::get<1>(interval)));
            else
                sequences.push_back(it->second.get_subseq(std::get<0>(interval), std::get<1>(interval)));
            seq_names.push_back(it->second.get_template_name() + "_" + std::to_string(std::get<0>(interval)) + "_" + std::to_string(std::get<1>(interval)));
            min_len = std::min(min_len, (int64_t) sequences[sequences.size() - 1].length());
        }

        // for all incomplete sequences, extract at most max_len bases; direction depending on whether start or end was estimated
        if (min_len > 1 && rescue_reads)
        {
            for (auto it : incomplete_intervals)
            {
                std::tuple<int32_t, int32_t> & interval = it.second;
                if (std::get<0>(interval) >= 0)
                {
                    int l = std::min(min_len, (int64_t) records[it.first].get_seq().length() - std::get<0>(interval) - 1);
                    if (l <= 200)
                        continue;
                    std::get<1>(interval) = std::get<0>(interval) + l;

                    if (records[it.first].get_primary_is_reverse())
                        sequences.push_back(records[it.first].get_reverse_subseq(std::get<0>(interval), std::get<1>(interval)));
                    else
                        sequences.push_back(records[it.first].get_subseq(std::get<0>(interval), std::get<1>(interval)));
                    seq_names.push_back(it.first + "_" + std::to_string(std::get<0>(interval)) + "_" + std::to_string(std::get<1>(interval)));
                }
                else if (std::get<1>(interval) >= 0)
                {
                    int l = std::min(min_len, (int64_t) std::get<1>(interval));
                    if (l <= 200)
                        continue;
                    std::get<0>(interval) = std::get<1>(interval) - l;

                    if (records[it.first].get_primary_is_reverse())
                        sequences.push_back(records[it.first].get_reverse_subseq(std::get<0>(interval), std::get<1>(interval)));
                    else
                        sequences.push_back(records[it.first].get_subseq(std::get<0>(interval), std::get<1>(interval)));
                    seq_names.push_back(it.first + "_" + std::to_string(std::get<0>(interval)) + "_" + std::to_string(std::get<1>(interval)));
                }
            }
        }

        // write sequences to fasta file and write file path + variant info (id) to out_file
        std::string fasta_name = out_dir + "/" + v_name + "_sequences.fa";
        std::ofstream fasta_out(fasta_name);
        if (!fasta_out.is_open())
        {
            std::cerr << "[ERROR] Could not open " << fasta_name << " for writing" << std::endl;
            return -1;
        }

        int w = 70;
        for (int j = 0; j < sequences.size(); ++j)
        {
            if (sequences[j].length() > 3500)
                std::cerr << "[WARNING] Variant " << variant_names[i] << ": Extracted sequence " << seq_names[j] << " exceeds expected size range (" << sequences[j].length() << ")" << std::endl;

            fasta_out << ">" <<  seq_names[j] << std::endl;
            int p = 0;
            while (p < sequences[j].length())
            {
                fasta_out << sequences[j].substr(p, w) << std::endl;
                p += w;
            }
        }

        fasta_out.close();

        // Determine phasing region
        int phasing_start_bin = ((int64_t) r.start / 100000);
        int phasing_end_bin = ((int64_t) r.end / 100000);

        std::string phasing_region_string = "";
        if (phasing_start_bin == phasing_end_bin)
            phasing_region_string = "chr1_" + std::to_string(phasing_start_bin * 100000) + "_" + std::to_string((phasing_start_bin+1) * 100000);
        else
            phasing_region_string = "chr1_" + std::to_string(phasing_start_bin * 100000 + 50000) + "_" + std::to_string((phasing_start_bin+1) * 100000 + 50000);

        f_out << fasta_name << "\t" << v_name << "\t" << phasing_region_string << std::endl;
    }

    f_out.close();

    return 0;
}
