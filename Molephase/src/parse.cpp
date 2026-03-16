# include <iostream>
# include <fstream>
# include "seqan/seq_io.h"
# include "seqan/sequence.h"
# include "seqan/arg_parse.h"
# include "parse.h"
// # include "split_io.h"
// # include "functions.h"
# include <time.h>

using namespace seqan;

seqan::ArgumentParser::ParseResult parsePhaseCommandLine(phaseOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("molephase");

    // Define arguments.

    // Define Options
    addOption(parser, seqan::ArgParseOption(
      "1", "readfile1", "Path to fastq read file.",
      seqan::ArgParseArgument::OUTPUT_FILE, "output_file[OUT]"));
    setRequired(parser, "1");
    addOption(parser, seqan::ArgParseOption(
      "2", "readfile2", "Path to fastq file containing second read in pair.",
      seqan::ArgParseArgument::OUTPUT_FILE, "output_file[OUT]"));
    addOption(parser, seqan::ArgParseOption(
      "o", "output_prefix", "Prefix of the phased output fastq file(s): <Prefix>.hap{1,2}.fq / <Prefix>.hap{1,2}.1.fq,<Prefix>.hap{1,2}.2.fq",
      seqan::ArgParseArgument::OUTPUT_FILE, "output_file[OUT]"));
    setRequired(parser, "o");
    addOption(parser, seqan::ArgParseOption(
      "k", "k-mer_size", "k-mer size used for kmer coverage analysis.",
    seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "77");
    setMinValue(parser, "k", "15");
    setMaxValue(parser, "k", "95");
    addOption(parser, seqan::ArgParseOption(
      "r", "reference", "Name of the reference file.",
      seqan::ArgParseArgument::INPUT_FILE, "input_file[IN]"));
    setRequired(parser, "r");
    addOption(parser, seqan::ArgParseOption(
      "m", "min-cov", "Minimium coverage of heterozygous k-mers.",
    seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "1");
    addOption(parser, seqan::ArgParseOption(
      "n", "cov-ratio", "Maximum ratio of ref-kmer and alt-kmer coverage.",
    seqan::ArgParseArgument::DOUBLE, "float"));
    setDefaultValue(parser, "n", "1.5");
    addOption(parser, seqan::ArgParseOption(
      "c", "chrom_name", "Contig name of the genomic region of interest.",
      seqan::ArgParseArgument::STRING, "string"));
    setRequired(parser, "c");
    addOption(parser, seqan::ArgParseOption(
      "s", "start", "start position of genomic region of interest.",
      seqan::ArgParseArgument::INTEGER, "unsigned"));
    setRequired(parser, "s");
    addOption(parser, seqan::ArgParseOption(
      "e", "end", "start position of genomic region of interest.",
      seqan::ArgParseArgument::INTEGER, "unsigned"));
    setRequired(parser, "e");


    seqan::addUsageLine(parser," -1 readfile -o [output_prefix] -r [reference.fa] -c [contig] -s [start] -e [end] [OPTIONS]");
    setShortDescription(parser, "Phase reads using heterozygous k-mer coverage.");
    // TODOTODO uncomment when connected to git repository
    // setVersion(parser, VERSION);
    // setDate(parser, DATE);
    addDescription(parser,
               "Takes reads from a region of interest. "
               "Reads are phased based on kmer distributions and pseudo-SNPs.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract argument and option values.
    getOptionValue(options.readfile1, parser, "1");
    getOptionValue(options.readfile2, parser, "2");
    options.paired = isSet(parser, "2");
    getOptionValue(options.output_prefix, parser, "o");
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.minCov, parser, "m");
    getOptionValue(options.covRat, parser, "n");
    getOptionValue(options.reference, parser, "r");
    getOptionValue(options.chrom, parser, "c");
    getOptionValue(options.start, parser, "s");
    getOptionValue(options.end, parser, "e");

    return seqan::ArgumentParser::PARSE_OK;
}

void printPhaseParseResult(phaseOptions & options, seqan::ArgumentParser::ParseResult res){
  std::cerr <<'\n' << "readfile1               \t" << options.readfile1 << '\n';
  if(options.paired){
    std::cerr << "readfile2               \t" << options.readfile2 << '\n';
  }
  std::cerr
    << "output_prefix           \t" << options.output_prefix << '\n'
    << "k                       \t" << options.k << '\n'
    << "minCov                  \t" << options.minCov << '\n'
    << "covRat                  \t" << options.covRat << '\n'
    << "reference               \t" << options.reference << '\n'
    << "chrom                   \t" << options.chrom << '\n'
    << "start                   \t" << options.start << '\n'
    << "end                     \t" << options.end << '\n'
    << '\n';
}
