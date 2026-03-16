#ifndef PARSE_H
#define PARSE_H

// # include <seqan/seq_io.h>
// # include <seqan/sequence.h>
# include "seqan/arg_parse.h"

struct phaseOptions{
  std::string readfile1;
  std::string readfile2;
  std::string output_prefix;

  bool paired;
  uint32_t k;
  uint32_t minCov;
  float covRat;
  std::string reference;
  std::string chrom;
  int32_t start;
  int32_t end;

  phaseOptions():
  minCov(1), covRat(1.5), k(77), paired(false)
  {}
};

seqan::ArgumentParser::ParseResult parsePhaseCommandLine(phaseOptions & options, int argc, char const ** argv);
void printPhaseParseResult(phaseOptions & options, seqan::ArgumentParser::ParseResult res);

#endif
