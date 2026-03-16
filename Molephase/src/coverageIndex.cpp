# include "seqan/seq_io.h"
# include "seqan/sequence.h"
# include "parse.h"
# include "call_index.h"
# include "coverageIndex.h"
using namespace seqan;

int buildCoverageIndex(coverageIndex_t & coverageIndex, std::string & readfile1_name, std::string & readfile2_name, uint32_t k, bool & paired){
  // load reads
  StringSet<Dna5String> reads1;
  StringSet<Dna5String> reads2;
  loadReads(reads1, reads2, readfile1_name, readfile2_name, paired);

  if(length(reads1)==0){
    return 1;
  }

  // fill index
  uint64_t index_size;
  if(paired){
    index_size = length(reads1)*(length(reads1[0])-k+1)*2;
  }else{
    index_size = lengthSum(reads1)-(length(reads1)*(k+1));
  }

  coverageIndex.resize(index_size);

  fillCoverageIndex(reads1, reads2, coverageIndex, k, paired);

  //filter index

  return 0;
}

void loadReads(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, std::string & readfile1_name, std::string & readfile2_name, bool & paired){
  // open files
  StringSet<CharString> ids;
  SeqFileIn file1(toCString(readfile1_name));
  readRecords(ids, reads1, file1);
  // read from readfiles
  if(paired){
    SeqFileIn file2(toCString(readfile2_name));
    readRecords(ids, reads2, file2);
  }
  // clear(ids);
  return;
}

void fillCoverageIndex(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex_t & coverageIndex, uint32_t k, bool & paired){
  fillIndexFromReads(reads1, coverageIndex, k);
  if(paired){
    fillIndexFromReads(reads2, coverageIndex, k);
  }
  return;
}

void fillIndexFromReads(StringSet<Dna5String> & reads, coverageIndex_t & coverageIndex, uint32_t k){
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  Iterator<Dna5String>::Type itread;
  hashId_t hashId;

  for (TStringSetIterator itreads = begin(reads); itreads!=end(reads); ++itreads){ //  for every read in reads

    // add first kmer to index
    Dna5String kmer = prefix(*itreads,k);
    hash_t hash(kmer);
    hashId = hash.bigger();
    coverageIndex.insert(hashId);

    // iterate over read to extract kmers
    for(itread=begin(*itreads)+k; itread<end(*itreads);itread++){
      // roll hash value
      hash.roll(*itread);
      // add kmer to Index
      hashId = hash.bigger();
      coverageIndex.insert(hashId);
    }
  }

  return;
}
