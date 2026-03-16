#ifndef COVERAGEINDEX_H
#define COVERAGEINDEX_H

# include "parse.h"
# include "call_index.h"
# include <algorithm>


using namespace seqan;

struct coverageIndexBucket_t{
  hashId_t id;
  uint32_t coverage;

  coverageIndexBucket_t(){
    coverage = 0;
  }

  bool operator<(coverageIndexBucket_t const& rhs){
      return coverage<rhs.coverage;
  }
};

struct coverageIndex_t{
  std::vector<coverageIndexBucket_t> table;
  uint64_t tableSize;

  coverageIndex_t(){}

  void clear(){
    tableSize = 0;
    table.clear();
    return;
  }

  void resize(uint64_t tableSize_p){
    tableSize = tableSize_p;
    table.resize(tableSize);
    return;
  }


  bool getpos(hashId_t & hashId, uint64_t & pos){ // look up hashId in index. modifies pos and returns whether kmer is present
    pos = hashId.i2 % tableSize;
    for(int d=0; d<500; d++){
      // first try in loop
      if(table[pos].id == hashId){
        return true;
      } //return if correct bucket found
      if(table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0){ // if empty field hash is not in index
        return false;
      }
      pos = (pos^(hashId.i1 >> (d & 31))) % tableSize; //first shift in loop
      // second try in loop
      if(table[pos].id == hashId){
        return true;
      } //return if correct bucket found
      if(table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0){ // if empty field hash is not in index
        return false;
      }
      pos = (pos^(hashId.i3 << (d & 31))) % tableSize; // second shift in loop
    }
    std::cerr << "Warning: kmer index is overfilled.\n";
    return false;
  }

  bool insert(hashId_t & hashId){ // insert kmer into hashtable or increment coverage for kmer
    uint64_t pos;
    pos = hashId.i2 % tableSize;
    for(int d=0; d<500; d++){
      // first try in loop
      table[pos].id = (table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        table[pos].coverage++;
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId.i1 >> (d & 31)) )% tableSize; // first shift in loop
      // second try in loop
      // #pragma omp atomic
      table[pos].id = (table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        // #pragma omp atomic
        table[pos].coverage++;
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId.i3 << (d & 31))) % tableSize; // second shift in loop
    }
    std::cerr << "Warning: kmer index is overfilled.\n";
    return false;
  }

  bool lookup(hashId_t & hashId, uint32_t & coverage){
    uint64_t pos;

    if(this->getpos(hashId, pos)){
      coverage = table[pos].coverage;
      return true;
    }

    coverage = 0;
    return false;
  }

  uint32_t median(std::string & reference, phaseOptions & options){
    // lookup reference kmer coverages and add to coverage list if coverage not 0
    std::vector<uint32_t> coverages; //vector of non zero coverages
    uint32_t coverage;
    coverages.reserve(options.end-options.start);

    hash_t hash(reference.substr(1,options.k));
    hashId_t hashId;

    for(std::string::iterator itr_ref = reference.begin()+options.k+1; itr_ref != reference.end(); itr_ref++){ // iterate over reference
      hashId=hash.bigger();
      this->lookup(hashId,coverage);
      if(coverage!=0){
        coverages.push_back(coverage);
      }
      hash.roll(*itr_ref);
    }

    // std::cerr << "\nCoverage distribution:\n";
    // for(int i=0;i<coverages.size();i++){
    //   std::cerr << coverages[i] << ", ";
    // }

    // partial sort
    std::nth_element(coverages.begin(), coverages.begin()+coverages.size()/2, coverages.end());
    // get median
    uint32_t median = coverages[coverages.size()/2];
    // std::cerr << "\n\nMedian coverage: " << median << "\n\n";
    return median;
  }

};

int buildCoverageIndex(coverageIndex_t & coverageIndex, std::string & readfile1_name, std::string & readfile2_name, uint32_t k, bool & paired);
void loadReads(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, std::string & readfile1_name, std::string & readfile2_name, bool & paired);
void fillCoverageIndex(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex_t & coverageIndex, uint32_t k, bool & paired);
void fillIndexFromReads(StringSet<Dna5String> & reads, coverageIndex_t & coverageIndex, uint32_t k);
#endif
