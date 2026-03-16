#ifndef CALL_INDEX_H
#define CALL_INDEX_H

#include <vector>
#include <cmath>
#include <stdint.h>
#include <string>
#include <iostream>
# include "seqan/sequence.h"

// #include "classes.h"
using namespace seqan;

struct hashId_t{
  uint64_t i1;
  uint64_t i2;
  uint64_t i3;

  hashId_t(){
    i1=0;
    i2=0;
    i3=0;
  }

  bool operator==(const hashId_t& rhs){
    return ( i1 == rhs.i1 && i2 == rhs.i2 && i3 == rhs.i3 );
  }

  bool operator>(const hashId_t& rhs){
    return ( i1>rhs.i1 || i1==rhs.i1 && i2>rhs.i2 || i1==rhs.i1 && i2==rhs.i2 && i3>rhs.i3);
  }
  //
  // bool operator<(const hashId_t& rhs){
  //   return ( i1<rhs.i1 || i1==rhs.i1 && i2<rhs.i2 || i1==rhs.i1 && i2==rhs.i2 && i3<rhs.i3);
  // }

};

inline bool operator<(const hashId_t& lhs, const hashId_t& rhs){
  return ( lhs.i1<rhs.i1 || lhs.i1==rhs.i1 && lhs.i2<rhs.i2 || lhs.i1==rhs.i1 && lhs.i2==rhs.i2 && lhs.i3<rhs.i3);
}

struct hash_t{ // hash kmers for 64 < k <= 96
  // binary represantation of kmer
  hashId_t fwd;
  // binary represantation of kmers reverse complement
  hashId_t rev;

  uint32_t k;
  uint64_t i1MaxVal; // maximum value for i1 depending on k
  uint64_t i1Shift;  // bits to shift to get i1s first nucleotide
  uint64_t i2MaxVal; // maximum value for i2 depending on k
  uint64_t i2Shift; // bits to shift to get i2s first nucleotide
  uint64_t i3MaxVal; // maximum value for i3 depending on k
  uint64_t i3Shift;  // bits to shift to get i3s first nucleotide
  typedef  void(hash_t::*hash_tMemFn)(char x);  // Please do this!
  hash_tMemFn roll_p; // function pointer to correct rolling function

  hash_t(){}

  hash_t(Dna5String & kmer){
    k = length(kmer);
    if(k<32){
      roll_p = &hash_t::roll1;
      i1MaxVal = 1;
      for(int i = 1; i<=k; i++){i1MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i1MaxVal -= 1;      i2MaxVal = 0;
      i3MaxVal = 0;
    }else if(k<64){
      roll_p = &hash_t::roll2;
      i1MaxVal = 18446744073709551615U; //(2^64-1)
      i2MaxVal = 1;
      for(int i = 1; i<=k-32; i++){i2MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i2MaxVal -= 1;
      i3MaxVal = 0;
    }else{
      roll_p = &hash_t::roll3;
      i1MaxVal = 18446744073709551615U;
      i2MaxVal = 18446744073709551615U;
      i3MaxVal = 1;
      for(int i = 1; i<=k-64; i++){i3MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i3MaxVal -= 1;
    }

    i1Shift = (k-1)*2-1;
    i2Shift = (k-33)*2;
    i3Shift = (k-65)*2;
    fwd.i1 = 0; fwd.i2 = 0; fwd.i3 = 0;
    rev.i1 = 0; rev.i2 = 0; rev.i3 = 0;

    // calculate kmer representation
    uint32_t kmax=32*(k>=32)+k*(k<32);
    for(Iterator<Dna5String>::Type it = begin(kmer); it < begin(kmer)+kmax; it++){
      fwd.i1 = (fwd.i1 << 2) | (((char)(*it) >> 1) & 3);
    }
    if(k>32){
      kmax=64*(k>=64)+k*(k<64);
      for(Iterator<Dna5String>::Type it = begin(kmer)+32; it < begin(kmer)+kmax; it++){
        fwd.i2 = (fwd.i2 << 2) | (((char)(*it) >> 1) & 3);
      }
      if(k>64){
        for(Iterator<Dna5String>::Type it = begin(kmer)+64; it < end(kmer); it++){
          fwd.i3 = (fwd.i3 << 2) | (((char)(*it) >> 1) & 3);
        }
      }
    }

    // calculate reverse complement representation
    kmax=33*(k>=32)+(k+1)*(k<32);
    for(Iterator<Dna5String>::Type it = end(kmer)-1; it != end(kmer)-kmax; it--){
      rev.i1 = (rev.i1 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
    }
    if(k>32){
      kmax=65*(k>=64)+(k+1)*(k<64);
      for(Iterator<Dna5String>::Type it = end(kmer)-33; it != end(kmer)-65; it--){
        rev.i2 = (rev.i2 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
      }
      if(k>64){
        for(Iterator<Dna5String>::Type it = end(kmer)-65; it != end(kmer)-(k+1); it--){
          rev.i3 = (rev.i3 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
        }
      }
    }
  }

  hash_t(std::string kmer){
    k = kmer.size();
    if(k<32){
      roll_p = &hash_t::roll1;
      i1MaxVal = 1;
      for(int i = 1; i<=k; i++){i1MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i1MaxVal -= 1;
      i2MaxVal = 0;
      i3MaxVal = 0;
    }else if(k<64){
      roll_p = &hash_t::roll2;
      i1MaxVal = 18446744073709551615U; //(2^64-1)
      i2MaxVal = 1;
      for(int i = 1; i<=k-32; i++){i2MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i2MaxVal -= 1;
      i2MaxVal = std::pow(2,(k-32)*2)-1;
      i3MaxVal = 0;
    }else{
      roll_p = &hash_t::roll3;
      i1MaxVal = 18446744073709551615U;
      i2MaxVal = 18446744073709551615U;
      i3MaxVal = 1;
      for(int i = 1; i<=k-64; i++){i3MaxVal *= 4;}  // std::pow gibt für k>26 das falsche ergebnis bei 4**k
      i3MaxVal -= 1;
    }

    i1Shift = (k-1)*2-1;
    i2Shift = (k-33)*2;
    i3Shift = (k-65)*2;
    fwd.i1 = 0; fwd.i2 = 0; fwd.i3 = 0;
    rev.i1 = 0; rev.i2 = 0; rev.i3 = 0;
    // calculate kmer representation
    uint32_t kmax=32*(k>=32)+k*(k<32);
    for(std::string::iterator it = kmer.begin(); it < kmer.begin()+kmax; it++){
      fwd.i1 = (fwd.i1 << 2) | (((*it) >> 1) & 3);
    }
    if(k>32){
      kmax=64*(k>=64)+k*(k<64);
      for(std::string::iterator it = kmer.begin()+32; it < kmer.begin()+kmax; it++){
        fwd.i2 = (fwd.i2 << 2) | (((*it) >> 1) & 3);
      }
      if(k>64){
        for(std::string::iterator it = kmer.begin()+64; it < kmer.end(); it++){
          fwd.i3 = (fwd.i3 << 2) | (((*it) >> 1) & 3);
        }
      }
    }
    // calculate reverse complement representation
    kmax=32*(k>=32)+k*(k<32);
    for(std::string::reverse_iterator it = kmer.rbegin(); it != kmer.rbegin()+kmax; it++){
      rev.i1 = (rev.i1 << 2) | ((((*it) >> 1) & 3) ^ 2);
    }
    if(k>32){
      kmax=64*(k>=64)+k*(k<64);
      for(std::string::reverse_iterator it = kmer.rbegin()+32; it != kmer.rbegin()+kmax; it++){
        rev.i2 = (rev.i2 << 2) | ((((*it) >> 1) & 3) ^ 2);
      }
      if(k>64){
        for(std::string::reverse_iterator it = kmer.rbegin()+64; it != kmer.rbegin()+k; it++){
          rev.i3 = (rev.i3 << 2) | ((((*it) >> 1) & 3) ^ 2);
        }
      }
    }
  }

  void rehash(Dna5String & kmer){
    // k = length(kmer);
    // i3MaxVal = std::pow(2,(k-64)*2)-1;
    // i3Shift = (k-65)*2;
    fwd.i1 = 0; fwd.i2 = 0; fwd.i3 = 0;
    rev.i1 = 0; rev.i2 = 0; rev.i3 = 0;
    // calculate kmer representation
    uint32_t kmax=32*(k>=32)+k*(k<32);
    for(Iterator<Dna5String>::Type it = begin(kmer); it < begin(kmer)+kmax; it++){
      fwd.i1 = (fwd.i1 << 2) | (((char)(*it) >> 1) & 3);
    }
    if(k>32){
      kmax=64*(k>=64)+k*(k<64);
      for(Iterator<Dna5String>::Type it = begin(kmer)+32; it < begin(kmer)+kmax; it++){
        fwd.i2 = (fwd.i2 << 2) | (((char)(*it) >> 1) & 3);
      }
      if(k>64){
        for(Iterator<Dna5String>::Type it = begin(kmer)+64; it < end(kmer); it++){
          fwd.i3 = (fwd.i3 << 2) | (((char)(*it) >> 1) & 3);
        }
      }
    }

    // calculate reverse complement representation
    kmax=33*(k>=32)+(k+1)*(k<32);
    for(Iterator<Dna5String>::Type it = end(kmer)-1; it != end(kmer)-kmax; it--){
      rev.i1 = (rev.i1 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
    }
    if(k>32){
      kmax=65*(k>=64)+(k+1)*(k<64);
      for(Iterator<Dna5String>::Type it = end(kmer)-33; it != end(kmer)-65; it--){
        rev.i2 = (rev.i2 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
      }
      if(k>64){
        for(Iterator<Dna5String>::Type it = end(kmer)-65; it != end(kmer)-(k+1); it--){
          rev.i3 = (rev.i3 << 2) | ((((char)(*it) >> 1) & 3) ^ 2);
        }
      }
    }
    return;
  }

  void roll(char newNuc){
    (*this.*roll_p)(newNuc);
  }

  void roll1(char newNuc){
    fwd.i1 = ((fwd.i1 << 2) & i1MaxVal) | ((newNuc >> 1) & 3);

    rev.i1 = ((rev.i1 >> 2) & i1MaxVal) | ((uint64_t)((newNuc & 6) ^ 4) << i1Shift);
    return;
  }

  void roll2(char newNuc){
    fwd.i1 = (fwd.i1 << 2)| (fwd.i2 >> i2Shift);
    fwd.i2 = ((fwd.i2 << 2) & i2MaxVal) | ((newNuc >> 1) & 3);

    rev.i2 = ((rev.i2 >> 2) & i2MaxVal) | (rev.i1 & 3) << i2Shift;
    rev.i1 = (rev.i1 >> 2) | ((uint64_t)((newNuc & 6) ^ 4) << 61);
    return;
  }

  void roll3(char newNuc){
    fwd.i1 = (fwd.i1 << 2) | (fwd.i2 >> 62);
    fwd.i2 = (fwd.i2 << 2) | (fwd.i3 >> i3Shift);
    fwd.i3 = ((fwd.i3 << 2) & i3MaxVal) | ((newNuc >> 1) & 3);

    rev.i3 = ((rev.i3 >> 2) & i3MaxVal) | ((rev.i2 & 3) << i3Shift);
    rev.i2 = (rev.i2 >> 2) | (rev.i1 & 3) << 62;
    rev.i1 = (rev.i1 >> 2) | ((uint64_t)((newNuc & 6) ^ 4) << 61);
    return;
  }

  hashId_t hamming(uint64_t pos, uint64_t alt){ // alt is element of {1,2,3} pos is element of {0, ... , 31}
    hashId_t fwd_buff =fwd;
    if(pos < 32){
      fwd_buff.i1 = fwd.i1 ^ (alt << (pos*2));
    }
    else if(pos < 64){
      fwd_buff.i2 = fwd.i2 ^ (alt << ((pos-32)*2));
    }
    else{
      fwd_buff.i3 = fwd.i3 ^ (alt << ((pos-64)*2));
    }
    // uint64_t fwd_buff = fwd ^ (alt << (pos*2));
    hashId_t rev_buff =rev;
    uint32_t k_2;
    if(pos < 32){
      k_2=k*(k<31)+31*(k>=31);
      rev_buff.i3 = rev.i3 ^ (alt << (k_2*2-2-(pos*2)));
    }
    else if(pos < 64){
      k_2=(k-32)*(k<31)+31*(k>=64);
      rev_buff.i2 = rev.i2 ^ (alt << (k_2*2-2-((pos-32)*2)));
    }
    else{
      k_2=k-64;
      rev_buff.i1 = rev.i1 ^ (alt << (k_2*2-2-((pos-64)*2)));
    }
    // uint64_t rev_buff = rev ^ (alt << (k*2-2-(pos*2)));

    if(fwd_buff > rev_buff){
      return fwd_buff;
    }
    return rev_buff;
  }

  hashId_t rollInPlace(char newNuc, hash_t & buff, bool & orientation){ // rolls from position of "this" but uses dummy to store intermediate results without constant redeclaration
    buff.fwd.i1 = (fwd.i1 << 2) | (fwd.i2 >> 62);
    buff.fwd.i2 = (fwd.i2 << 2) | (fwd.i3 >> i3Shift);
    buff.fwd.i3 = ((fwd.i3 << 2) & i3MaxVal) | ((newNuc >> 1) & 3);

    buff.rev.i3 = ((rev.i3 >> 2) & i3MaxVal) | ((rev.i2 & 3) << i3Shift);
    buff.rev.i2 = (rev.i2 >> 2) | (rev.i1 & 3) << 62;
    buff.rev.i1 = (rev.i1 >> 2) | ((uint64_t)((newNuc & 6) ^ 4) << 61);
    return buff.bigger(orientation);
  }

  hashId_t rollBackInPlace(char newNuc, hash_t & buff, bool & orientation){ // rolls backwards from position of "this" but uses dummy to store intermediate results without constant redeclaration
    buff.fwd.i3 = ((fwd.i3 >> 2) & i3MaxVal) | ((fwd.i2 & 3) << i3Shift);
    buff.fwd.i2 = (fwd.i2 >> 2) | (fwd.i1 & 3) << 62;
    buff.fwd.i1 = (fwd.i1 >> 2) | ((uint64_t)(newNuc & 6) << 61);

    buff.rev.i1 = (rev.i1 << 2) | (rev.i2 >> 62);
    buff.rev.i2 = (rev.i2 << 2) | (rev.i3 >> i3Shift);
    buff.rev.i3 = ((rev.i3 << 2) & i3MaxVal) | (((newNuc ^ 4) >> 1) & 3);
    return buff.bigger(orientation);
  }

  hashId_t bigger(){
    if(fwd > rev){
      return fwd;
    }
    return rev;
  }

  hashId_t bigger(bool & orientation){
    if(fwd > rev){
      orientation = false;
      return fwd;
    }
    orientation = true;
    return rev;
  }

  void print(std::string name){
    std::cerr <<  name << ":\ti1: " << fwd.i1 << " i2: " << fwd.i2 << " i3: " << fwd.i3 << "\n";
    return;
  }

  void printr(std::string name){
    std::cerr <<  name << ":\tr1: " << rev.i1 << " r2: " << rev.i2 << " r3: " << rev.i3 << "\n";
    return;
  }

  void printKmer(std::string name){
    std::string kmer="";
    uint64_t bin;
    uint32_t lim;
    if(k<=32){lim=k-1;}else{lim=31;}
    for(int i=lim; i>=0; i--){
      bin = fwd.i1>>(i*2) & 3;
      if(bin==0){kmer+="A";}
      if(bin==1){kmer+="C";}
      if(bin==2){kmer+="T";}
      if(bin==3){kmer+="G";}
    }
    if(k<=64){lim=k-33;}else{lim=31;}
    for(int i=lim; i>=0; i--){
      bin = fwd.i2>>(i*2) & 3;
      if(bin==0){kmer+="A";}
      if(bin==1){kmer+="C";}
      if(bin==2){kmer+="T";}
      if(bin==3){kmer+="G";}
    }
    for(int i=k-65; i>=0; i--){
      bin = fwd.i3>>(i*2) & 3;
      if(bin==0){kmer+="A";}
      if(bin==1){kmer+="C";}
      if(bin==2){kmer+="T";}
      if(bin==3){kmer+="G";}
    }
    std::cerr << name << ":\t" << kmer << "\n";
  }
};

struct hashTableField_t{
  hashId_t id;
  std::vector<std::string> unitigs;  // unitigs in which the kmer appears
  std::vector<uint32_t> positions; // position in unitigs
  std::vector<bool> orientations;

  hashTableField_t(){}
};

struct hashTable_t{
  std::vector<hashTableField_t> table;
  uint64_t tableSize;

  hashTable_t(){
  }

  void resize(uint64_t tableSize_p){
    tableSize = tableSize_p;
    table.resize(tableSize);
    return;
  }

  bool getpos(hashId_t & hashId, uint64_t & pos){ // look up hashId in index. modifies pos and returns whether kmer is present
    pos = hashId.i2 % tableSize;
    for(int d=0; d<500; d++){
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

  bool insert(hashId_t & hashId, std::string & unitig, uint32_t & position, bool orientation){ // insert kmer into hashtable or expand list of unitigs for kmer
    uint64_t pos;

    pos = hashId.i2 % tableSize;
    for(int d=0; d<500; d++){
      // first try in loop
      table[pos].id = (table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        table[pos].unitigs.push_back(unitig);
        table[pos].positions.push_back(position);
        table[pos].orientations.push_back(orientation);
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId.i1 >> (d & 31)) )% tableSize; // first shift in loop
      // second try in loop
      table[pos].id = (table[pos].id.i1==0 && table[pos].id.i2==0 && table[pos].id.i3==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        table[pos].unitigs.push_back(unitig);
        table[pos].positions.push_back(position);
        table[pos].orientations.push_back(orientation);
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId.i3 << (d & 31))) % tableSize; // second shift in loop
    }
    return false;
  }

  bool lookup(hashId_t & hashId, std::vector<std::string> & unitigs, std::vector<uint32_t> & positions, std::vector<bool> & orientations){
    uint64_t pos;
    if(this->getpos(hashId, pos)){
      unitigs = table[pos].unitigs;
      positions = table[pos].positions;
      orientations = table[pos].orientations;

      return true;
    }
    return false;
  }


};

#endif
