# include <iostream>
# include <fstream>
# include <vector>
# include <random>
# include "molephase.h"
# include "parse.h"
# include "coverageIndex.h"
# include "seqan/seq_io.h"
# include "seqan/sequence.h"
# include <bitset>
# include "call_index.h"

using namespace seqan;

void loadReadData(ReadData & readData, phaseOptions & options);
void fillIndex32FromReads(StringSet<Dna5String> & reads, coverageIndex32_t & coverageIndex, uint32_t k);
void fillCoverageIndex32(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex32_t & coverageIndex, uint32_t k);
bool loadReference(Dna5String & reference, phaseOptions & options);
// void selectHeterozygousKmers(Dna5String & reference, coverageIndex32_t & coverageIndex, std::vector<uint64_t> & refKmers, std::vector<uint64_t> & altKmers, std::set<uint64_t> & hetKmers, phaseOptions & options);
void selectHeterozygousKmers(Dna5String & reference, coverageIndex_t & coverageIndex, std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::set<hashId_t> & hetKmers, std::vector<uint32_t> & kmersPos, phaseOptions & options); //[TEST2]
void buildKmerToBarcodes(std::set<hashId_t> & hetKmers, ReadData & readData, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, phaseOptions & options);
void reduceKmerSet(std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::vector<uint32_t> & kmersPos, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, phaseOptions & options); //[TEST2]
void initializeHaplotypes(std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::vector<hashId_t> & hap1, std::vector<hashId_t> & hap2, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport);
void countBarcodeSupport(std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::vector<hashId_t> & hap1, std::vector<hashId_t> & hap2);
void phaseKmers(std::vector<hashId_t> & hap1, std::vector<hashId_t> & hap2, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, phaseOptions & options);
void phaseBarcodes(std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2);
void touchOutput(phaseOptions & options);
void writeToUnphased(ReadData & readData, phaseOptions & options);
void writePhasedOutput(std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2, ReadData & readData, phaseOptions & options);
void writePhasedOutput_paired(std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2, ReadData & readData, phaseOptions & options);

void countConnectedComponents(std::map<std::string, std::set<uint64_t>> & barcodeToKmers, std::set<hashId_t> & hetKmers);

int main(int argc, char const **argv){

  // Argument parser #############################################################################
  phaseOptions options;
  seqan::ArgumentParser::ParseResult res = parsePhaseCommandLine(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) return res;

  if(options.end <= options.start){ std::cerr << "\nError: Parameter start (option -s ["<< options.start <<"]) must hold smaller value than parameter end (option -e [" << options.end << "])\n"; return 1;}

  printPhaseParseResult(options, res);

  // load readfiles
  std::cerr << "Loading readfiles...";

  ReadData readData;
  loadReadData(readData, options);

  if(seqan::empty(readData.reads1) && seqan::empty(readData.reads2)){ // if no reads in file: Write empty output.
    touchOutput(options);
    std::cerr << "\nReadfiles are empty. Writing empty output files.\n";
    return 0;
  }

  std::cerr << "............done.\n";

  std::cerr << "Building coverage index...";

  coverageIndex_t coverageIndex;
  readData.getLength();
  uint64_t index_size = readData.seqLength-((length(readData.reads1)+length(readData.reads2))*(options.k+1));
  index_size *= 1.2;
  coverageIndex.resize(index_size);
  fillIndexFromReads(readData.reads1, coverageIndex, options.k);
  if(options.paired){
    fillIndexFromReads(readData.reads2, coverageIndex, options.k);
  }

  std::cerr << "......done.\n";

  std::cerr << "Loading reference.........";

  // load reference sequence
  Dna5String reference;
  if(!loadReference(reference, options)){
    std::cerr << "\nERROR: Failed to load reference, shuting down!\nMake sure that the provided reference file (-r option) is valid and accessible.\n";
    return 1;
  }

  std::cerr << "......done.\n";

  std::cerr << "Phasing barcodes...";
  // select heterozygous kmers (ref/non-ref kmer pairs with hamming-dist=1 and similar coverage)
  std::set<hashId_t> hetKmers;
  std::vector<hashId_t> refKmers;
  std::vector<hashId_t> altKmers;
  std::vector<uint32_t> kmersPos; //[TEST2]
  refKmers.reserve(10000);
  altKmers.reserve(10000);
  kmersPos.reserve(10000); //[TEST2]
  selectHeterozygousKmers(reference, coverageIndex, refKmers, altKmers, hetKmers, kmersPos, options); //[TEST2]

  if(refKmers.empty()){ // if no reads in file: Write empty output.
    writeToUnphased(readData, options);
    std::cerr << "\nNo phasable kmers found. Writing all output to unphased.\n";
    return 0;
  }

  std::map<hashId_t, std::set<std::string>> kmerToBarcodes;
  std::map<std::string, std::pair<uint32_t, uint32_t>> barcodeSupport;
  buildKmerToBarcodes(hetKmers, readData, kmerToBarcodes, barcodeSupport, options); // build index that maps hetkmers to barcodes

  reduceKmerSet(refKmers, altKmers, kmersPos, kmerToBarcodes, options); //[TEST2]

  std::vector<hashId_t> hap1;
  std::vector<hashId_t> hap2;
  hap1.reserve(refKmers.size());
  hap2.reserve(refKmers.size());
  initializeHaplotypes(refKmers, altKmers, kmerToBarcodes, hap1, hap2, barcodeSupport);

  phaseKmers(hap1, hap2, kmerToBarcodes, barcodeSupport, options);

  // TODO: potentially adjust to add barcodes with comparable support for both haplotypes  to unphased barcodes
  std::set<std::string> barcodesHap1;
  std::set<std::string> barcodesHap2;
  phaseBarcodes(barcodeSupport, barcodesHap1, barcodesHap2); // add barcode to haplotype with higher barcodeSupport

  std::cerr << ".............done.\n";
  std::cerr << "Writing output...";

  if(options.paired){
    writePhasedOutput_paired(barcodesHap1, barcodesHap2, readData, options);
  }else{
    writePhasedOutput(barcodesHap1, barcodesHap2, readData, options);
  }

  std::cerr << "...............done.\n";

  return 0;
}

void fillIndex32FromReads(StringSet<Dna5String> & reads, coverageIndex32_t & coverageIndex, uint32_t k){
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  Iterator<Dna5String>::Type itread;
  uint64_t hashId;

  for (TStringSetIterator itreads = begin(reads); itreads!=end(reads); ++itreads){ //  for every read in reads

    // add first kmer to index
    Dna5String kmer = prefix(*itreads,k);
    hash32_t hash(kmer);
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

void fillCoverageIndex32(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex32_t & coverageIndex, uint32_t k){
  fillIndex32FromReads(reads1, coverageIndex, k);
  fillIndex32FromReads(reads2, coverageIndex, k);
  return;
}

// IO functions
bool loadReference(Dna5String & reference, phaseOptions & options){ // load ref.fa.fai and extract region of interest from ref.fa
  // load ref.fa.fai or create it if it doesnt exist.
  // std::cerr << "Loading ref.fai...";
  FaiIndex faiIndex;
  if (!open(faiIndex, toCString(options.reference))){
    std::cerr << "..........................failed.\nBuilding ref.fai..";
    if (!build(faiIndex, toCString(options.reference))){
      std::cerr << "\nERROR: FASTA index could not be loaded or built.\n";
      return 1;
    }
    if (!save(faiIndex)) // Name is stored from when reading.
    {
      std::cerr << "\nWARNING: FASTA index could not be written to disk.\n";
    }
  }else{
    // std::cerr << "............................done.\n";
  }

  // look up chromosome position in fai
  unsigned idx = 0;
  if (!getIdByName(idx, faiIndex, options.chrom)){
    std::cout << "\n\n\tERROR: FAI index has no entry for:\t" << options.chrom << "\n";
    return false;
  }

  // extract region of interest
  IupacString reference_IupacString;
  readRegion(reference_IupacString, faiIndex, idx, options.start-1, options.end);
  reference = reference_IupacString;

  return true;
}

uint_fast8_t getBarcodeLength(std::string id){
  std::string barcode;
  barcode=id.substr(id.find("BX:Z:")+5,id.find(' ',id.find("BX:Z:")+5)-id.find("BX:Z:")-5);              // determine BX:Z: entry
  return barcode.size();
}

std::string getBarcode(std::string id1, uint_fast8_t barcode_length){ //retreive the barcode from linked reads
  std::size_t pos=id1.find(" ");
  std::string barcode;
  if(pos<1000){
    id1=id1.substr(id1.find(" "),10000);
    barcode=id1.substr(id1.find("BX:Z:")+5,barcode_length);
  }else{
    barcode="BAD_BARCODE_____";
  }
  return barcode;
}

std::string getReadName(std::string id1, uint_fast8_t barcode_length){ //retreive the readName (for long reads)
  return id1.substr(0, id1.find(" "));
}

std::string getReadName(std::string id1){ //retreive the readName (for long reads)
  return id1.substr(0, id1.find(" "));
}

// Functions for Barcode Phasing
bool covRelation(uint32_t & ref_coverage, uint32_t & alt_coverage){
  if(alt_coverage < 2){ // best: 1 tried: 2,3
    return false;
  }
  if((float)(alt_coverage/ref_coverage) > 1.5 || (float)(ref_coverage/alt_coverage) > 1.5){ // best:1.5 tried: 1.2, 2, 3
    return false;
  }
  // if(alt_coverage > 10 || ref_coverage > 10){
  //   return false;
  // }
  // std::cerr << ref_coverage << " : " << alt_coverage << "\n";
  return true; // coverage relation fits
}

// void selectHeterozygousKmers(Dna5String & reference, coverageIndex32_t & coverageIndex, std::vector<uint64_t> & refKmers, std::vector<uint64_t> & altKmers, std::set<uint64_t> & hetKmers, phaseOptions & options){
void selectHeterozygousKmers(Dna5String & reference, coverageIndex_t & coverageIndex, std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::set<hashId_t> & hetKmers, std::vector<uint32_t> & kmersPos, phaseOptions & options){ //[TEST2]

  uint32_t ref_coverage;
  uint32_t alt_coverage;

  // initialize kmer lookup
  Dna5String kmer = infix(reference, 0, options.k);
  hash_t hash(kmer);
  hashId_t hashId; // hash value of reference kmer

  // iterate over reference kmers
  for(int i = options.k; i<length(reference); i++){ // iterate over reference
    // lookup coverage of kmer in index
    hashId=hash.bigger();
    coverageIndex.lookup(hashId,ref_coverage);

    if(ref_coverage > 1){
      // check for kmers with hamming distance 1 to hetKmer
      bool candidateFound = false;
      bool failed = false;
      hashId_t hashIdBuff;
      hashId_t hashIdCandidate; // hash value of potential alt kmer
      // define hashBuff as every possible altered version of hashId
      for(int pos = 0; pos!=options.k; pos++){ // for every position in kmer
        for(int alt = 1; alt!=4; alt++){
          hashIdBuff=hash.hamming(pos,alt);
          coverageIndex.lookup(hashIdBuff,alt_coverage);
          if(covRelation(ref_coverage, alt_coverage)){ // if coverage relation fits: mark as candidate
            if(candidateFound == true){failed=true; break;} // if more than 1 candidate for this position is found
            hashIdCandidate = hashIdBuff;
            candidateFound = true;
          }
        }
        if(failed==true){break;}
      }
      if(failed==false && candidateFound==true){
        // insert kmer to "random" haplotype
        hetKmers.insert(hashId);
        hetKmers.insert(hashIdCandidate);
        refKmers.push_back(hashId);
        altKmers.push_back(hashIdCandidate);
        kmersPos.push_back(i); //[TEST2]
      }
    } // if refKmer valid
    hash.roll(reference[i]); // roll to next refKmer
  } // for kmer in reference
  return;
}

void fillKmerToBarcodes(std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::set<hashId_t> & hetKmers, hash_t & hash, std::string & identifier){
  hashId_t hashId = hash.bigger();
  if(hetKmers.count(hashId)){       // if kmer (hashId) in hetKmers: update maps
    if(kmerToBarcodes.count(hashId)){ // add barcode to barcode list  of kmer
      (kmerToBarcodes[hashId]).insert(identifier);
    }else{ // create entry for kmer with barcode as first element
      kmerToBarcodes[hashId]={identifier};
    }
  }
  return;
}

void buildKmerToBarcodes(std::set<hashId_t> & hetKmers, ReadData & readData, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, phaseOptions & options){
  typedef Iterator<StringSet<Dna5String> >::Type ReadsIterator;
  typedef Iterator<StringSet<CharString> >::Type IdsIterator;
  typedef Iterator<Dna5String>::Type ReadIterator;
  ReadsIterator itReads2 = begin(readData.reads2);
  IdsIterator itrIds = begin(readData.ids1);

  uint_fast8_t barcode_length;
  std::string (*getIdentifier)(std::string id1, uint_fast8_t barcode_length); // function pointer
  if(options.paired){
    barcode_length = getBarcodeLength(toCString(*itrIds));
    getIdentifier = &getBarcode;
  }else{
    barcode_length = 0;
    getIdentifier = &getReadName;
  }

  std::string identifier; // readname for long reads/ barcode for linked reads

  // iterate over reads
  for (ReadsIterator itReads1 = begin(readData.reads1); itReads1!=end(readData.reads1); ++itReads1){ //  for every readpair in reads1 and reads2
    // extract barcode from id
    identifier = getIdentifier(toCString(*itrIds), barcode_length);
    itrIds++;
    barcodeSupport[identifier]=std::make_pair(0,0);
    // initialize hash for read1
    Dna5String kmer = prefix(*itReads1,options.k);
    hash_t hash(kmer);
    fillKmerToBarcodes(kmerToBarcodes, hetKmers, hash, identifier);

    // iterate over kmers in read 1
    for(ReadIterator itRead=begin(*itReads1)+options.k;itRead!=end(*itReads1);itRead++){
      hash.roll(*itRead);
      fillKmerToBarcodes(kmerToBarcodes, hetKmers, hash, identifier);
    }

    if(!options.paired){continue;}
    // reinitialize hash for read2
    Dna5String kmer2 = prefix(*itReads2,options.k);
    hash_t hash2(kmer2);
    fillKmerToBarcodes(kmerToBarcodes, hetKmers, hash2, identifier);

    // iterate over kmers in read 2
    for(ReadIterator itRead=begin(*itReads2)+options.k;itRead!=end(*itReads2);itRead++){
      hash2.roll(*itRead);
      fillKmerToBarcodes(kmerToBarcodes, hetKmers, hash2, identifier);
    }

    itReads2++;
  }
  return;
}

void loadReadData(ReadData & readData, phaseOptions & options){
  SeqFileIn file1(toCString(options.readfile1));
  readRecords(readData.ids1, readData.reads1, readData.quals1, file1);
  close(file1);
  if(options.paired){
    SeqFileIn file2(toCString(options.readfile2));
    readRecords(readData.ids2, readData.reads2, readData.quals2, file2);
    close(file2);
  }
  return;
}

void reduceKmerSet(std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::vector<uint32_t> & kmersPos, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, phaseOptions & options){ //[TEST2]
  // check for duplicate Kmers
  std::set<hashId_t> seenKmers;
  std::set<hashId_t> duplicateKmers;
  for(std::vector<hashId_t>::iterator it = refKmers.begin(); it!= refKmers.end(); it++){
    if(seenKmers.count(*it)){ // if kmer already seen
      duplicateKmers.insert(*it);
    }else{
      seenKmers.insert(*it);
    }
  }
  for(std::vector<hashId_t>::iterator it = altKmers.begin(); it!= altKmers.end(); it++){
    if(seenKmers.count(*it)){ // if kmer already seen
      duplicateKmers.insert(*it);
    }else{
      seenKmers.insert(*it);
    }
  }

  // std::cerr << "\nduplicteKmers.size(): " << duplicateKmers.size() << "\n";
  // std::cerr << "before:\naltKmers.size(): " << altKmers.size() << "\t refKmers.size(): " << refKmers.size() << "\n";

  // remove duplicate Kmers
  if(duplicateKmers.size()){
    std::vector<hashId_t> refKmersBuff; refKmersBuff.reserve(refKmers.size());
    std::vector<hashId_t> altKmersBuff; altKmersBuff.reserve(altKmers.size());
    std::vector<uint32_t> kmersPosBuff; kmersPosBuff.reserve(kmersPos.size());
    std::vector<hashId_t>::iterator itAlt = altKmers.begin();
    std::vector<uint32_t>::iterator itPos = kmersPos.begin();
    for(std::vector<hashId_t>::iterator itRef = refKmers.begin(); itRef!=refKmers.end(); itRef++){
      if( !(duplicateKmers.count(*itRef) || duplicateKmers.count(*itAlt)) ){
        refKmersBuff.push_back(*itRef);
        altKmersBuff.push_back(*itAlt);
        kmersPosBuff.push_back(*itPos);
      }
      itAlt++; itPos++;
    }
    refKmers=refKmersBuff;
    altKmers=altKmersBuff;
    kmersPos=kmersPosBuff;
  }
  // std::cerr << "after:\naltKmers.size(): " << altKmers.size() << "\t refKmers.size(): " << refKmers.size() << "\n";

  // merge adjacent Kmers
  std::vector<hashId_t> refKmersBuff; refKmersBuff.reserve(refKmers.size());
  std::vector<hashId_t> altKmersBuff; altKmersBuff.reserve(altKmers.size());
  std::vector<hashId_t>::iterator itRef = refKmers.begin();
  std::vector<hashId_t>::iterator itAlt = altKmers.begin();
  std::vector<uint32_t>::iterator itPos = kmersPos.begin();
  hashId_t refKmer = *itRef; hashId_t altKmer=*itAlt; uint32_t kmerPos=*itPos;
  refKmersBuff.push_back(*itRef);
  altKmersBuff.push_back(*itAlt);
  itRef++; itAlt++; itPos++;
  for(itRef; itRef!=refKmers.end(); itRef++){
    if(*itPos - kmerPos < options.k){ // if new kmer within x bases of old kmer
      //merge kmers (update kmerToBarcodes)
      kmerToBarcodes[refKmer].insert(kmerToBarcodes[*itRef].begin(),kmerToBarcodes[*itRef].end());
      kmerToBarcodes[altKmer].insert(kmerToBarcodes[*itAlt].begin(),kmerToBarcodes[*itAlt].end());
    }else{
      // add new kmer
      refKmer=*itRef;
      altKmer=*itAlt;
      kmerPos=*itPos;
      refKmersBuff.push_back(*itRef);
      altKmersBuff.push_back(*itAlt);
      // std::cerr << "kmerPos: " <<s *itPos << "\n";
    }
    itAlt++; itPos++;
  }
  refKmers=refKmersBuff;
  altKmers=altKmersBuff;

  // std::cerr << "merged:\naltKmers.size(): " << altKmers.size() << "\t refKmers.size(): " << refKmers.size() << "\n";


  return;
}

int32_t getPhaseScore(hashId_t & kmer1, hashId_t & kmer2, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::map<hashId_t,std::set<std::string>> & kmerToBarcodes){
  int32_t score = 0;
  for(std::set<std::string>::iterator itr_bcSet=kmerToBarcodes[kmer1].begin(); itr_bcSet!=kmerToBarcodes[kmer1].end(); itr_bcSet++){ // for every barcode of kmer 1
    // sum up the support of that barcode within the current haplotype configuration
    score+=barcodeSupport[*itr_bcSet].first;
    score-=barcodeSupport[*itr_bcSet].second;
  }
  for(std::set<std::string>::iterator itr_bcSet=kmerToBarcodes[kmer2].begin(); itr_bcSet!=kmerToBarcodes[kmer2].end(); itr_bcSet++){ // for every barcode that supports hap2
    // sum up the support of that barcode within the current haplotype configuration
    score-=barcodeSupport[*itr_bcSet].first;
    score+=barcodeSupport[*itr_bcSet].second;
  }
  return score;
}

void adjustBarcodeSupport(std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, hashId_t kmer1, hashId_t kmer2){
  // update barcode support map for the swap of haplotype pair [kmer1, kmer2]
  std::set<std::string>::iterator itr_bcSet;
  for(itr_bcSet=kmerToBarcodes[kmer1].begin(); itr_bcSet!=kmerToBarcodes[kmer1].end(); itr_bcSet++){ // for every barcode that supports hap1
    barcodeSupport[*itr_bcSet].first--;
    barcodeSupport[*itr_bcSet].second++;
  }
  for(itr_bcSet=kmerToBarcodes[kmer2].begin(); itr_bcSet!=kmerToBarcodes[kmer2].end(); itr_bcSet++){ // for every barcode that supports hap2
    barcodeSupport[*itr_bcSet].first++;
    barcodeSupport[*itr_bcSet].second--;
  }
  return;
}

void extendBarcodeSupport(std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, hashId_t & kmer1, hashId_t & kmer2){
  // insert new haplotype pair [kmer1,kmer2] into the barcode support map
  std::set<std::string>::iterator itr_bcSet;
  for(itr_bcSet=kmerToBarcodes[kmer1].begin(); itr_bcSet!=kmerToBarcodes[kmer1].end(); itr_bcSet++){ // for barcodes that support kmer
    if(barcodeSupport.count(*itr_bcSet)){
      barcodeSupport[*itr_bcSet].first++;
    }else{
      barcodeSupport[*itr_bcSet]=std::make_pair(1,0);
    }
  }
  for(itr_bcSet=kmerToBarcodes[kmer2].begin(); itr_bcSet!=kmerToBarcodes[kmer2].end(); itr_bcSet++){ // for barcodes that support kmer
    if(barcodeSupport.count(*itr_bcSet)){
      barcodeSupport[*itr_bcSet].second++;
    }else{
      barcodeSupport[*itr_bcSet]=std::make_pair(0,1);
    }
  }
  return;
}

void initializeHaplotypes(std::vector<hashId_t> & refKmers, std::vector<hashId_t> & altKmers, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::vector<hashId_t> & hap1, std::vector<hashId_t> & hap2, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport){
  // insert first kmer into haps
  hap1.push_back(refKmers[0]);
  hap2.push_back(altKmers[0]);
  extendBarcodeSupport(barcodeSupport, kmerToBarcodes, refKmers[0], altKmers[0]);
  refKmers[0]=refKmers.back();
  altKmers[0]=altKmers.back();
  refKmers.pop_back();
  altKmers.pop_back();

  // while there are still kmers left
  while(refKmers.size()>0){
    // calculate phasing scores for all remaining kmers
    std::vector<hashId_t>::iterator it_refK = refKmers.begin();
    std::vector<hashId_t>::iterator it_altK = altKmers.begin();
    int32_t bestPos=0, bestScore=0, score;
    for(it_refK; it_refK!=refKmers.end(); it_refK++){
      score = getPhaseScore(*it_refK, *it_altK, barcodeSupport, kmerToBarcodes);
      if(std::abs(score) > std::abs(bestScore)){
        bestScore = score;
        bestPos = std::distance(refKmers.begin(),it_refK);
      }
      it_altK++;
    }
    // insert kmers with best scores into haps
    if(bestScore == 0){ // if kmer not phasable: insert randomly
      score = ((std::rand() % 2) *2) - 1; // returns 1 or -1
    }
    if(bestScore > 0){
      hap1.push_back(refKmers[bestPos]);
      hap2.push_back(altKmers[bestPos]);
      extendBarcodeSupport(barcodeSupport, kmerToBarcodes, refKmers[bestPos], altKmers[bestPos]);
    }else{
      hap1.push_back(altKmers[bestPos]);
      hap2.push_back(refKmers[bestPos]);
      extendBarcodeSupport(barcodeSupport, kmerToBarcodes, altKmers[bestPos], refKmers[bestPos]);
    }
    refKmers[bestPos]=refKmers.back();
    altKmers[bestPos]=altKmers.back();
    refKmers.pop_back();
    altKmers.pop_back();
  }
  return;
}

void phaseKmers(std::vector<hashId_t> & hap1, std::vector<hashId_t> & hap2, std::map<hashId_t, std::set<std::string>> & kmerToBarcodes, std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, phaseOptions & options){
  bool flipped=true;

  while(flipped){
    int32_t bestScore=0, bestPos=0, score;
    flipped = false;

    std::vector<hashId_t>::iterator itr_hap1, itr_hap2=hap2.begin();
    // determine most benefitial haplotype flip
    for(itr_hap1=hap1.begin(); itr_hap1!=hap1.end(); itr_hap1++){ // for every phasing kmer
      score = getPhaseScore(*itr_hap1, *itr_hap2, barcodeSupport, kmerToBarcodes);
      if(score < bestScore){
        bestScore = score;
        bestPos = std::distance(hap1.begin(),itr_hap1);
      }
      itr_hap2++;
    }

    // if bestScore is negative: flip haplotype; else: phasing finished
    if(bestScore < 0){
      flipped = true;
      // update barcodeSupport map
      adjustBarcodeSupport(barcodeSupport, kmerToBarcodes, hap1[bestPos], hap2[bestPos]);
      // flip hap
      std::swap(hap1[bestPos], hap2[bestPos]);
    }
  }
  return;
}

void phaseBarcodes(std::map<std::string, std::pair<uint32_t, uint32_t>> & barcodeSupport, std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2){
  std::map<std::string, std::pair<uint32_t, uint32_t>>::iterator itr_map = barcodeSupport.begin();
  for(itr_map; itr_map!=barcodeSupport.end();itr_map++){
    // std::cerr << itr_map->second.first << "\t" << itr_map->second.second << "\n";
    if(itr_map->second.first > (itr_map->second.second)){
      barcodesHap1.insert(itr_map->first);
    }
    if(itr_map->second.second > (itr_map->second.first)){
      barcodesHap2.insert(itr_map->first);
    }
  }
  return;
}

void touchOutput(phaseOptions & options){
  if(options.paired){
    SeqFileOut out1_hap1(toCString(options.output_prefix+".hap1.1.fq"));
    SeqFileOut out2_hap1(toCString(options.output_prefix+".hap1.2.fq"));
    SeqFileOut out1_hap2(toCString(options.output_prefix+".hap2.1.fq"));
    SeqFileOut out2_hap2(toCString(options.output_prefix+".hap2.2.fq"));
    SeqFileOut out1_unphased(toCString(options.output_prefix+".unphased.1.fq"));
    SeqFileOut out2_unphased(toCString(options.output_prefix+".unphased.2.fq"));
  }else{
    SeqFileOut out_hap1(toCString(options.output_prefix+".hap1.fq"));
    SeqFileOut out_hap2(toCString(options.output_prefix+".hap2.fq"));
    SeqFileOut out_unphased(toCString(options.output_prefix+".unphased.fq"));
  }
  return;
}

void writeToUnphased(ReadData & readData, phaseOptions & options){
  std::set<std::string> barcodesHap1;
  std::set<std::string> barcodesHap2;
  if(options.paired){
    writePhasedOutput_paired(barcodesHap1, barcodesHap2, readData, options);
  }else{
    writePhasedOutput(barcodesHap1, barcodesHap2, readData, options);
  }
  return;
}

void writePhasedOutput(std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2, ReadData & readData, phaseOptions & options){
  // setup iterators
  typedef Iterator<StringSet<Dna5String> >::Type ReadSetIt_t;
  typedef Iterator<StringSet<CharString> >::Type IdSetIt_t;
  ReadSetIt_t itrReads = begin(readData.reads1);
  IdSetIt_t itrIds = begin(readData.ids1);
  IdSetIt_t itrQuals = begin(readData.quals1);

  // open output files
  SeqFileOut out_hap1(toCString(options.output_prefix+".hap1.fq"));
  SeqFileOut out_hap2(toCString(options.output_prefix+".hap2.fq"));
  SeqFileOut out_unphased(toCString(options.output_prefix+".unphased.fq"));

  // iterate over reads
  std::string ReadName;

  for(itrReads; itrReads!=end(readData.reads1); itrReads++){
    // determine Barcode
    ReadName = getReadName(toCString(*itrIds));
    // write to correct files
    if(barcodesHap1.count(ReadName)){
      writeRecord(out_hap1,*itrIds,*itrReads,*itrQuals);
      // if(barcodesHap2.count(ReadName)){
      //   writeRecord(out1_hap2,*itrIds,*itrReads,*itrQuals);
      // }
    }else if(barcodesHap2.count(ReadName)){
      writeRecord(out_hap2,*itrIds,*itrReads,*itrQuals);
    }else{
      writeRecord(out_unphased,*itrIds,*itrReads,*itrQuals);
    }

    ++itrIds;
    ++itrQuals;
  }

  // close files
  close(out_hap1);
  close(out_hap2);
  close(out_unphased);

  // if one hap is empty, write dummy read
  if(barcodesHap1.size()==0){
    std::ifstream  src1(options.output_prefix+".hap2.fq", std::ios::binary);
    std::ofstream  dst1(options.output_prefix+".hap1.fq",   std::ios::binary);
    dst1 << src1.rdbuf();
    src1.close(); dst1.close();
  }
  if(barcodesHap2.size()==0){
    std::ifstream  src1(options.output_prefix+".hap1.fq", std::ios::binary);
    std::ofstream  dst1(options.output_prefix+".hap2.fq",   std::ios::binary);
    dst1 << src1.rdbuf();
    src1. close(); dst1.close();
  }

  return;
}

void writePhasedOutput_paired(std::set<std::string> & barcodesHap1, std::set<std::string> & barcodesHap2, ReadData & readData, phaseOptions & options){
  // setup iterators
  typedef Iterator<StringSet<Dna5String> >::Type ReadSetIt_t;
  typedef Iterator<StringSet<CharString> >::Type IdSetIt_t;
  ReadSetIt_t itrReads1 = begin(readData.reads1);
  ReadSetIt_t itrReads2 = begin(readData.reads2);
  IdSetIt_t itrIds1 = begin(readData.ids1);
  IdSetIt_t itrIds2 = begin(readData.ids2);
  IdSetIt_t itrQuals1 = begin(readData.quals1);
  IdSetIt_t itrQuals2 = begin(readData.quals2);

  // open output files
  SeqFileOut out1_hap1(toCString(options.output_prefix+".hap1.1.fq"));
  SeqFileOut out2_hap1(toCString(options.output_prefix+".hap1.2.fq"));
  SeqFileOut out1_hap2(toCString(options.output_prefix+".hap2.1.fq"));
  SeqFileOut out2_hap2(toCString(options.output_prefix+".hap2.2.fq"));
  SeqFileOut out1_unphased(toCString(options.output_prefix+".unphased.1.fq"));
  SeqFileOut out2_unphased(toCString(options.output_prefix+".unphased.2.fq"));

  // iterate over reads
  std::string Barcode;
  uint_fast8_t barcode_length = getBarcodeLength(toCString(*itrIds1));

  for(itrReads1; itrReads1!=end(readData.reads1); itrReads1++){
    // determine Barcode
    Barcode = getBarcode(toCString(*itrIds1), barcode_length);
    // write to correct files
    if(barcodesHap1.count(Barcode)){
      writeRecord(out1_hap1,*itrIds1,*itrReads1,*itrQuals1);
      writeRecord(out2_hap1,*itrIds2,*itrReads2,*itrQuals2);
      // if(barcodesHap2.count(Barcode)){
      //   writeRecord(out1_hap2,*itrIds1,*itrReads1,*itrQuals1);
      //   writeRecord(out2_hap2,*itrIds2,*itrReads2,*itrQuals2);
      // }
    }else if(barcodesHap2.count(Barcode)){
      writeRecord(out1_hap2,*itrIds1,*itrReads1,*itrQuals1);
      writeRecord(out2_hap2,*itrIds2,*itrReads2,*itrQuals2);
    }else{
      writeRecord(out1_unphased,*itrIds1,*itrReads1,*itrQuals1);
      writeRecord(out2_unphased,*itrIds2,*itrReads2,*itrQuals2);
    }

    ++itrReads2;
    ++itrIds1;
    ++itrIds2;
    ++itrQuals1;
    ++itrQuals2;
  }

  // close files
  close(out1_hap1);
  close(out2_hap1);
  close(out1_hap2);
  close(out2_hap2);
  close(out1_unphased);
  close(out2_unphased);

  // if one hap is empty, write dummy read
  if(barcodesHap1.size()==0){
    std::ifstream  src1(options.output_prefix+".hap2.1.fq", std::ios::binary);
    std::ofstream  dst1(options.output_prefix+".hap1.1.fq",   std::ios::binary);
    dst1 << src1.rdbuf();
    src1.close(); dst1.close();
    std::ifstream  src2(options.output_prefix+".hap2.2.fq", std::ios::binary);
    std::ofstream  dst2(options.output_prefix+".hap1.2.fq",   std::ios::binary);
    dst2 << src2.rdbuf();
    src2.close(); dst2.close();
  }
  if(barcodesHap2.size()==0){
    std::ifstream  src1(options.output_prefix+".hap1.1.fq", std::ios::binary);
    std::ofstream  dst1(options.output_prefix+".hap2.1.fq",   std::ios::binary);
    dst1 << src1.rdbuf();
    src1. close(); dst1.close();
    std::ifstream  src2(options.output_prefix+".hap1.2.fq", std::ios::binary);
    std::ofstream  dst2(options.output_prefix+".hap2.2.fq",   std::ios::binary);
    dst2 << src2.rdbuf();
    src2.close(); dst2.close();
  }

  return;
}
