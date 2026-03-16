CXX=g++
CC=g++

INCLUDE_PATH=./include

ifdef INCLUDE_PATH
	override CXXFLAGS += -I$(INCLUDE_PATH)
endif

override CXXFLAGS+=-std=c++23 -lhts -lz -lspoa -labpoa -fopenmp -lpthread -g

all: create_consensus extract_region_sequences evaluate_consensus_seqs

create_consensus: src/create_consensus.cpp
	@$(CXX) src/create_consensus.cpp $(CXXFLAGS) -o $@

determine_consensus_times: src/create_consensus_time_benchmark.cpp
	@$(CXX) src/create_consensus_time_benchmark.cpp $(CXXFLAGS) -o $@

extract_region_sequences: src/extract_region_sequences.cpp src/record.cpp
	@$(CXX) src/extract_region_sequences.cpp src/record.cpp $(CXXFLAGS) -o $@

evaluate_consensus_seqs: src/evaluate_consensus_seqs.cpp
	@$(CXX) src/evaluate_consensus_seqs.cpp $(CXXFLAGS) -o $@

clean:
	if [ -f create_consensus ]; then rm create_consensus; fi
	if [ -f determine_consensus_times ]; then rm determine_consensus_times; fi
	if [ -f extract_region_sequences ]; then rm extract_region_sequences; fi
	if [ -f evaluate_consensus_seqs ]; then rm evaluate_consensus_seqs; fi
