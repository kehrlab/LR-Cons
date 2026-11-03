CXX=g++
CC=g++

INCLUDE_PATH=./include
BINARY_PATH=bin

ifdef INCLUDE_PATH
	override CXXFLAGS += -I$(INCLUDE_PATH)
endif

override CXXFLAGS+=-std=c++23 -lhts -lz -lspoa -labpoa -fopenmp -lpthread -g

all: create_consensus extract_region_sequences evaluate_consensus_seqs extract_region_sv_stats

create_consensus: src/create_consensus.cpp
	@$(CXX) src/create_consensus.cpp $(CXXFLAGS) -o $(BINARY_PATH)/$@

extract_region_sequences: src/extract_region_sequences.cpp src/record.cpp
	@$(CXX) src/extract_region_sequences.cpp src/record.cpp $(CXXFLAGS) -o $(BINARY_PATH)/$@

evaluate_consensus_seqs: src/evaluate_consensus_seqs.cpp
	@$(CXX) src/evaluate_consensus_seqs.cpp $(CXXFLAGS) -o $(BINARY_PATH)/$@

extract_region_sv_stats: src/extract_region_sv_stats.cpp
	@$(CXX) src/extract_region_sv_stats.cpp $(CXXFLAGS) -o $(BINARY_PATH)/$@

clean:
	if [ -f $(BINARY_PATH)/create_consensus ]; then rm $(BINARY_PATH)/create_consensus; fi
	if [ -f $(BINARY_PATH)/determine_consensus_times ]; then rm $(BINARY_PATH)/determine_consensus_times; fi
	if [ -f $(BINARY_PATH)/extract_region_sequences ]; then rm $(BINARY_PATH)/extract_region_sequences; fi
	if [ -f $(BINARY_PATH)/evaluate_consensus_seqs ]; then rm $(BINARY_PATH)/evaluate_consensus_seqs; fi
