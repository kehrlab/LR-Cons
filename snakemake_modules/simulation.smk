include: "simulation/simulate_data.smk"
include: "simulation/preprocess_data.smk"

include: "simulation/main/preprocess.smk"
include: "simulation/main/consensus.smk"

include: "simulation/kmer/consensus_k.smk"

include: "simulation/uncertainty/preprocess_uncertainty.smk"
include: "simulation/uncertainty/consensus_uncertainty.smk"

include: "simulation/giab/giab_preprocess.smk"
include: "simulation/giab/giab_consensus.smk"
include: "simulation/giab/giab_simulate.smk"
