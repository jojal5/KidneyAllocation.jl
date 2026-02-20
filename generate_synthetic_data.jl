# In terminal, run the following command
# julia --project=. --threads=8 generate_synthetic_data.jl

using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, build_donor_registry, simulate_initial_state_indexed

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
donors = build_donor_registry(donor_filepath)

# Load fitted decision model
@load "src/SyntheticData/TreeDecisionModel.jld2"


using Base.Threads

output_filename = "initial_waiting_lists_indexed.jld2"
output_dir = joinpath(@__DIR__, "src", "SyntheticData")
mkpath(output_dir)
output_path = joinpath(output_dir, output_filename)

num_simulations = 2
base_seed = 1234 

origin_date = Date(2000,1,1)

# Each entry i holds the indexed waiting list produced by simulation i
waiting_indices = Vector{Vector{UInt16}}(undef, num_simulations)
waiting_day_offsets = Vector{Vector{Int16}}(undef, num_simulations)

Threads.@threads for sim_id in 1:num_simulations
    rng = MersenneTwister(base_seed + sim_id) # Ensure that all simulation are independent

    idx, dates = simulate_initial_state_indexed(recipients, donors, fm, u; origin_date=origin_date, rng=rng)

    waiting_indices[sim_id] = UInt16.(idx)
    waiting_day_offsets[sim_id] = Int16.(Dates.value.(dates .- origin_date))
end

jldsave(output_path; waiting_indices, waiting_day_offsets, origin_date)

