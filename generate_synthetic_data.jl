# In terminal, run the following command
# julia --project=. --threads=8 generate_synthetic_data.jl

# using Pkg
# Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2

using KidneyAllocation

import KidneyAllocation: build_donor_registry, build_recipient_registry
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold



recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
donors = build_donor_registry(donor_filepath)

# Fit decision model
data = retrieve_decision_data(donor_filepath, recipient_filepath)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

fm = glm(model, data, Bernoulli(), LogitLink())

u = fit_decision_threshold(fm)


using Base.Threads

outdir = joinpath(@__DIR__, "src", "SyntheticData")
mkpath(outdir)

Threads.@threads for i in 1:1000

    recipient_list, donor_list =
        KidneyAllocation.generate_pseudo_history(recipients, donors, fm, u)

    filename = string("sim_", lpad(i, 4, '0') ,".jld2")
    filepath = joinpath(outdir, filename)

    jldsave(filepath; recipient_list, donor_list) 
end