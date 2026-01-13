using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames


import KidneyAllocation: build_recipient_registry, load_recipient

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

build_recipient_registry(recipient_filepath, cpra_filepath)


# Estimate the recipient arrival rate
df = load_recipient(recipient_filepath)
df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)
λᵣ = length(unique(df2.CAN_ID)) / 6






import KidneyAllocation: load_donor, build_donor_registry

filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df = load_donor(filepath)

dropmissing!(df)

donors = build_donor_registry(filepath)

# Estimate the donor arrival rate
df = load_donor(filepath)
df2 = filter(row -> 2014 ≤ year(row.DON_DEATH_TM) < 2020, df)
filter!(row -> row.DECISION =="Acceptation", df2)
λₒ = length(unique(df2.CAN_ID)) / 6