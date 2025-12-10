using KidneyAllocation
using Test

using Dates, CSV, DataFrames

df = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv", DataFrame)

unique(df.OUTCOME)

df_outcome = unique(select(df, :OUTCOME))

CSV.write("outcome.csv", df_outcome)