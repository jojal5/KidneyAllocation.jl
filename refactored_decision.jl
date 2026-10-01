using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, build_donor_registry

recipients_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
donors_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipient_dict = build_recipient_registry(recipients_filepath, cpra_filepath)

df = KidneyAllocation.load_recipient(recipients_filepath)

recipient_dict[3]

import KidneyAllocation: check_df_column_constant, check_df_columns
