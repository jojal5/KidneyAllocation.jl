using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Gadfly

using KidneyAllocation

## Fonctions utiles pas encore introduites dans la librairie

import KidneyAllocation: check_df_columns, check_df_column_constant, parse_abo

"""
    filter_first_waitlist_episode(df_recipient) -> DataFrame

Return records from the first waiting-list episode of one recipient, ordered by
`:UPDATE_TM`. Rows after the first exit outcome are removed, while the exit
record itself is retained. If no exit outcome is recorded, all rows are kept.
"""
function filter_first_waitlist_episode(df_recipient::AbstractDataFrame)::DataFrame
    check_df_column_constant(df_recipient, :CAN_ID)
    check_df_columns(df_recipient, :OUTCOME, :UPDATE_TM)

    isempty(df_recipient) && return DataFrame(df_recipient)

    exit_outcomes = ("TX", "TX VIVANT", "X", "DCD")
    is_exit(outcome) = uppercase(strip(string(outcome))) ∈ exit_outcomes

    sorted_df = sort(df_recipient, :UPDATE_TM)
    exit_index = findfirst(is_exit, sorted_df.OUTCOME)
    last_index = something(exit_index, nrow(sorted_df))

    return reverse(sorted_df[1:last_index, :])
end


## Chargement des données brutes

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients = KidneyAllocation.load_recipient(recipient_filepath)
df_donors = KidneyAllocation.load_donor(donor_filepath)
df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])


## Récupération des candidats dont au moins un score a été calculé pour la liste 5 entre 2016 et 2020
# Reconstruction des recveurs par le fichier des donneurs

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2024,12,31), filtered_donors)
G = groupby(filtered_donors, :DON_ID)



G[1]

df = DataFrame(DON_AGE = Int[], DON_BLOOD = String[], ETHNICITY = Int[], HYPERTENSION = Bool[], DIABETES = Bool[],
    CVA = Bool[], CREATININE = Int[], HCV = Bool[], DCD = Bool[], KDRI = Float64[])
# allowmissing!(df)

for g in G
    r = first(g)

    hypertension = r.HYPERTENSION == 1
    diabetes = r.DIABETES == 1
    cva = r.DEATH ∈ (4, 16)
    creatinine = ismissing(r.CREATININE) ? missing : KidneyAllocation.creatinine_mgdl(r.CREATININE)
    dcd = r.DCD == 1

    try
        # If non-missing values for kdri computation
        kdri = KidneyAllocation.evaluate_kdri(r.DON_AGE, r.HEIGHT, r.WEIGHT, hypertension, diabetes, cva, creatinine, dcd)
        push!(df, [r.DON_AGE, r.DON_BLOOD, r.ETHNICITY, r.HYPERTENSION, r.DIABETES, cva, r.CREATININE, r.HCV, dcd, kdri])
    catch
        continue
    end

end


X = KidneyAllocation.construct_feature_matrix_from_df(df, features)
y = df.KDRI


model = DecisionTree.build_tree(y, X)
# apply learned model
# apply_tree(model, [-0.9,3.0,5.1,1.9,0.0])
# run 3-fold cross validation, returns array of coefficients of determination (R^2)
n_folds = 3
r2 = nfoldCV_tree(y, X, n_folds)

