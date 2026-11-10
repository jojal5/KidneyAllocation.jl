using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using KidneyAllocation


## Chargement des données 

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

## build_donor_registry

donor_registry, given_kidneys = KidneyAllocation.build_donor_registry(df_donors)
recipient_registry = KidneyAllocation.build_recipient_registry(df_recipients)


count(skipmissing(df_recipients.SCORE_MATCHING_FLAG))
count(skipmissing(.!(df_recipients.SCORE_MATCHING_FLAG)))
count(ismissing.(df_recipients.SCORE_MATCHING_FLAG))



## Récupération des décisions

df_offers = KidneyAllocation.load_donor(donor_filepath)
filter!(row -> row.ATT_TYPE == 5, df_offers)

df_offers_filtered = similar(df_offers, 0)

for df_offer in groupby(df_offers, :DON_ID)

    filtered_df = KidneyAllocation.filter_offers_until_last_acceptance(df_offer)
    append!(df_offers_filtered, filtered_df)

end


df_offers_filtered

## Time to tx
# Computed in elapsed days since dialysis

time_to_transplant = Dict{Int, Int}()

for row in eachrow(df_donors)
    if haskey(recipient_registry, row.CAN_ID)
        recipient = recipient_registry[row.CAN_ID]
        elapsed_days = Dates.days(Dates.Date(row.DON_DEATH_TM) - recipient.dialysis)
        time_to_transplant[row.CAN_ID] = elapsed_days
    end
end


## Time to first offer
# Computed in elapsed days since dialysis

time_to_first_offer = Dict{Int, Int}()

df_offers_filtered

for g in groupby(df_offers_filtered, :CAN_ID)

    CAN_ID = first(g.CAN_ID)

    if haskey(recipient_registry, CAN_ID)
        recipient = recipient_registry[CAN_ID]
        elapsed_days = Dates.days(Dates.Date(minimum(g.DON_DEATH_TM)) - recipient.dialysis)
        time_to_first_offer[CAN_ID] = elapsed_days
    end

end
