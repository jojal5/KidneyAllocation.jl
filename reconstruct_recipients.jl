using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using KidneyAllocation


## Chargement des données 

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

# On ne retient que la première inscription des candidats
df_donors = unique( sort(df_donors, :DON_DEATH_TM), :CAN_ID)
# Ça fait qu'on ignore les patients qui reviennent et les arrivées des reins qui leurs sont attribués

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







for row in eachrow(df)
    CAN_ID = row.CAN_ID
    DON_ID = row.DON_ID

    if haskey(recipient_registry, CAN_ID) && haskey(donor_registry, DON_ID)

        recipient = recipient_registry[row.CAN_ID]
        donor = donor_registry[row.DON_ID]


        donor_age = donor.age
        kdri = donor.kdri
        recipient_age = KidneyAllocation.years_between(recipient.birth, donor.arrival)
        recipient_blood = recipient.blood
        waittime = KidneyAllocation.fractionalyears_between(recipient.dialysis, donor.arrival)
        s = row.DON_CAN_SCORE
        mismatch = KidneyAllocation.mismatch_count(recipient, donor)
        cpra = recipient.cpra
        decision = row.DECISION == "Acceptation"

        push!(df_decision, [DON_ID, donor_age, kdri, recipient_age, recipient_blood, waittime, s, cpra, mismatch, decision])

    end
end