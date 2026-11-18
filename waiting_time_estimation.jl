# In terminal, run the following command
# julia --project=. --threads=16 waiting_time_estimation.jl

using Base.Threads

using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

## Load data

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

# On ne retient que les patients majeurs
filter!(row -> KidneyAllocation.years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df_recipients)

# On ne retient que la première inscription des candidats
# df_donors = unique( sort(df_donors, :DON_DEATH_TM), :CAN_ID)
# Ça fait qu'on ignore les patients qui reviennent et les arrivées des reins qui leurs sont attribués

donor_registry, given_kidneys_registry = KidneyAllocation.build_donor_registry(df_donors)
recipient_registry = KidneyAllocation.build_recipient_registry(df_recipients)


## Load decision model

# @load "src/SyntheticData/TreeDecisionModel.jld2"
@load "src/SyntheticData/GLMDecisionModel.jld2"


## Estimation of the recipient arrival rate

start_time = Date(2016,1,1)
end_time = Date(2020,1,1)



import KidneyAllocation: fractionalyears_between, count_donor_arrivals, count_recipient_arrivals


recipient_arrival_rate = count_recipient_arrivals(df_recipients, start_time, end_time) / fractionalyears_between(start_time, end_time)


## Estimation of the donor arrival rate


donor_arrival_rate = count_donor_arrivals(df_donors, start_time, end_time) / fractionalyears_between(start_time, end_time)

## Retrieve transplanted recipients for which waiting time has to be estimated

# CPRA < 80
# Registered between 2016 and 2020
# Transplanted before 2020

# TX_ids = Int64[]

# df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
# dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
# filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
# filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
# filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
# filter!(row->row.CPRA < 80, df_TX)

# TX_ids = df_TX.CAN_ID

## Retrieve untransplanted recipients for which waiting time has to be estimated

# CPRA < 80
# Registered between 2016 and 2020
# Removed before 2020

TX_ids = Int64[]

removed_status = ["X", "TX VIVANT", "DCD"]

df_TX = filter(row-> uppercase(strip(row.OUTCOME)) ∈ removed_status, df_recipients)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

filter!(row -> row.CAN_ID ∈ keys(recipient_registry), df_TX)

TX_ids = df_TX.CAN_ID


## Selection of a recipient and simulation

nyears = 5
nsim = 1000

Threads.@threads for id in TX_ids
# for id in TX_ids
# id = 13707

    recipient = recipient_registry[id]

    # offers = KidneyAllocation.simulate_recipient_offers(
    #     recipient,
    #     recipient_filepath,
    #     recipient_registry,
    #     donor_filepath,
    #     donor_registry,
    #     dm,
    #     nyears,
    #     nsim;
    #     recipient_arrival_rate=recipient_arrival_rate,
    #     donor_arrival_rate=donor_arrival_rate
    #     )

    offers = KidneyAllocation.simulate_recipient_offers(recipient::Recipient,
        df_recipients,
        recipient_registry,
        donor_registry,
        given_kidneys_registry,
        dm::AbstractDecisionModel;
        nyears = nyears,
        nsim = nsim,
        recipient_arrival_rate=recipient_arrival_rate,
        donor_arrival_rate=donor_arrival_rate
        )

    df = KidneyAllocation.offers_to_dataframes(recipient, offers, dm)
    filename = string("/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/",id,".csv")
    CSV.write(filename, df)

end


