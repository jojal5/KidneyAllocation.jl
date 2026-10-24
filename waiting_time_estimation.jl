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

donor_registry, given_kidneys_registry = KidneyAllocation.build_donor_registry(df_donors) 
recipient_registry, time_to_first_offer_registry, time_to_transplant_registry = KidneyAllocation.build_recipient_registry(df_recipients)


## Load decision model

@load "src/SyntheticData/TreeDecisionModel.jld2"


## Estimation of the recipient arrival rate

start_time = Date(2016,1,1)
end_time = Date(2020,1,1)



import KidneyAllocation: check_df_columns, fractionalyears_between
"""
    count_recipient_arrivals(df, start_time, end_time) -> Int

Count recipients whose listing time lies in `[start_time, end_time)`.
"""
function count_recipient_arrivals(df::AbstractDataFrame, start_time::Date, end_time::Date)::Int
    check_df_columns(df, :CAN_ID, :CAN_LISTING_DT)

    start_time ≤ end_time ||
        throw(ArgumentError("`start_time` must not be after `end_time`"))

    n_recipients = 0

    for recipient_df in groupby(df, :CAN_ID)
        listing_time = first(recipient_df.CAN_LISTING_DT)

        if start_time ≤ listing_time < end_time
            n_recipients += 1
        end
    end

    return n_recipients
end

"""
    count_donor_arrivals(df, start_time, end_time) -> Int

Count donors whose death time lies in `[start_time, end_time)`.
"""
function count_donor_arrivals(df::AbstractDataFrame, start_time::Date, end_time::Date)::Int
    check_df_columns(df, :DON_ID, :DON_DEATH_TM)

    start_time ≤ end_time ||
        throw(ArgumentError("`start_time` must not be after `end_time`"))

    n_donors = 0

    for donor_df in groupby(df, :DON_ID)
        donor_arrival_time = first(donor_df.DON_DEATH_TM)

        if start_time ≤ donor_arrival_time < end_time
            n_donors += 1
        end
    end

    return n_donors
end


recipient_arrival_rate = count_recipient_arrivals(df_recipients, start_time, end_time) / fractionalyears_between(start_time, end_time)


## Estimation of the donor arrival rate


donor_arrival_rate = count_donor_arrivals(df_donors, start_time, end_time) / fractionalyears_between(start_time, end_time)

## Retrieve transplanted recipients for which waiting time has to be estimated

# CPRA < 80
# Registered between 2016 and 2020
# Transplanted before 2020

TX_ids = Int64[]

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

TX_ids = df_TX.CAN_ID


## Selection of a recipient and simulation

    nyears = 5
    nsim = 5

Threads.@threads for id in TX_ids

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

    # filename = string("/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/",id,".jld2")
    # jldsave(filename; offers)

    df = KidneyAllocation.offers_to_dataframes(recipient, offers, dm)
    filename = string("/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/",id,".csv")
    CSV.write(filename, df)

end

