# In terminal, run the following command
# julia --project=. --threads=12 waiting_time_estimation.jl

using Base.Threads

using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

import KidneyAllocation: count_donor_arrivals, count_recipient_arrivals, fractionalyears_between

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

## Load decision model

@load "src/SyntheticData/TreeDecisionModel.jld2"


## Estimation of the recipient arrival rate

start_time = DateTime(2012,1,1,0,0,0)
end_time = DateTime(2020,1,1,0,0,0)

df_recipient = KidneyAllocation.load_recipient(recipient_filepath)
recipient_arrival_rate = count_recipient_arrivals(df_recipient, start_time, end_time) / fractionalyears_between(start_time, end_time)


## Estimation of the donor arrival rate

start_time = DateTime(2012,1,1,0,0,0)
end_time = DateTime(2020,1,1,0,0,0)

df_donor = KidneyAllocation.load_donor(donor_filepath)
donor_arrival_rate = count_donor_arrivals(df_donor, start_time, end_time) / fractionalyears_between(start_time, end_time)


## Build recipient registry by CAN_ID

recipient_registry = KidneyAllocation.build_recipient_registry(recipient_filepath, cpra_filepath)


## Build donor registry by DON_ID

donor_registry = KidneyAllocation.build_donor_registry(donor_filepath)

## Dictionary of recipient last status

last_status = Dict{Int64, Tuple{DateTime, String}}()

for g in groupby(df_recipient, :CAN_ID)
    last_status[g.CAN_ID[1]] = KidneyAllocation.get_last_update(g)
end

last_status[3]

## Retrieve transplanted recipients for which waiting time has to be estimated

# CPRA < 80
# Registered between 2015 and 2020
# Transplanted before 2020
# Only for their first registration if past transplant

recipient_ids = Int64[]

for id in keys(recipient_registry)
    r = recipient_registry[id]
    date = first(last_status[id])
    status = uppercase(strip(last(last_status[id])))

    if r.arrival ≥ Date(2015,1,1) && r.arrival < Date(2020,1,1)
        if r.cpra ≤ 80
            if status == "TX" && date < Date(2020,1,1)
                push!(recipient_ids, id)
            end
        end
    end  
end

## Selection of a recipient and simulation

    nyears = 10
    nsim = 1000

Threads.@threads for id in recipient_ids

    recipient = recipient_registry[id]

    offers = KidneyAllocation.simulate_recipient_offers(
        recipient,
        recipient_filepath,
        recipient_registry,
        donor_filepath,
        donor_registry,
        dm,
        nyears,
        nsim;
        recipient_arrival_rate=recipient_arrival_rate,
        donor_arrival_rate=donor_arrival_rate
        )

    df = KidneyAllocation.offers_to_dataframes(recipient, offers)

    filename = string("/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/",id,".csv")

    CSV.write(filename, df)

end




