using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Gadfly, JLD2, Random 

using KidneyAllocation

import KidneyAllocation: reconstruct_recipients, build_recipient_registry, load_recipient, build_donor_registry

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipient = load_recipient(recipient_filepath)

df2 = filter(row -> year(row.CAN_LISTING_DT) < 2013, df)
cand_before_2013 = unique(df2.CAN_ID)

df2 = filter(row -> year(row.CAN_LISTING_DT) < 2020, df)
cand_before_2020 = unique(df2.CAN_ID)

new_recipients = setdiff(cand_before_2020, cand_before_2013)

λᵣ = length(new_recipients)/6












function get_arrival_departure(df::AbstractDataFrame, future_date::Date=Date(2024,1,1))

    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"

    # Sort the dataframe lines so that the most recent is on top
    idx = sortperm(df.UPDATE_TM; rev=true)
    outcomes = uppercase.(String.(df.OUTCOME[idx]))
    updates = df.UPDATE_TM[idx]

    arrival = df.CAN_LISTING_DT[1]

    if outcomes[1] == "1"
        departure = future_date # Arbitrary date after the end of the historic period
    else
        departure = updates[1] # Si transplanté ou retiré
    end

    return arrival, departure
end


G = groupby(df_recipient, :CAN_ID)

arrival, departure = get_arrival_departure(G[1])

arrival = Vector{Date}(undef, length(G))
departure= Vector{Date}(undef, length(G))
can_id = Vector{Int64}(undef, length(G))

for (i, g) in enumerate(G)
    cand_id[i] = g.CAN_ID
    arrival[i] , departure[i] = get_arrival_departure(g)
end

n_actives = Int64[]
for iYear in 2014:2019
    push!(n_actives, count(arrival .< Date(iYear,1,1) .< departure))
end
# Ça ne concorde pas avec le nombre obtenu avec build_recipient_registry dans simulation.jl

n_arrivals = Int64[]
for iYear in 2014:2019
    push!(n_arrivals, count(Date(iYear,1,1) .≤ arrival .≤ Date(iYear,12,31)))
end

n_departures = Int64[]
for iYear in 2014:2019
    push!(n_departures, count(Date(iYear,1,1) .≤ departure .≤ Date(iYear,12,31)))
end

df_donors = KidneyAllocation.load_donor(donor_filepath)
dropmissing!(df_donors, :DON_DEATH_TM)

G = groupby(df_donors, :DON_ID)

date_of_death = Vector{Date}(undef, length(G))
for (i, g) in enumerate(G)
    date_of_death[i] = g.DON_DEATH_TM[1]
end

n_donors = Int64[]
for iYear in 2014:2019
    push!(n_donors, count(Date(iYear,1,1) .≤ date_of_death .≤ Date(iYear,12,31)))
end

df = DataFrame(Year = 2014:2019, Actives = n_actives, Arrivals = n_arrivals, Departures = n_departures, Donors = n_donors)



recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

# expiration = Union{Date, Nothing}[]
expiration = Date[]
for r in recipients
    if !isnothing(r.expiration_date)
    push!(expiration, r.expiration_date)
    end
end

count( year.(expiration) .== 2014)
count( year.(expiration) .== 2015)
count( year.(expiration) .== 2016)
count( year.(expiration) .== 2017)
count( year.(expiration) .== 2018)
count( year.(expiration) .== 2019)