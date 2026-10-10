using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Gadfly, Statistics

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

import KidneyAllocation: retrieve_observed_waiting_list, load_recipient, load_donor

## Candidates on the waiting list

t = Date(2012,1,1):Month(1):Date(2024,1,1)
n = Vector{Int64}(undef, length(t))

for (i,tᵢ) in enumerate(t)
    ids = retrieve_observed_waiting_list(recipient_filepath, tᵢ)
    n[i] = length(ids)
end

df = DataFrame(Date = t, Candidates = n)

plot(df, x=:Date, y=:Candidates)


## Monthly donor arrivals

df_donors = load_donor(donor_filepath)

df_donor_arrival = DataFrame(DON_ID = Int64[], DON_DEATH_TM = DateTime[])

for g in groupby(df_donors, :DON_ID)
    push!(df_donor_arrival, [first(g.DON_ID), first(g.DON_DEATH_TM)])
end
# groupby(df_donors, :DON_DEATH_TM => yearmonth)

df_donor_arrival.year = year.(df_donor_arrival.DON_DEATH_TM)
df_donor_arrival.month = month.(df_donor_arrival.DON_DEATH_TM)

df = combine(groupby(df_donor_arrival, [:year, :month]), :DON_DEATH_TM => length => :Donors)

d = [Date(r.year, r.month, 1) for r in eachrow(df) ]
df.Date = d

plot(df, x=:Date, y=:Donors)


## Monthly candidate arrivals

df_recipients = load_recipient(recipient_filepath)

df_recipient_arrival = DataFrame(CAN_ID = Int64[], CAN_LISTING_DT = DateTime[])

for g in groupby(df_recipients, :CAN_ID)
    push!(df_recipient_arrival, [first(g.CAN_ID), first(g.CAN_LISTING_DT)])
end

df_recipient_arrival.year = year.(df_recipient_arrival.CAN_LISTING_DT)
df_recipient_arrival.month = month.(df_recipient_arrival.CAN_LISTING_DT)

df = combine(groupby(df_recipient_arrival, [:year, :month]), :CAN_LISTING_DT => length => :Candidates)

d = [Date(r.year, r.month, 1) for r in eachrow(df) ]
df.Date = d

plot(df, x=:Date, y=:Candidates)



## Retrieve transplanted recipients for prediction

df = filter(row -> row.OUTCOME == "TX" &&
    row.CAN_LISTING_DT ≥ Date(2016,1,1) &&
    row.UPDATE_TM < Date(2020,1,1),
    df_recipients)


tx_ids = df.CAN_ID

df.time_to_transplant = Dates.days.(df.UPDATE_TM - df.CAN_LISTING_DT)

time_to_first_offer = Vector{Union{Int64, Missing}}(undef, nrow(df))

for (i,r) in enumerate(eachrow(df))
    
    println(i)
    
    # extract the first offer after listing_date

    df_offers = filter(row-> row.CAN_ID == r.CAN_ID &&
    row.DON_DEATH_TM ≥ r.CAN_LISTING_DT &&
    row.DON_DEATH_TM ≤ r.UPDATE_TM, df_donors)

    if !isempty(df_offers)
        time_to_first_offer[i] = Dates.days(minimum(df_offers.DON_DEATH_TM) - r.CAN_LISTING_DT)
    end

end

count(ismissing.(time_to_first_offer))
# Pour 40 patients transplantés, on ne retrouve pas les offres dans le fichier donors.

df.time_to_first_offer = time_to_first_offer

# On remplace le temps avant la première offre avec le UPDATE_TM du fichier des receveurs.
# for (i,r) in eachrow(df)
#     if ismissing(r.time_to_first_offer)
#         df
# end

df[!, :time_to_first_offer] = coalesce.(df[!, :time_to_first_offer], df[!, :time_to_transplant])

df

combine(groupby(df, :CAN_BLOOD), :time_to_first_offer => mean => :mean_waiting_time)





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