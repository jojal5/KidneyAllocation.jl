using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Gadfly

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

donor_registry, given_kidneys_registry = KidneyAllocation.build_donor_registry(df_donors) 
recipient_registry, time_to_first_offer_registry, time_to_transplant_registry = KidneyAllocation.build_recipient_registry(df_recipients)

## TX recipients 

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

TX_ids = df_TX.CAN_ID


## Verify if the computed score matches the recorded one. If not, it indicates wrong recipients entries in the original file data.


score_matches = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # id = TX_ids[1]

    println(id)

    df_CAN_DON = filter(row -> row.CAN_ID == id, df_donors)
    filter!(row -> row.DECISION == "Acceptation", df_CAN_DON)

    if isempty(df_CAN_DON)
        continue
    end
    
    don_id = first(df_CAN_DON.DON_ID)
    recorded_score = first(df_CAN_DON.DON_CAN_SCORE)

    # Vérification du score
    computed_score = KidneyAllocation.score(donor_registry[don_id], recipient_registry[id])

    if isapprox(recorded_score, computed_score, atol = .1)
        score_matches[i] = true
    end

end

count(score_matches)


id = 3655
df_CAN_DON = filter(row -> row.CAN_ID == id, df_donors)
filter!(row -> row.DECISION == "Acceptation", df_CAN_DON)


# Le candidat 3655 a été transplanté d'un rein de la liste 4.

## Check recipient by recipient

# id = 3748
# id = 10615
id = TX_ids[3]
id = TX_ids[1]

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/"
filepath = joinpath(sim_folder, "$id.csv")

res = CSV.read(filepath, DataFrame)

df = combine(groupby(res, :simulation_id), :elapsed_days_since_dialysis => minimum => :elapsed_days)

plot(df, y=:elapsed_days, Geom.boxplot)

time_to_first_offer_registry[id]
time_to_transplant_registry[id]


filter(row -> row.CAN_ID == id, df_recipients)
df2 = filter(row -> row.CAN_ID == id, df_donors)
filter!(row -> row.DECISION == "Acceptation", df2)

don_id = first(df2.DON_ID)
id

# Vérification du score
s = KidneyAllocation.score(donor_registry[don_id], recipient_registry[id])
isapprox(first(df2.DON_CAN_SCORE), s; atol = .1) 


# Le score enregistré est de 10.23, on dirait qu'il y a erreur dans les entrées du fichier Candidate
KidneyAllocation.score(donor_registry[5382], recipient_registry[1678])







## Check calibration for all recipients - First offer

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    println(i)

    # df_recipient = filter(row -> row.CAN_ID == id, df_recipients)
    # df_donor = filter(row -> row.CAN_ID == id, df_donors)

    # if isempty(df_recipient) || isempty(df_donor)
    #     n_empty +=1
    #     continue
    # end

    obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    # obs_waiting_time_before_first_offer = time_to_transplant_registry[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)


## Check calibration for all recipients - Time to transplant

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    println(i)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = time_to_transplant_registry[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    ind_accepted = res.acceptance_probability .≥ rand(nrow(res))
    accepted_offers = res[ind_accepted, :]

    t1 = first_qualifying_offer_quantile_given_offer(accepted_offers, .25)
    t2 = first_qualifying_offer_quantile_given_offer(accepted_offers, .5)
    t3 = first_qualifying_offer_quantile_given_offer(accepted_offers, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)


## Check calibration for all recipients - Time to first offer compared with observed time to transplant

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    println(i)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = time_to_transplant_registry[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    # ind_accepted = res.acceptance_probability .≥ rand(nrow(res))
    # accepted_offers = res[ind_accepted, :]

    t1 = first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)


σ = sqrt(p₂*(1. .-  p₂)/length(q2))
p₂ + 1.96*σ

σ = sqrt(p₃*(1. .-  p₃)/length(q3))
p₃ + 1.96*σ


## Check calibration for permanently removed patients

TX_ids = Int64[]

removed_status = ["X", "TX VIVANT", "DCD"]

df_TX = filter(row-> uppercase(strip(row.OUTCOME)) ∈ removed_status, df_recipients)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

filter!(row -> row.CAN_ID ∈ keys(recipient_registry), df_TX)

TX_ids = df_TX.CAN_ID

id = TX_ids[3]

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/"
filepath = joinpath(sim_folder, "$id.csv")

res = CSV.read(filepath, DataFrame)

df = combine(groupby(res, :simulation_id), :elapsed_days => minimum => :elapsed_days)

plot(df, y=:elapsed_days, Geom.boxplot)

time_to_first_offer_registry[id]
time_to_transplant_registry[id]


filter(row -> row.CAN_ID == id, df_recipients)
filter(row -> row.CAN_ID == id, df_donors)














## Time to transplant

filepath = joinpath(sim_folder, "$id.csv")
res = CSV.read(filepath, DataFrame)

ind_accepted = res.acceptance_probability .≥ rand(nrow(res))
accepted_offers = res[ind_accepted, :]

df = combine(groupby(accepted_offers, :simulation_id), :elapsed_days => minimum => :elapsed_days)

plot(df, y=:elapsed_days, Geom.boxplot)

time_to_first_offer_registry[id]
time_to_transplant_registry[id]



## First transplant only

df = CSV.read(recipient_filepath, DataFrame, missingstring=["", "-", "NULL"])
exit_outcomes = ["TX", "TX VIVANT", "X", "DCD"]

# df_first_transplant = similar(df, 0)

first_transplant_ids = Int64[]

for g in groupby(df, :CAN_ID)
    n_outcome = 0
    for r in eachrow(g)
        if uppercase(strip(r.OUTCOME)) ∈ exit_outcomes
            n_outcome += 1
        end
    end
    if n_outcome ≤ 1
        push!(first_transplant_ids, first(g.CAN_ID))
    end
end

first_transplant_ids

TX_ids_filtered = Int64[]

for id in TX_ids
    if id ∈ first_transplant_ids
        push!(TX_ids_filtered, id)
    end
end

TX_ids_filtered

g = filter(row -> row.CAN_ID == 284, df)




println("fin")



# """
#     count_recipient_arrivals(df, start_time, end_time) -> Int

# Return the number of distinct recipients whose listing date lies in
# `[start_time, end_time)`.
# """
# function count_recipient_arrivals(df::AbstractDataFrame, start_time::DateTime, end_time::DateTime)::Int

#     check_df_columns(df, :CAN_ID, :CAN_LISTING_DT)

#     start_time ≤ end_time || throw(ArgumentError("`start_time` must not be after `end_time`"))

#     n_recipients = 0

#     for recipient_df in groupby(df, :CAN_ID)
#         listing_time = DateTime(first(recipient_df.CAN_LISTING_DT))

#         if start_time ≤ listing_time < end_time
#             n_recipients += 1
#         end
#     end

#     return n_recipients
# end


# """
#     count_donor_arrivals(df, start_time, end_time) -> Int

# Return the number of distinct donors whose death time lies in
# `[start_time, end_time)`.
# """
# function count_donor_arrivals(
#     df::AbstractDataFrame,
#     start_time::DateTime,
#     end_time::DateTime,
# )::Int

#     check_df_columns(df, :DON_ID, :DON_DEATH_TM)

#     start_time ≤ end_time ||
#         throw(ArgumentError("`start_time` must not be after `end_time`"))

#     n_donors = 0

#     for donor_df in groupby(df, :DON_ID)
#         donor_arrival_time = DateTime(first(donor_df.DON_DEATH_TM))

#         if start_time ≤ donor_arrival_time < end_time
#             n_donors += 1
#         end
#     end

#     return n_donors
# end





















donor_registry, given_kidneys = KidneyAllocation.build_donor_registry(df_donors) 
recipient_registry, time_to_first_offer_registry, time_to_transplant_registry = KidneyAllocation.build_recipient_registry(df_recipients)



res = CSV.read("/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/11634.csv", DataFrame)

nsim = 1000
nyears = 10

"""
    km_quantile(times, events, p) -> Union{Float64,Nothing}

Estimate the `p`-th quantile of an event time from right-censored observations
using the Kaplan–Meier estimator.

Return `nothing` if the estimated event probability does not reach `p`.
"""
function km_quantile(
    times::AbstractVector{<:Real},
    events::AbstractVector{Bool},
    p::Real,
)::Union{Float64,Nothing}

    length(times) == length(events) ||
        throw(ArgumentError("`times` and `events` must have the same length"))
    0 < p < 1 ||
        throw(ArgumentError("`p` must lie strictly between 0 and 1"))
    all(t -> t ≥ 0, times) ||
        throw(ArgumentError("`times` must be non-negative"))

    event_times = sort(unique([
        times[i] for i in eachindex(times) if events[i]
    ]))

    survival = 1.0
    target_survival = 1 - p

    for time in event_times
        n_at_risk = count(t -> t ≥ time, times)
        n_events = count(
            i -> events[i] && times[i] == time,
            eachindex(times),
        )

        survival *= 1 - n_events / n_at_risk

        survival ≤ target_survival && return Float64(time)
    end

    return nothing
end

import KidneyAllocation.check_df_columns

"""
    first_qualifying_offer_quantile(df, p; max_kdri, n_simulations,
                                    censoring_time) -> Union{Float64,Nothing}

Estimate the `p`th quantile of time to the first offer with KDRI at most
`max_kdri`. Simulations without such an offer are right-censored at
`censoring_time`.
"""
function first_qualifying_offer_quantile(
    df::AbstractDataFrame,
    p::Real;
    max_kdri::Real=Inf,
    n_simulations::Integer=1000,
    censoring_time::Real=Inf,
)::Union{Float64,Nothing}

    check_df_columns(df, :simulation_id, :elapsed_days, :kdri)

    n_simulations > 0 ||
        throw(ArgumentError("`n_simulations` must be positive"))
    censoring_time ≥ 0 ||
        throw(ArgumentError("`censoring_time` must be non-negative"))
    0 < p < 1 ||
        throw(ArgumentError("`p` must lie strictly between 0 and 1"))
    all(id -> id isa Integer && 1 ≤ id ≤ n_simulations, df.simulation_id) ||
        throw(ArgumentError("`:simulation_id` must lie in 1:$n_simulations"))
    all(t -> 0 ≤ t ≤ censoring_time, df.elapsed_days) ||
        throw(ArgumentError("`:elapsed_days` must lie in [0, censoring_time]"))

    times = fill(Float64(censoring_time), n_simulations)
    events = falses(n_simulations)

    qualifying_df = filter(row -> row.kdri ≤ max_kdri, df)

    for simulation_df in groupby(qualifying_df, :simulation_id)
        simulation_id = first(simulation_df.simulation_id)

        times[simulation_id] = minimum(simulation_df.elapsed_days)
        events[simulation_id] = true
    end

    return km_quantile(times, events, p)
end




function time_to_first_offer(df::AbstractDataFrame, p::Real; max_kdri::Real=Inf, n_simulations::Integer=1000)
    check_df_columns(df, :simulation_id, :elapsed_days, :kdri)

    all(id -> id isa Integer && 1 ≤ id ≤ n_simulations, df.simulation_id) ||
        throw(ArgumentError("`:simulation_id` must lie in 1:$n_simulations"))

    0 < p < 1 ||
        throw(ArgumentError("`p` must lie strictly between 0 and 1"))

    filtered_df = filter(row->row.kdri < max_kdri, df)

    if isempty(filtered_df)
        return Inf
    else

        df_min_time = combine(groupby(filtered_df, :simulation_id), :elapsed_days => minimum => :elapsed_days)

        t = fill(Inf, n_simulations)

        for g in groupby(df_min_time, :simulation_id)
            t[first(g.simulation_id)] = first(g.elapsed_days)
        end

        return quantile(t, p)

    end

end

time_to_first_offer(res, 0.5, max_kdri = .8)
first_qualifying_offer_quantile(res, 0.25)

"""
    probability_qualifying_offer_by(df, t; max_kdri, n_simulations, censoring_time) -> Float64

Return the simulated probability of receiving an offer with KDRI at most
`max_kdri` by time `t`.
"""
function probability_qualifying_offer_by(
    df::AbstractDataFrame,
    t::Integer;
    max_kdri::Real=Inf,
    n_simulations::Integer=1000,
)::Float64

    check_df_columns(df, :simulation_id, :elapsed_days, :kdri)

    n_simulations > 0 ||
        throw(ArgumentError("`n_simulations` must be positive"))

    all(id -> id isa Integer && 1 ≤ id ≤ n_simulations, df.simulation_id) ||
        throw(ArgumentError("`:simulation_id` must lie in 1:$n_simulations"))


    has_qualifying_offer = falses(n_simulations)

    for simulation_df in groupby(df, :simulation_id)
        simulation_id = first(simulation_df.simulation_id)

        has_qualifying_offer[simulation_id] = any(
            row -> row.elapsed_days ≤ t && row.kdri ≤ max_kdri,
            eachrow(simulation_df),
        )
    end

    return mean(has_qualifying_offer)
end

@time probability_qualifying_offer_by(res, 50; max_kdri = .85, n_simulations=1000)

"""
    first_qualifying_offer_quantile_given_offer(df, p; max_kdri=Inf)

Return the `p`th quantile of time to first qualifying offer, conditional on
receiving at least one qualifying offer during the simulation.
"""
function first_qualifying_offer_quantile_given_offer(
    df::AbstractDataFrame,
    p::Real;
    max_kdri::Real=Inf,
)::Union{Float64,Nothing}
    check_df_columns(df, :simulation_id, :elapsed_days, :kdri)

    0 < p < 1 ||
        throw(ArgumentError("`p` must lie strictly between 0 and 1"))

    qualifying_df = filter(row -> row.kdri ≤ max_kdri, df)
    isempty(qualifying_df) && return nothing

    first_times = combine(
        groupby(qualifying_df, :simulation_id),
        :elapsed_days => minimum => :elapsed_days,
    )

    return quantile(first_times.elapsed_days, p)
end



id = 3748
# id = 10615

df_recipients = KidneyAllocation.load_recipient(recipient_filepath)
df_donors = KidneyAllocation.load_donor(donor_filepath)

df_recipient = filter(row -> row.CAN_ID == id, df_recipients)
obs_waiting_time_before_transplantation = Dates.days(KidneyAllocation.get_last_update(df_recipient)[1] - first(df_recipient.CAN_LISTING_DT))

df_donor = filter(row -> row.CAN_ID == id, df_donors)
obs_waiting_time_before_first_offer = Dates.days(minimum(df_donor.DON_DEATH_TM) - first(df_recipient.CAN_LISTING_DT))

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations/"
filepath = joinpath(sim_folder, "$id.csv")

res = CSV.read(filepath, DataFrame)

time_to_first_offer(res, .25)

df = combine(groupby(res, :simulation_id), :elapsed_days => minimum => :elapsed_days)

plot(df, y=:elapsed_days, Geom.boxplot)




paths = readdir(sim_folder; join = true)
files = filter(isfile, paths)

q1 = falses(length(files))
q2 = falses(length(files))
q3 = falses(length(files))

n_empty = 0

for (i, file) in enumerate(files)

    println(i)
    id = parse(Int64, first.(splitext.(basename.(file))))

    df_recipient = filter(row -> row.CAN_ID == id, df_recipients)
    df_donor = filter(row -> row.CAN_ID == id, df_donors)

    if isempty(df_recipient) || isempty(df_donor)
        n_empty +=1
        continue
    end

    obs_waiting_time_before_first_offer = Dates.days(minimum(df_donor.DON_DEATH_TM) - first(df_recipient.CAN_LISTING_DT))

    res = CSV.read(file, DataFrame)

    t1 = first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

end


n_empty

p₁ = count(q1)/(length(q1) - n_empty)
p₂ = count(q2)/(length(q2) - n_empty)
p₃ = count(q3)/(length(q3) - n_empty)




mean(q1)




import KidneyAllocation: retrieve_observed_waiting_list

t = Date(2012,1,1):Month(1):Date(2020,1,1)
n = Vector{Int64}(undef, length(t))

for (i,tᵢ) in enumerate(t)
    ids = retrieve_observed_waiting_list(recipient_filepath, tᵢ)
    n[i] = length(ids)
end

df = DataFrame(Date = t, Candidates = n)

plot(df, x=:Date, y=:Candidates)




df_donors


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



df_recipients

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




"""
    indices_in_registry(recipients, registry) -> Vector{Int}

Return the index of each recipient in `registry`, or `0` if not found.
"""
function indices_in_registry(
    recipients::AbstractVector{Recipient},
    registry::AbstractVector{Recipient},
)
    idx = Vector{Int}(undef, length(recipients))
    for (i, r) in enumerate(recipients)
        j = findfirst(==(r), registry)
        idx[i] = j === nothing ? 0 : j
    end
    return idx
end

@testset "indices_in_registry()" begin
    import KidneyAllocation.indices_in_registry
    # Registry (tiny and fakes)
    registry = [
        Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), A, 69, 2403, 7, 35, 4, 103, 10),
        Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),
    ]

    # Two first recipients in registry, but not the third
    recipients = vcat(registry[2], registry[1], Recipient(Date(1965, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),)

    @test indices_in_registry(recipients, registry) == [2, 1, 0]
end


"""
    active_recipient_ids(recipient_filepath, date) -> Vector{Int}

Return the `CAN_ID`s of recipients active on `date` after removing rows with
missing values in required fields.
"""
function active_recipient_ids(recipient_filepath::String, date::Date)
    df = load_recipient(recipient_filepath)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]
    dropmissing!(df, required)

    out = Int[]
    for g in groupby(df, :CAN_ID)
        a, d = recipient_arrival_departure(g)
        if a ≤ date < d
            push!(out, g.CAN_ID[1])
        end
    end

    return out
end

# Registry (tiny and fakes)
registry = [
    Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
    Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), A, 69, 2403, 7, 35, 4, 103, 10),
    Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),
]


"""
    recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids) -> Vector{Recipient}

Load recipient history and return one `Recipient` per `CAN_ID` in `can_ids`.

Throws an error if a requested `CAN_ID` is not found.
"""
function recipients_from_can_ids(
    recipient_filepath::String,
    cpra_filepath::String,
    can_ids::AbstractVector{<:Int},
)::Vector{Recipient}

    df_recipient = load_recipient(recipient_filepath)
    cpra_by_can_id = build_last_cpra_registry(cpra_filepath)

    # Group once
    G = groupby(df_recipient, :CAN_ID)

    # Build mapping CAN_ID -> SubDataFrame
    group_by_id = Dict(first(g.CAN_ID) => g for g in G)

    recipients = Vector{Recipient}(undef, length(can_ids))

    for (i, id) in enumerate(can_ids)

        haskey(group_by_id, id) ||
            throw(ArgumentError("CAN_ID $id not found in recipient file"))

        g = group_by_id[id]

        expiration_date = infer_recipient_expiration_date(g)

        sort!(g, :UPDATE_TM, rev=true)

        cpra = get(cpra_by_can_id, id, 0)

        recipients[i] = recipient_from_row(first(g), cpra, expiration_date)
    end

    return recipients
end


recipients_from_can_ids(recipient_filepath, cpra_filepath, [1, 3])


function get_active_recipients(recipient_filepath::String, cpra_filepath::String, date::Date)

    can_ids = active_recipient_ids(recipient_filepath::String, date::Date)

    recipients = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)

    return recipients
end

can_ids = active_recipient_ids(recipient_filepath, Date(2014, 1, 1,))
r = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)

r = get_active_recipients(recipient_filepath, cpra_filepath, Date(2014, 1, 1))


"""
    transplant_dates_by_recipient(df_donors) -> Dict{Int,Date}

Return a dictionary mapping each transplanted `CAN_ID` to its transplant date
(`DON_DEATH_TM`), based on accepted donor–recipient pairs.
"""
function transplant_dates_by_recipient(df_donors::AbstractDataFrame)

    df = unique(df_donors, [:CAN_ID, :DON_ID])
    filter!(row -> row.DECISION == "Acceptation", df)

    transplant_date_by_id = Dict{Int,Date}()

    for r in eachrow(df)
        transplant_date_by_id[r.CAN_ID] = r.DON_DEATH_TM
    end

    return transplant_date_by_id
end

df_donors = load_donor(donor_filepath)
transplant_dates_by_recipient(df_donors)






















## preprocessing raw data

using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Gadfly

using KidneyAllocation

import KidneyAllocation.parse_abo


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"


df_recipients = KidneyAllocation.load_recipient(recipient_filepath)
df_donors = KidneyAllocation.load_donor(donor_filepath)
df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])


CAN_ID = 14183

df_recipient = filter(row->row.CAN_ID==CAN_ID, df_recipients)
df_donor = filter(row->row.CAN_ID == CAN_ID, df_donors)
df_cpra_recipient = filter(row -> row.CAN_ID == CAN_ID, df_cpra)

if isempty(df_cpra_recipient)
    cpra = 0
else
    target = r.UPDATE_TM
    nearest_index = argmin(abs.(df_cpra_recipient.UPDATE_TM .- target))
    cpra = Int(df_cpra_recipient.CAN_CPRA[nearest_index])
end

r = first(df_recipient)
d = last(df_donor)

R = Recipient(r.CAN_BTH_DT,
    # d.DON_DEATH_TM,
    r.CAN_LISTING_DT,
    r.CAN_LISTING_DT,
    parse_abo(String(r.CAN_BLOOD)),
    r.CAN_A1,
    r.CAN_A2,
    r.CAN_B1,
    r.CAN_B2,
    r.CAN_DR1,
    r.CAN_DR2,
    cpra;
    expiration_date = nothing,
    active_waiting_proportion = 1.)

D = KidneyAllocation.donor_from_row(last(df_donor))

KidneyAllocation.score(D, R)





# Le patient CAN_ID = 14183 a un CPRA de 95, alors c'est pour ça que sa date de dialyse est inconnue.
# Le patient CAN_ID = 142 a reçu plusieurs offres de la liste 5, mais a accepté une offre de la liste 10. Il était compétitif jusqu'à ce moment.


# On va retirer les séjours multiples des patients, mais il faut les compter pour bien calculer la compétition


# On récupére les offres qui se sont produites entre 2016 et 2020
# filtered_donors = filter(row -> row.ATT_TYPE == 5, df_donors)
filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)

# On regroupe les patients qui ont reçu une de ces offres
G = groupby(filtered_donors, :CAN_ID)

g = G[1000]

df_recipient= filter(row -> row.CAN_ID == first(g.CAN_ID), df_recipients)

filter!()





filter(row->row.CAN_ID == 142, df_donors)




## Récupération des candidats qui ont reçu au moins une offre de la liste 5 entre 2016 et 2020

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)
G = groupby(filtered_donors, :DON_ID)

df_offers_recipients = similar(filtered_donors, 0)
given_kidneys_registry = Dict{Int, Int}()

for df_offer in groupby(filtered_donors, :DON_ID)
    df_list5 = filter(row -> row.ATT_TYPE == 5, df_offer)
    if !isempty(df_list5)
        append!(df_offers_recipients, df_list5)
        given_kidneys_registry[first(df_list5.DON_ID)] = count(skipmissing(df_list5.STATUS .== "TX"))
    end
end

df_offers_recipients

offered_recipients = unique(df_offers_recipients.CAN_ID)

# Pour chacun de ces patients, on va récupérer l'historique associé à ces offres et on va faire des vérifications de score

# CAN_ID = offered_recipients[29]
CAN_ID = 7181

df_offers_recipient = filter(row-> row.CAN_ID == CAN_ID, df_offers_recipients)
df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)

df_recipient_first_registration = filter_first_inscription(df_recipient)
r = first(df_recipient_first_registration)

df_cpra_recipient = filter(row -> row.CAN_ID == CAN_ID, df_cpra)

if isempty(df_cpra_recipient)
    cpra = 0
else
    target = r.UPDATE_TM
    nearest_index = argmin(abs.(df_cpra_recipient.UPDATE_TM .- target))
    cpra = Int(df_cpra_recipient.CAN_CPRA[nearest_index])
end

# cpra = 46

R = Recipient(r.CAN_BTH_DT,
    r.CAN_DIAL_DT,
    r.CAN_LISTING_DT,
    parse_abo(String(r.CAN_BLOOD)),
    r.CAN_A1,
    r.CAN_A2,
    r.CAN_B1,
    r.CAN_B2,
    r.CAN_DR1,
    r.CAN_DR2,
    cpra;
    expiration_date = nothing,
    active_waiting_proportion = 1.)

D = KidneyAllocation.donor_from_row(last(df_offers_recipient))

KidneyAllocation.score(D, R)





several_listing = Int64[]
for (i, CAN_ID) in enumerate(offered_recipients)
    println(i)
    df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)
    filter!(row -> row.UPDATE_TM ≥ Date(2016,1,1) && row.UPDATE_TM ≤ Date(2019,12,31), df_recipient)
    try 
        KidneyAllocation.check_df_at_most_one_exit_outcome(df_recipient)
    catch 
        push!(several_listing, CAN_ID)
    end
end

# Le patient 7387 revient sur la liste tout de suite après sa TX 
# Le patient 284 est bon exemple où toutes les données sont mélangées lors de la deuxième inscription
# Le patient 7181, on a que la transplantation est enregistrée le 2017-07-14 dans le fichiers des receveurs, mais le donneur est décédé le 2016-01-14

# Le patient 5373, il y a échec de la greffe, on devrait se concentrer que sur la première.


CAN_ID = several_listing[5]
df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)


import KidneyAllocation: check_df_columns, check_df_column_constant

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



filter_first_inscription(df_recipient)




## Récupération des candidats qui ont reçu au moins une offre de la liste 5 entre 2016 et 2020
# Reconstruction des recveurs par le fichier des donneurs

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)
G = groupby(filtered_donors, :DON_ID)

df_offers_recipients = similar(filtered_donors, 0)
given_kidneys_registry = Dict{Int, Int}()

for df_offer in groupby(filtered_donors, :DON_ID)
    df_list5 = filter(row -> row.ATT_TYPE == 5, df_offer)
    if !isempty(df_list5)
        append!(df_offers_recipients, df_list5)
        given_kidneys_registry[first(df_list5.DON_ID)] = count(skipmissing(df_list5.STATUS .== "TX"))
    end
end

offered_recipients = unique(df_offers_recipients.CAN_ID)

filtered_recipients = DataFrame(CAN_ID = Int[], CAN_BTH_DT = Date[], CAN_BLOOD = String[],
    CAN_A1 = Int[], CAN_A2 = Int[], CAN_B1 = Int[], CAN_B2 = Int[], CAN_DR1 = Int[], CAN_DR2 = Int[],
    CAN_LISTING_DT = DateTime[], CAN_DIAL_DT = Date[], 
    OUTCOME = String[], UPDATE_TM = DateTime[],
    CPRA = Int[], active_waiting_proportion = Float64[])

diff_score = DataFrame(CAN_ID = Int[], DON_ID = Int[], REGISTERED_SCORE = Float64[], COMPUTED_SCORE = Float64[])


for CAN_ID in offered_recipients
# CAN_ID = offered_recipients[2]

    df_offers_recipient = sort(filter(row-> row.CAN_ID == CAN_ID, df_offers_recipients), :DON_DEATH_TM, rev = true)

    # Récupération des informations sur le receveur
    df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)
    if ismissing(last(df_recipient.CAN_DIAL_DT))
        continue
    end
    df_recipient_first_registration = filter_first_waitlist_episode(df_recipient)
    p = KidneyAllocation.recipient_active_waiting_proportion(df_recipient_first_registration)

    # Récupération de la ligne correspondante à la première transplantation
    is_tx(outcome) = uppercase(strip(string(outcome))) == "TX"
    ind = findfirst(is_tx, skipmissing(df_offers_recipient.STATUS))

    # Si receveur non transplanté, on récupére la ligne de la dernière offre (1ère ligne dans le sorted dataframe)
    if isnothing(ind)
        ind = 1
    end

    # Date à laquelle on calcule le score
    target_date = df_offers_recipient.DON_DEATH_TM[ind]

    # Récupération du CPRA au moment du calcul du score
    df_cpra_recipient = sort(filter(row -> row.CAN_ID == CAN_ID, df_cpra), :UPDATE_TM, rev=true)
    latest_index = findfirst(t -> t ≤ tx_date, df_cpra_recipient.UPDATE_TM)
    cpra = isnothing(latest_index) ? 0 : Int(df_cpra_recipient.CAN_CPRA[latest_index])

    # Transformation en Recipient et Donor structures pour le calcul du score
    r = first(df_recipient_first_registration)

    R = Recipient(r.CAN_BTH_DT, r.CAN_DIAL_DT, r.CAN_LISTING_DT, parse_abo(String(r.CAN_BLOOD)),
            r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
            cpra; expiration_date = nothing, active_waiting_proportion = 1.)

    D = KidneyAllocation.donor_from_row(df_offers_recipient[ind, :])

    computed_score = KidneyAllocation.score(D, R)
    registered_score = df_offers_recipient.DON_CAN_SCORE[ind]

    DON_ID = df_offers_recipient.DON_ID[ind]

    if !isapprox(registered_score, computed_score, atol=.1)
        @warn "Registered score $registered_score and computed score $computed_score mismatches between CAN_ID = $CAN_ID and DON_ID = $DON_ID."
        push!(diff_score, [CAN_ID, DON_ID, registered_score, computed_score])
    else
        push!(filtered_recipients,
                [r.CAN_ID, r.CAN_BTH_DT, r.CAN_BLOOD,
                r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
                r.CAN_LISTING_DT, r.CAN_DIAL_DT,
                r.OUTCOME, tx_date,
                cpra, p])
    end
end




function load_data(recipient_filepath::AbstractString, cpra_filepath::AbstractString, donor_filepath::AbstractString)

    df_recipients = load_recipient(recipient_filepath)
    df_donors = load_donor(donor_filepath)
    df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])

    G = groupby(df_recipients, :CAN_ID)

    filtered_df_recipients = similar(df_recipients, 0)
    filtered_df_donors = similar(df_donors, 0)

    cpra = Vector{Int64}(undef, length(G))
    time_to_first_offer = Vector{Union{Int64,Missing}}(undef, length(G))
    time_to_transplant = Vector{Union{Int64,Missing}}(undef, length(G))
    active_waiting_proportion = Vector{Union{Float64,Missing}}(undef, length(G))

    exit_status = ["TX", "TX VIVANT", "X", "DCD"]

    for (i, df_recipient) in enumerate(G)

        r = last(sort(df_recipient, :UPDATE_TM))
        push!(filtered_df_recipients, r)

        CAN_ID = r.CAN_ID
        listing_date = r.CAN_LISTING_DT
        dialysis_date = r.CAN_DIAL_DT
        exit_date = uppercase(strip(r.OUTCOME)) ∈ exit_status ? r.UPDATE_TM : nothing


        df_cpra_recipient = filter(row -> row.CAN_ID == CAN_ID, df_cpra)

        if isempty(df_cpra_recipient)
            cpra[i] = 0
        else
            target = r.UPDATE_TM
            nearest_index = argmin(abs.(df_cpra_recipient.UPDATE_TM .- target))
            cpra[i] = Int(df_cpra_recipient.CAN_CPRA[nearest_index])
        end


        if uppercase(strip(r.OUTCOME)) == "TX"
            # time_to_transplant[i] = Dates.days(r.UPDATE_TM - r.CAN_LISTING_DT)
            time_to_transplant[i] = Dates.days(r.UPDATE_TM - dialysis_date)
        else
            time_to_transplant[i] = missing
        end

        active_waiting_proportion[i] = KidneyAllocation.recipient_active_waiting_proportion(df_recipient)


        df_offers = filter(row -> row.CAN_ID == CAN_ID && row.DON_DEATH_TM ≥ listing_date, df_donors)

        if !isnothing(exit_date)
            filter!(row -> row.DON_DEATH_TM ≤ exit_date, df_offers)
        end

        # Filtre les df_donors pour ne garder que les lignes qui correspondent à la bonne date de dialyse des receveurs
        append!(filtered_df_donors, df_offers)

        if !isempty(df_offers)
            # time_to_first_offer[i] = Dates.days(minimum(df_offers.DON_DEATH_TM) - listing_date)
            time_to_first_offer[i] = Dates.days(minimum(df_offers.DON_DEATH_TM) - dialysis_date)
        else
            time_to_first_offer[i] = missing
        end

    end

    filtered_df_recipients.CPRA = cpra
    filtered_df_recipients.active_waiting_proportion = active_waiting_proportion
    filtered_df_recipients.time_to_first_offer = time_to_first_offer
    filtered_df_recipients.time_to_transplant = time_to_transplant

    # Si time_to_first_offer == missing mais qu'on a time_to_transplant != missing, time_to_first_offer est fixé à time_to_transplant
    filtered_df_recipients[!, :time_to_first_offer] = coalesce.(filtered_df_recipients[!, :time_to_first_offer], filtered_df_recipients[!, :time_to_transplant])

    return filtered_df_recipients, filtered_df_donors

end