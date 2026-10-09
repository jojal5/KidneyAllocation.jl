using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Random, Statistics, Test

using Gadfly

using KidneyAllocation


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"


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






