"""
    first_qualifying_offer_quantile_given_offer(df, p; max_kdri=Inf)

Return the `p`th quantile of time to first qualifying offer, conditional on
receiving at least one qualifying offer during the simulation. Simulation dataframe should be generated according to waiting_time_estimation.jl.
"""
function first_qualifying_offer_quantile_given_offer(
    df::AbstractDataFrame,
    p::Real;
    max_kdri::Real=Inf,
)::Union{Float64,Nothing}
    check_df_columns(df, :simulation_id, :elapsed_days_since_arrival, :kdri)

    0 < p < 1 ||
        throw(ArgumentError("`p` must lie strictly between 0 and 1"))

    qualifying_df = filter(row -> row.kdri ≤ max_kdri, df)
    isempty(qualifying_df) && return nothing

    first_times = combine(
        groupby(qualifying_df, :simulation_id),
        :elapsed_days_since_dialysis => minimum => :elapsed_days_since_dialysis,
    )

    return quantile(first_times.:elapsed_days_since_dialysis, p)
end

"""
    generate_arrivals(ids, arrival_rate; origin, nyears, rng) ->
        (sampled_ids, arrival_dates)

Generate a Poisson number of arrivals with mean `arrival_rate * nyears`.
Return identifiers sampled with replacement from `ids` and uniformly sampled
arrival dates over the simulation window.
"""
function generate_arrivals(
    ids::AbstractVector{<:Int},
    arrival_rate::Real;
    origin::Date,
    nyears::Int,
    rng::AbstractRNG=Random.default_rng(),
)
    arrival_rate ≥ 0 ||
        throw(ArgumentError("`arrival_rate` must be non-negative"))
    nyears ≥ 0 ||
        throw(ArgumentError("`nyears` must be non-negative"))

    n_arrivals = rand(rng, Poisson(arrival_rate * nyears))
    n_arrivals > 0 && isempty(ids) &&
        throw(ArgumentError("`ids` cannot be empty when generating arrivals"))

    sim_end = origin + Year(nyears)
    arrival_dates = sample_days(origin, sim_end, n_arrivals; rng=rng)
    sampled_ids = rand(rng, ids, n_arrivals)

    return sampled_ids, arrival_dates
end

"""
    generate_arrivals(registry::Dict{Int,Recipient}, arrival_rate; ...) ->
        Vector{Recipient}

Generate recipient arrivals by sampling registry templates and shifting their
timelines to simulated arrival dates.
"""
function generate_arrivals(
    registry::Dict{Int,Recipient},
    arrival_rate::Real;
    origin::Date,
    nyears::Int,
    rng::AbstractRNG=Random.default_rng(),
)::Vector{Recipient}

    ids = collect(keys(registry))
    sampled_ids, sampled_arrivals = generate_arrivals(
        ids, arrival_rate; origin=origin, nyears=nyears, rng=rng,
    )

    return reconstruct_recipients(registry, sampled_ids, sampled_arrivals)
end

"""
    generate_arrivals(registry::Dict{Int,Donor}, arrival_rate; ...) -> Vector{Donor}

Generate donor arrivals by sampling registry templates and shifting their
timelines to simulated arrival dates.
"""
function generate_arrivals(
    registry::Dict{Int,Donor},
    kidney_by_id::Dict{Int,Int},
    arrival_rate::Real;
    origin::Date,
    nyears::Int,
    rng::AbstractRNG=Random.default_rng(),
)::Vector{Donor}

    ids = collect(keys(registry))
    sampled_ids, sampled_arrivals = generate_arrivals(
        ids, arrival_rate; origin=origin, nyears=nyears, rng=rng,
    )

    multiple_sampled_ids = Int64[]
    multiple_arrival_dates = Date[]

    for (i, (id, arrival)) in enumerate(zip(sampled_ids, sampled_arrivals))
        for j = 1:kidney_by_id[id]
            push!(multiple_sampled_ids, id)
            push!(multiple_arrival_dates, arrival)
        end
    end
    return reconstruct_donors(registry, multiple_sampled_ids, multiple_arrival_dates)
end

"""
    probability_qualifying_offer_by(df, t; max_kdri, n_simulations, censoring_time) -> Float64

Return the simulated probability of receiving an offer with KDRI at most
`max_kdri` by time `t`. Simulation dataframe should be generated according to waiting_time_estimation.jl.
"""
function probability_qualifying_offer_by(
    df::AbstractDataFrame,
    t::Integer;
    max_kdri::Real=Inf,
    n_simulations::Integer=1000,
)::Float64

    check_df_columns(df, :simulation_id, :elapsed_days_since_dialysis, :kdri)

    n_simulations > 0 ||
        throw(ArgumentError("`n_simulations` must be positive"))

    all(id -> id isa Integer && 1 ≤ id ≤ n_simulations, df.simulation_id) ||
        throw(ArgumentError("`:simulation_id` must lie in 1:$n_simulations"))


    has_qualifying_offer = falses(n_simulations)

    for simulation_df in groupby(df, :simulation_id)
        simulation_id = first(simulation_df.simulation_id)

        has_qualifying_offer[simulation_id] = any(
            row -> row.elapsed_days_since_dialysis ≤ t && row.kdri ≤ max_kdri,
            eachrow(simulation_df),
        )
    end

    return mean(has_qualifying_offer)
end

"""
    reconstruct_recipients(recipient_registry, ids, arrival_dates) -> Vector{Recipient}

Return the recipients identified by `ids` in the `recipient_registry`, with each timeline shifted so that
its arrival date equals the corresponding value in `arrival_dates`.
"""
function reconstruct_recipients(
    recipient_registry::Dict{Int,Recipient},
    ids::AbstractVector{<:Integer},
    arrival_dates::AbstractVector{Date},
)::Vector{Recipient}

    length(ids) == length(arrival_dates) ||
        throw(ArgumentError("`ids` and `arrival_dates` must have the same length"))

    reconstructed = Vector{Recipient}(undef, length(ids))

    for (i, (recipient_id, arrival_date)) in enumerate(zip(ids, arrival_dates))
        reconstructed[i] = shift_recipient_timeline(recipient_registry[recipient_id], arrival_date)
    end

    return reconstructed
end

"""
    reconstruct_donors(donors_registry, ids, arrival_dates) -> Vector{Donor}

Return the donors identified by `ids` in the `donor_registry`, with arrival dates shifted to
`arrival_dates`.
"""
function reconstruct_donors(
    donor_registry::Dict{Int,Donor},
    ids::AbstractVector{<:Integer},
    arrival_dates::AbstractVector{Date},
)::Vector{Donor}

    length(ids) == length(arrival_dates) ||
        throw(ArgumentError("`ids` and `arrival_dates` must have the same length"))

    reconstructed = Vector{Donor}(undef, length(ids))

    for (i, (id, arrival_date)) in enumerate(zip(ids, arrival_dates))
        reconstructed[i] = set_donor_arrival(donor_registry[id], arrival_date)
    end

    return reconstructed
end


"""
    simulate_recipient_offers(...)

Run `nsim` kidney-allocation simulations over `nyears` years and return, for
each simulation, the donors whose offers reach `recipient`.

The focal recipient is evaluated against the observed waiting list at their
arrival date plus simulated recipient arrivals, and is assumed to refuse every
offer.
"""
function simulate_recipient_offers(recipient::Recipient,
    df_recipients::AbstractDataFrame,
    recipient_registry::Dict{Int, Recipient},
    donor_registry::Dict{Int, Donor},
    given_kidneys_registry::Dict{Int, Int},
    dm::AbstractDecisionModel;
    nyears::Real,
    nsim::Int,
    recipient_arrival_rate::Real,
    donor_arrival_rate::Real,
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
    )

    # Get the initial list
    initial_recipient_ids = KidneyAllocation.retrieve_observed_waiting_list(df_recipients, recipient.arrival - Day(1))

    initial_recipients = Recipient[]

    for id in initial_recipient_ids
        if id in keys(recipient_registry)
            push!(initial_recipients, recipient_registry[id])
        end
    end

    
    offers = Vector{Vector{Donor}}(undef, nsim)

    for i in 1:nsim

        new_recipients = KidneyAllocation.generate_arrivals(recipient_registry, recipient_arrival_rate; origin=recipient.arrival, nyears=nyears, rng=rng)

        waiting_recipients = vcat(initial_recipients, new_recipients)


        ## Generate donor arrivals for the next nyears

        donors = KidneyAllocation.generate_arrivals(donor_registry, given_kidneys_registry, donor_arrival_rate, origin = recipient.arrival,nyears=nyears, rng=rng)


        ## Retrieve all the donors offered to recipient

        offered_donors = KidneyAllocation.get_recipient_offers(recipient, donors, waiting_recipients, dm)

        offers[i] = unique(offered_donors) 

    end

    return offers

end


"""
    offers_to_dataframes(focal_recipient, offers_by_simulation) ->
        Tuple{DataFrame,DataFrame}

Convert simulated offers for `focal_recipient` into a long-format offer table
and a simulation-summary table.

The offer table has one row per offer. The summary table has one row per
simulation, including simulations with no offers.
"""
function offers_to_dataframes(
    focal_recipient::Recipient,
    offers_by_simulation::Vector{Vector{Donor}},
    dm::AbstractDecisionModel,
)
    offers_df = DataFrame(
        simulation_id=Int[],
        offer_number=Int[],
        donor_arrival=Date[],
        elapsed_days_since_arrival=Int[],
        elapsed_days_since_dialysis=Int[], 
        kdri=Float64[],
        acceptance_probability = Float64[],
    )

    for (simulation_id, offered_donors) in enumerate(offers_by_simulation)
        for (offer_number, donor) in enumerate(offered_donors)
            donor_date = Date(donor.arrival)

            p = acceptance_probability(dm, focal_recipient, donor)

            push!(offers_df, (
                simulation_id=simulation_id,
                offer_number=offer_number,
                donor_arrival=donor_date,
                elapsed_days_since_arrival = Dates.value(donor_date - focal_recipient.arrival),
                elapsed_days_since_dialysis = Dates.value(donor_date - focal_recipient.dialysis),
                kdri=KidneyAllocation.get_kdri(donor),
                acceptance_probability = p,
            ))
        end
    end

    return offers_df
end

