
# """
#     sample_arrival_dates(origin, sim_end, n) -> Vector{Date}

# Sample `n` dates uniformly in [`origin`, `sim_end`].
# """
# sample_arrival_dates(origin::Date, sim_end::Date, n::Integer) =
#     KidneyAllocation.sample_days(origin, sim_end, n)

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
    origin,
    nyears,
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
        ids, arrival_rate; origin, nyears, rng,
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
        ids, arrival_rate; origin, nyears, rng,
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
    simulate_initial_state_indexed(donors, recipients, dm; 
        start_date, nyears, donor_rate, recipient_rate, origin_date, rng
    ) -> (final_indices, shifted_arrival_dates)

Simulate arrivals and allocations over the time horizon and return the registry
indices and shifted arrival dates of recipients still active at the end.

# Keyword Arguments
- `start_date::Date`: Start of the simulation window.
- `nyears::Int`: Length of the simulation horizon in years.
- `donor_rate::Real`: Mean annual donor arrival rate.
- `recipient_rate::Real`: Mean annual recipient arrival rate.
- `origin_date::Date`: Reference date used to shift arrival times.
- `rng::AbstractRNG`: Random number generator.
"""
function simulate_initial_state_indexed(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel;
    start_date::Date=Date(2014, 1, 1),
    nyears::Int=10,
    donor_rate::Real=148.0,
    recipient_rate::Real=272.83,
    origin_date::Date=Date(2000, 1, 1),
    rng::AbstractRNG=Random.default_rng(),
)

    simulation_end = start_date + Year(nyears)

    # Recipients active at start_date: registry indices + their arrival dates
    active_at_start_mask = is_active.(recipients, start_date)
    waiting_registry_indices = findall(active_at_start_mask)
    waiting_arrival_dates = KidneyAllocation.get_arrival.(recipients[waiting_registry_indices])

    # New recipient arrivals: registry indices + simulated arrival dates
    sampled_recipient_indices, sampled_recipient_arrival_dates =
        generate_arrivals(eachindex(recipients), recipient_rate;
            origin=start_date, nyears=nyears, rng=rng)

    append!(waiting_registry_indices, sampled_recipient_indices)
    append!(waiting_arrival_dates, sampled_recipient_arrival_dates)

    # New donor arrivals (used only internally for allocation)
    sampled_donor_indices, sampled_donor_arrival_dates =
        generate_arrivals(eachindex(donors), donor_rate;
            origin=start_date, nyears=nyears, rng=rng)

    # Reconstruct temporary objects for allocation only
    waiting_recipients =
        reconstruct_recipients(recipients, waiting_registry_indices, waiting_arrival_dates)

    arriving_donors =
        reconstruct_donors(donors, sampled_donor_indices, sampled_donor_arrival_dates)

    # Allocate donors (positions refer to waiting_recipients)
    allocated_positions = allocate(arriving_donors, waiting_recipients, dm)

    # Filter non-attributed organ (i.e. ind ==0)
    filter!(>(0), allocated_positions)
    sort!(allocated_positions)

    # Remove transplanted recipients from the indexed representation
    deleteat!(waiting_registry_indices, allocated_positions)
    deleteat!(waiting_arrival_dates, allocated_positions)
    deleteat!(waiting_recipients, allocated_positions)

    # Keep recipients active at end of simulation
    active_at_end_mask = is_active.(waiting_recipients, simulation_end)

    final_recipient_indices = waiting_registry_indices[active_at_end_mask]
    final_arrival_dates = waiting_arrival_dates[active_at_end_mask]

    # Shift arrivals so that simulation_end maps to origin_date
    time_waited = simulation_end .- final_arrival_dates
    shifted_arrival_dates = origin_date .- time_waited

    return final_recipient_indices, shifted_arrival_dates
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
    recipient_filepath::AbstractString,
    recipient_registry::Dict{Int, Recipient},
    donor_filepath::AbstractString,
    donor_registry::Dict{Int, Donor},
    dm::AbstractDecisionModel,
    nyears::Real,
    nsim::Int;
    recipient_arrival_rate::Real,
    donor_arrival_rate::Real,
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
    )

    # Get the initial list
    initial_recipient_ids = KidneyAllocation.retrieve_observed_waiting_list(recipient_filepath, recipient.arrival - Day(1))

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

        kidney_by_id = KidneyAllocation.kidneys_given_by_donor(KidneyAllocation.load_donor(donor_filepath))

        donors = KidneyAllocation.generate_arrivals(donor_registry, kidney_by_id, donor_arrival_rate, origin = recipient.arrival,nyears=nyears, rng=rng)


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
)
    offers_df = DataFrame(
        simulation_id=Int[],
        offer_number=Int[],
        donor_arrival=Date[],
        elapsed_days=Int[],
        kdri=Float64[],
    )

    for (simulation_id, offered_donors) in enumerate(offers_by_simulation)
        for (offer_number, donor) in enumerate(offered_donors)
            donor_date = Date(donor.arrival)

            push!(offers_df, (
                simulation_id=simulation_id,
                offer_number=offer_number,
                donor_arrival=donor_date,
                elapsed_days=Dates.value(donor_date - focal_recipient.arrival),
                kdri=KidneyAllocation.get_kdri(donor), # or donor.kdri, depending on your API
            ))
        end
    end

    return offers_df
end