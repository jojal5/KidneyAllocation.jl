
"""
    sample_arrival_dates(origin, sim_end, n) -> Vector{Date}

Sample `n` dates uniformly in [`origin`, `sim_end`].
"""
sample_arrival_dates(origin::Date, sim_end::Date, n::Integer) =
    KidneyAllocation.sample_days(origin, sim_end, n)

"""
    generate_arrivals(indices, arrival_rate; origin, nyears, rng) -> (sampled_indices, arrival_dates)

Sample a Poisson number of arrivals and return resampled `indices` with uniformly sampled arrival dates over the simulation window.
"""
function generate_arrivals(indices::AbstractVector{<:Int}, arrival_rate::Real;
    origin::Date = Date(2000, 1, 1),
    nyears::Int = 10,
    rng::AbstractRNG = Random.default_rng())

    sim_end = origin + Year(nyears)

    n_arrivals = rand(rng, Poisson(arrival_rate * nyears))

    arrival_dates = KidneyAllocation.sample_arrival_dates(origin, sim_end, n_arrivals)
    sampled_indices = rand(rng, indices, n_arrivals)

    return sampled_indices, arrival_dates

end

"""
    reconstruct_recipients(recipients, indices, arrival_dates) -> Vector{Recipient}

Return `recipients[indices[i]]` with timelines shifted to `arrival_dates[i]`.
"""
function reconstruct_recipients(
    recipients::AbstractVector{Recipient},
    indices::AbstractVector{<:Integer},
    arrival_dates::AbstractVector{<:Date},
)::Vector{Recipient}

    length(indices) == length(arrival_dates) ||
        throw(ArgumentError("`indices` and `arrival_dates` must have the same length"))

    n = length(indices)
    reconstructed = Vector{Recipient}(undef, n)

    @inbounds for i in 1:n
        reconstructed[i] = shift_recipient_timeline(recipients[indices[i]], arrival_dates[i])
    end

    return reconstructed
end

"""
    reconstruct_donors(donors, indices, arrival_dates) -> Vector{Donor}

Return `donors[indices[i]]` with arrival set to `arrival_dates[i]`.
"""
function reconstruct_donors(
    donors::AbstractVector{Donor},
    indices::AbstractVector{<:Integer},
    arrival_dates::AbstractVector{<:Date},
)::Vector{Donor}

    length(indices) == length(arrival_dates) ||
        throw(ArgumentError("`indices` and `arrival_dates` must have the same length"))

    n = length(indices)
    reconstructed = Vector{Donor}(undef, n)

    @inbounds for i in 1:n
        reconstructed[i] = set_donor_arrival(donors[indices[i]], arrival_dates[i])
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
    start_date::Date = Date(2014, 1, 1),
    nyears::Int = 10,
    donor_rate::Real = 242.0,
    recipient_rate::Real = 272.83,
    origin_date::Date = Date(2000, 1, 1),
    rng::AbstractRNG = Random.default_rng(),
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
