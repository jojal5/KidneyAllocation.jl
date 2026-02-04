
"""
    sample_arrival_dates(origin::Date, sim_end::Date, n::Integer) -> Vector{Date}

Sample `n` arrival dates uniformly over [`origin`, `sim_end`].

This is a small convenience wrapper around `KidneyAllocation.sample_days`.
"""
sample_arrival_dates(origin::Date, sim_end::Date, n::Integer) =
    KidneyAllocation.sample_days(origin, sim_end, n)

"""
    generate_recipient_arrivals(
        recipients::Vector{Recipient};
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        arrival_rate::Real = 272.83,
        rng::AbstractRNG = Random.default_rng()
    ) -> Vector{Recipient}

Generate a synthetic sequence of recipient arrivals over a given time horizon.

Arrivals are modeled as a Poisson process with mean `arrival_rate * nyears`.
Recipient profiles are resampled from `recipients` and their timelines are
shifted to match simulated arrival dates within [`origin`, `origin + Year(nyears)`].

# Keyword arguments
- `rng`: Random number generator used for reproducible simulations.

# See also
- [`generate_donor_arrivals`](@ref)
- [`generate_pseudo_history`](@ref)
"""
function generate_recipient_arrivals(
    recipients::Vector{Recipient};
    origin::Date = Date(2000, 1, 1),
    nyears::Int = 10,
    arrival_rate::Real = 272.83,
    rng::AbstractRNG = Random.default_rng(),
)::Vector{Recipient}

    sim_end = origin + Year(nyears)

    n_arrivals = rand(rng, Poisson(arrival_rate * nyears))
    arrival_dates = sample_arrival_dates(origin, sim_end, n_arrivals)

    sampled = rand(rng, recipients, n_arrivals)
    return KidneyAllocation.shift_recipient_timeline.(sampled, arrival_dates)
end

"""
    generate_donor_arrivals(
        donors::Vector{Donor};
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        arrival_rate::Real = 242.0,
        rng::AbstractRNG = Random.default_rng()
    ) -> Vector{Donor}

Generate a synthetic sequence of donor arrivals over a given time horizon.

Arrivals are modeled as a Poisson process with mean `arrival_rate * nyears`.
Donor profiles are resampled from `donors` and their arrival dates are updated
within [`origin`, `origin + Year(nyears)`].

# Keyword arguments
- `rng`: Random number generator used for reproducible simulations.

# See also
- [`generate_recipient_arrivals`](@ref)
- [`generate_pseudo_history`](@ref)
"""
function generate_donor_arrivals(
    donors::Vector{Donor};
    origin::Date = Date(2000, 1, 1),
    nyears::Int = 10,
    arrival_rate::Real = 242.0,
    rng::AbstractRNG = Random.default_rng(),
)::Vector{Donor}

    sim_end = origin + Year(nyears)

    n_arrivals = rand(rng, Poisson(arrival_rate * nyears))
    arrival_dates = sample_arrival_dates(origin, sim_end, n_arrivals)

    sampled = rand(rng, donors, n_arrivals)
    return set_donor_arrival.(sampled, arrival_dates)
end

"""
    simulate_initial_recipient_list(
        recipients::Vector{Recipient},
        donors::Vector{Donor},
        fm,
        u;
        origin::Date = Date(2014,1,1),
        nyears::Int = 10,
        donor_rate::Real = 242.0,
        recipient_rate::Real = 272.83,
        reference_date::Date = Date(2000,1,1),
        rng::AbstractRNG = Random.default_rng()
    ) -> Vector{Recipient}

Simulate the evolution of a kidney transplant waiting list over a fixed horizon
and return an updated "initial" recipient population.

Workflow:
1. Keep recipients active at `origin`.
2. Simulate recipient arrivals and donor arrivals over the next `nyears` years.
3. Allocate donors sequentially using `allocate(waiting_recipients, new_donors, fm, u)`.
4. Remove transplanted recipients and then remove recipients inactive at the end.
5. Shift the timelines of remaining recipients so that the end of the simulation
   period maps to `reference_date`.

# Keyword arguments
- `reference_date`: Date used as the new time origin after shifting timelines.
- `rng`: Random number generator used for reproducibility.

# Notes
- This function uses resampling from the provided reference pools of recipients
  and donors to generate synthetic arrivals.
- If you run large Monte Carlo experiments, consider index-based and buffer-reuse
  variants to reduce allocations.

# See also
- [`generate_pseudo_history`](@ref)
- [`allocate`](@ref)
"""
function simulate_initial_recipient_list(
    recipients::Vector{Recipient},
    donors::Vector{Donor},
    fm,
    u;
    origin::Date = Date(2014, 1, 1),
    nyears::Int = 10,
    donor_rate::Real = 242.0,
    recipient_rate::Real = 272.83,
    reference_date::Date = Date(2000, 1, 1),
    rng::AbstractRNG = Random.default_rng(),
)::Vector{Recipient}

    sim_end = origin + Year(nyears)

    # Waiting list at origin
    is_waiting = is_active.(recipients, origin)
    waiting_recipients = recipients[is_waiting]

    # Add new recipients
    new_recipients = generate_recipient_arrivals(
        recipients;
        origin=origin,
        nyears=nyears,
        arrival_rate=recipient_rate,
        rng=rng,
    )
    append!(waiting_recipients, new_recipients)

    # Generate donors
    new_donors = generate_donor_arrivals(
        donors;
        origin=origin,
        nyears=nyears,
        arrival_rate=donor_rate,
        rng=rng,
    )

    # Allocate donors (uses package allocation policy)
    allocation_indices = allocate(waiting_recipients, new_donors, fm, u)

    # Remove transplanted recipients (indices refer to `waiting_recipients`)
    filter!(>(0), allocation_indices)
    sort!(allocation_indices)
    deleteat!(waiting_recipients, allocation_indices)

    # Keep only active at end
    is_active_at_end = is_active.(waiting_recipients, sim_end)
    active_recipients = waiting_recipients[is_active_at_end]

    # Shift timelines so sim_end becomes reference_date
    waited_times = sim_end .- KidneyAllocation.get_arrival.(active_recipients)
    shifted_arrivals = reference_date .- waited_times
    return shift_recipient_timeline.(active_recipients, shifted_arrivals)
end

"""
    generate_pseudo_history(
        recipients::Vector{Recipient},
        donors::Vector{Donor},
        fm,
        u;
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        donor_rate::Real = 242.0,
        recipient_rate::Real = 272.83,
        rng::AbstractRNG = Random.default_rng()
    ) -> (Vector{Recipient}, Vector{Donor})

Generate a synthetic ("pseudo") history consisting of:
- an initial waiting list (synthetic baseline recipients), and
- a donor stream over the simulation period.

The initial waiting list is produced by [`simulate_initial_recipient_list`](@ref).
New recipients and donors are generated over the interval
[`origin`, `origin + Year(nyears)`] using Poisson arrival models and resampling.

# Returns
- `(waiting_recipients, new_donors)` where:
  - `waiting_recipients` includes the synthetic initial recipients plus new
    recipient arrivals over the period,
  - `new_donors` are donors with simulated arrival dates over the period.

# Keyword arguments
- `rng`: Random number generator for reproducibility.

# See also
- [`simulate_initial_recipient_list`](@ref)
- [`generate_recipient_arrivals`](@ref)
- [`generate_donor_arrivals`](@ref)
"""
function generate_pseudo_history(
    recipients::Vector{Recipient},
    donors::Vector{Donor},
    fm,
    u;
    origin::Date = Date(2000, 1, 1),
    nyears::Int = 10,
    donor_rate::Real = 242.0,
    recipient_rate::Real = 272.83,
    rng::AbstractRNG = Random.default_rng(),
)
    initial_recipients = simulate_initial_recipient_list(
        recipients, donors, fm, u;
        rng=rng,
    )

    new_recipients = generate_recipient_arrivals(
        recipients;
        origin=origin,
        nyears=nyears,
        arrival_rate=recipient_rate,
        rng=rng,
    )

    waiting_recipients = vcat(initial_recipients, new_recipients)

    new_donors = generate_donor_arrivals(
        donors;
        origin=origin,
        nyears=nyears,
        arrival_rate=donor_rate,
        rng=rng,
    )

    return waiting_recipients, new_donors
end