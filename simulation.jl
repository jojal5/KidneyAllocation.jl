using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, load_recipient, is_active, is_expired, is_abo_compatible
import KidneyAllocation: load_donor, build_donor_registry
import KidneyAllocation: shift_recipient_timeline, set_donor_arrival
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision
import KidneyAllocation: score, years_between, fractionalyears_between
import KidneyAllocation: allocate_one_donor, allocate


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

# Estimate the recipient arrival rate
df = load_recipient(recipient_filepath)
df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)
λᵣ = length(unique(df2.CAN_ID)) / 6

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

donors = build_donor_registry(donor_filepath)

# Estimate the donor arrival rate
df = load_donor(donor_filepath)
df2 = filter(row -> 2014 ≤ year(row.DON_DEATH_TM) < 2020, df)
filter!(row -> row.DECISION == "Acceptation", df2)
λₒ = length(unique(df2.CAN_ID)) / 6


# Fit decision model
data = retrieve_decision_data(donor_filepath, recipient_filepath)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

fm = glm(model, data, Bernoulli(), LogitLink())

u = fit_decision_threshold(fm)
# ------------------------------------------------------------------------------------

# Retrieve the waiting recipients at January 1st, 2014
ind_active = is_active.(recipients, Date(2013, 12, 31))
waiting_recipients = recipients[ind_active]

# Generate new recipients for the next 10 years
nᵣ = rand(Poisson(λᵣ * 10))                                                             # Number of recipients
tᵣ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nᵣ)             # Arrival dates
sampled_recipients = rand(recipients, nᵣ)                                               # Sampled recipients
new_recipients = KidneyAllocation.shift_recipient_timeline.(sampled_recipients, tᵣ)     # Adjust the sampled recipient timelines

append!(waiting_recipients, new_recipients)                                             # Add the new recipients to the original list


# Generate new donors for the next 10 years
nₒ = rand(Poisson(λₒ * 10))                                                             # Number of recipients
tₒ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nₒ)             # Arrival dates
sampled_donors = rand(donors, nₒ)                                                       # Sampled recipients
new_donors = set_donor_arrival.(sampled_donors, tₒ)                                     # Adjust the arrival dates

# KidneyAllocation.get_arrival.(new_donors)


donor = new_donors[1]
arrival = donor.arrival
eligible_index = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_recipients = waiting_recipients[eligible_index]
chosen_recipient_idx = allocate_one_donor(eligible_recipients, donor, fm, u)
chosen_recipient = eligible_recipients[chosen_recipient_idx]

@time allocate_one_donor(eligible_recipients, donor, fm, u)

# Sanity checks
score.(donor, chosen_recipient)
get_decision(donor, chosen_recipient, fm, u)


@time ind = allocate(waiting_recipients, new_donors, fm, u)

# Sanity checks
score(new_donors[100], waiting_recipients[ind[100]])
get_decision(new_donors[100], waiting_recipients[ind[100]], fm, u)


@time ind = allocate(waiting_recipients, new_donors, fm, u, until = 1)

findlast(ind .!= 0)


@time r = KidneyAllocation.simulate_initial_recipient_list(recipients, donors, fm, u)

"""
    generate_pseudo_history(
        recipients::Vector{Recipient},
        donors::Vector{Donor},
        fm,
        u;
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        donor_rate::Real = 242.0,
        recipient_rate::Real = 272.83
    ) -> (Vector{Recipient}, Vector{Donor})

Generate a synthetic ("pseudo") transplant history over a given time horizon.

This function constructs a simulated waiting list and donor stream by:

1. Generating an initial recipient population using
   [`simulate_initial_recipient_list`](@ref).
2. Simulating new recipient arrivals over the period
   [`origin`, `origin + Year(nyears)`] using a Poisson process with rate
   `recipient_rate`.
3. Simulating new donor arrivals over the same period using a Poisson process
   with rate `donor_rate`.
4. Resampling recipient and donor profiles and shifting their timelines to
   match simulated arrival dates.

The resulting output represents a synthetic historical dataset that can be
used as input for downstream allocation simulations.

# Arguments
- `recipients::Vector{Recipient}`: Reference pool of recipient profiles.
- `donors::Vector{Donor}`: Reference pool of donor profiles.
- `fm`: Fitted statistical model used in allocation decisions.
- `u`: Decision parameter (e.g., acceptance threshold).
- `origin::Date`: Starting date of the simulated history.
- `nyears::Int`: Length of the simulated period in years.
- `donor_rate::Real`: Mean annual arrival rate of donors.
- `recipient_rate::Real`: Mean annual arrival rate of recipients.

# Returns
- `(Vector{Recipient}, Vector{Donor})`:
    - `waiting_recipients`: Synthetic waiting list at the end of the
      initialization phase.
    - `new_donors`: Synthetic donor arrivals over the simulation period.

# Modeling Assumptions
- Recipient and donor arrivals follow independent Poisson processes.
- New individuals are generated by resampling existing profiles and shifting
  their timelines.
- Initial recipients reflect the output of
  `simulate_initial_recipient_list`.

# Performance Notes
- This function allocates new vectors when concatenating recipient lists and
  generating arrivals.
- For large-scale Monte Carlo studies, repeated calls may benefit from
  buffer reuse and in-place variants.

# Reproducibility
- Results depend on the global random number generator state.
  For reproducible simulations, set a seed before calling this function:

      Random.seed!(1234)

# See also
- [`simulate_initial_recipient_list`](@ref)
- [`allocate`](@ref)
- [`allocate_one_donor`](@ref)
- [`get_decision`](@ref)
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
)
    initial_recipients = KidneyAllocation.simulate_initial_recipient_list(recipients, donors, fm, u)

    sim_end = origin + Year(nyears)

    # Generate new recipients over the simulation period
    n_rec = rand(Poisson(recipient_rate * nyears))
    rec_arrivals = KidneyAllocation.sample_days(origin, sim_end, n_rec)
    sampled_recipients = rand(recipients, n_rec)
    new_recipients = KidneyAllocation.shift_recipient_timeline.(sampled_recipients, rec_arrivals)
    waiting_recipients = vcat(initial_recipients, new_recipients)

    # Generate new donors over the simulation period
    n_don = rand(Poisson(donor_rate * nyears))
    don_arrivals = KidneyAllocation.sample_days(origin, sim_end, n_don)
    sampled_donors = rand(donors, n_don)
    new_donors = set_donor_arrival.(sampled_donors, don_arrivals)

    return waiting_recipients, new_donors

end

r,d = generate_pseudo_history(recipients, donors, fm, u)


"""
    generate_recipient_arrivals(
        recipients::Vector{Recipient};
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        arrival_rate::Real = 272.83
    ) -> Vector{Recipient}

Generate a synthetic sequence of recipient arrivals over a given time horizon.

This function simulates new recipient arrivals between `origin` and
`origin + Year(nyears)` by:

1. Drawing the total number of arrivals from a Poisson distribution with
   mean `arrival_rate * nyears`.
2. Sampling recipient profiles from the reference population `recipients`.
3. Assigning each sampled recipient a random arrival date within the
   simulation period.
4. Shifting recipient timelines to match the simulated arrival dates.

The output can be used to augment an existing waiting list in allocation
simulations.

# Arguments
- `recipients::Vector{Recipient}`: Reference pool of recipient profiles used
  for resampling.
- `origin::Date`: Starting date of the simulation period.
- `nyears::Int`: Length of the simulation horizon in years.
- `arrival_rate::Real`: Mean annual arrival rate of recipients.

# Returns
- `Vector{Recipient}`: Vector of newly generated recipients with updated
  arrival times.

# Modeling Assumptions
- Recipient arrivals follow a Poisson process.
- Arrival times are uniformly distributed over the simulation period.
- New recipients are generated by resampling existing profiles and shifting
  their timelines.

# Performance Notes
- Broadcasting in `shift_recipient_timeline.(...)` allocates temporary arrays.
- For large-scale simulations, consider buffer-reusing variants.

# See also
- [`simulate_initial_recipient_list`](@ref)
- [`generate_pseudo_history`](@ref)
- [`allocate`](@ref)
"""
function generate_recipient_arrivals(recipients::Vector{Recipient}; origin::Date=Date(2000,1,1), nyears::Int=10, arrival_rate::Real=272.83)

    sim_end = origin + Year(nyears)

    n_arrivals = rand(Poisson(arrival_rate * nyears))
    arrival_dates = KidneyAllocation.sample_days(origin, sim_end, n_arrivals)
    sampled_recipients = rand(recipients, n_arrivals)
    new_recipients = KidneyAllocation.shift_recipient_timeline.(sampled_recipients, arrival_dates)
    
    return new_recipients

end

generate_recipient_arrivals(recipients)


"""
    generate_donor_arrivals(
        donors::Vector{Donor};
        origin::Date = Date(2000,1,1),
        nyears::Int = 10,
        arrival_rate::Real = 242.0
    ) -> Vector{Donor}

Generate a synthetic sequence of donor arrivals over a given time horizon.

This function simulates new donor arrivals between `origin` and
`origin + Year(nyears)` by:

1. Drawing the total number of arrivals from a Poisson distribution with
   mean `arrival_rate * nyears`.
2. Sampling donor profiles from the reference population `donors`.
3. Assigning each sampled donor a random arrival date within the simulation
   period.
4. Updating donor timelines to match the simulated arrival dates.

The output can be used as input for kidney allocation simulations.

# Arguments
- `donors::Vector{Donor}`: Reference pool of donor profiles used for resampling.
- `origin::Date`: Starting date of the simulation period.
- `nyears::Int`: Length of the simulation horizon in years.
- `arrival_rate::Real`: Mean annual arrival rate of donors.

# Returns
- `Vector{Donor}`: Vector of newly generated donors with updated arrival times.

# Modeling Assumptions
- Donor arrivals follow a Poisson process.
- Arrival times are uniformly distributed over the simulation period.
- New donors are generated by resampling existing profiles and shifting their
  timelines.

# Performance Notes
- Broadcasting in `set_donor_arrival.(...)` allocates temporary arrays.
- For large-scale simulations, consider buffer-reusing or in-place variants.

# See also
- [`generate_recipient_arrivals`](@ref)
- [`generate_pseudo_history`](@ref)
- [`simulate_initial_recipient_list`](@ref)
- [`allocate`](@ref)
"""
function generate_donor_arrivals(donors::Vector{Donor}; origin::Date=Date(2000,1,1), nyears::Int=10, arrival_rate::Real=242.0)

    sim_end = origin + Year(nyears)

    n_arrivals = rand(Poisson(arrival_rate * nyears))
    arrival_dates = KidneyAllocation.sample_days(origin, sim_end, n_arrivals)
    sampled_donors = rand(donors, n_arrivals)
    new_donors = set_donor_arrival.(sampled_donors, arrival_dates)
    
    return new_donors

end

generate_donor_arrivals(donors)


using JLD2





@save "pseudo_history_0001.jld2" r

@time @load "pseudo_history_0001.jld2"