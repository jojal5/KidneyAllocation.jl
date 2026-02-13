using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2, Random

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

dm = GLMDecisionModel(fm)

threshold = fit_decision_threshold(dm)
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


donor = new_donors[1000]
arrival = donor.arrival
# eligible_index = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
# eligible_recipients = waiting_recipients[eligible_index]
# chosen_recipient_idx = allocate_one_donor(donor, eligible_recipients, fm, u)
# chosen_recipient = eligible_recipients[chosen_recipient_idx]
# @time allocate_one_donor(eligible_recipients, donor, fm, u)

eligible_mask = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_index = findall(eligible_mask)

ranked_indices = KidneyAllocation.rank_eligible_indices_by_score(donor,waiting_recipients, eligible_index )

@time chosen_index = allocate_one_donor(donor, waiting_recipients, dm, threshold)

chosen_recipient = waiting_recipients[chosen_index]

# Sanity checks
score.(donor, chosen_recipient)
KidneyAllocation.acceptance_probability(dm, chosen_recipient, donor)
KidneyAllocation.decide(dm, threshold,chosen_recipient, donor)


@time ind = allocate(new_donors, waiting_recipients, dm, threshold)


# Sanity checks
score(new_donors[100], waiting_recipients[ind[100]])
KidneyAllocation.acceptance_probability(dm, waiting_recipients[ind[100]], new_donors[100])
KidneyAllocation.decide(dm, threshold, waiting_recipients[ind[100]], new_donors[100])

@time ind = allocate(new_donors, waiting_recipients, dm, threshold, until = 1)

findlast(ind .!= 0)


using BenchmarkTools
@btime KidneyAllocation.acceptance_probability(dm, waiting_recipients[ind[100]], new_donors[100])
@btime KidneyAllocation.rank_eligible_indices_by_score(donor,waiting_recipients, eligible_index )






import KidneyAllocation: generate_arrivals, reconstruct_donors, reconstruct_recipients, simulate_initial_state_indexed

@time ind, d = generate_arrivals(eachindex(recipients), 100)
@time r = recipients[ind]
@time shift_recipient_timeline.(r, d)
@time reconstruct_recipients(recipients, ind, d)

@time ind, d = generate_arrivals(eachindex(donors), 100)
@time r = donors[ind]
@time set_donor_arrival.(r, d)
@time reconstruct_donors(donors, ind, d)

@time ind, d = simulate_initial_state_indexed(recipients, donors, fm, u)
reconstruct_recipients(recipients, ind, d)

# 15863 - 0.093 years
r = Recipient(Date(1945,03,11),Date(1980,11,2),Date(2000,1,1), O,
    29, 29, 44, 44, 7, 7,
    0)
# r = shift_recipient_timeline(r, Date(2000,1,1))


# 15472 - 0.063 years
r = Recipient(Date(1931,09,17),Date(1999,10,14),Date(2000,1,1), O,
    1, 2, 35, 61, 103, 13,
    0)
# r = shift_recipient_timeline(r, Date(2000,1,1))


@load "src/SyntheticData/initial_waiting_lists_indexed.jld2"

waiting_time = Vector{Float64}(undef, 1000)

# for iSim in 1:1
for iSim in 1:length(waiting_indices)
    arrival_dates = origin_date .+ Day.(waiting_day_offsets[iSim])

    waiting_recipients = reconstruct_recipients(recipients, waiting_indices[iSim], arrival_dates)

    ind, d = generate_arrivals(eachindex(recipients), λᵣ)
    append!(waiting_recipients, reconstruct_recipients(recipients, ind, d))

    ind, d = generate_arrivals(eachindex(donors), λₒ)
    new_donors = reconstruct_donors(donors, ind, d)

    # Add the recipient of interest
    pushfirst!(waiting_recipients, r)

    ind = allocate(waiting_recipients, new_donors, fm, u; until=1)
    ind[findlast(ind .!= 0)]

    waiting_time[iSim] = fractionalyears_between(Date(2000,1,1), new_donors[findlast(ind .!= 0)].arrival)

end

# using Gadfly

plot(y=waiting_time, Geom.boxplot)