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


