using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

## Estimation recipient arrival rate

import KidneyAllocation.load_recipient

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
# recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

df_recipient = load_recipient(recipient_filepath)
select!(df_recipient, Not(:NB_PAST_TRANS))
dropmissing!(df_recipient)

G = groupby(df_recipient, :CAN_ID)

n = 0
for g in G
    r = first(g)
    if  Date(2012,12,31) < r.CAN_LISTING_DT < Date(2020,1,1)
        n+=1
    end
end

recipient_arrival_rate = n/6

## Build recipient registry by CAN_ID

recipient_by_CAN_ID = KidneyAllocation.build_recipient_registry(recipient_filepath, cpra_filepath)

## Retrieve active waiting recipients at the origin

import KidneyAllocation: recipient_arrival_departure

origin = Date(2014,1,1)

G = groupby(df_recipient, :CAN_ID)

origin = Date(2014,1,1)
active_can_id = Int[]

for g in G
    arrival, departure = recipient_arrival_departure(g)
    if arrival ≤ origin < departure
        push!(active_can_id, g.CAN_ID[1])
    end
end

initial = [recipient_registry_by_can_id[i] for i in active_can_id]


## Generate recipient arrivals for the next nyears

import KidneyAllocation.shift_recipient_timeline

nyears = 10

# Number of recipients
nᵣ = rand(Poisson(recipient_arrival_rate * 10)) 
# Arrival dates                             
tᵣ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nᵣ)
# Sampled CAN_ID
sampled_can_id = rand(keys(recipient_registry_by_can_id), nᵣ)
# Sampled recipients with the adjusted timeline 
new_recipients = Vector{Recipient}(undef, nᵣ)
for (i,id) in enumerate(sampled_can_id)
    sampled_recipient = recipient_registry_by_can_id[id]
    new_recipients[i] = shift_recipient_timeline(sampled_recipient, tᵣ[i])
end

# Sanity check
# KidneyAllocation.get_arrival.(new_recipients) == tᵣ

waiting_recipients = vcat(initial, new_recipients)

       
## Estimate the donor arrival rate

import KidneyAllocation: load_donor

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_donors = load_donor(donor_filepath)

G = groupby(df_donors, :DON_ID)

n = 0
for g in G
    r = first(g)
    if  Date(2012,12,31) < r.DON_DEATH_TM < Date(2020,1,1)
        n+=1
    end
end

donor_arrival_rate = n/6

## Build donor registry by DON_ID

import KidneyAllocation: donor_from_row

donor_registry_by_don_id = Dict{Int, Donor}()

for g in G
    r = first(g)
    donor_registry_by_don_id[r.DON_ID] = donor_from_row(r)
end








import KidneyAllocation: kidneys_given_by_donor


## Load decision model

@load "src/SyntheticData/TreeDecisionModel.jld2"









                                        # Add the new recipients to the original list


# Generate new donors for the next 10 years
nₒ = rand(Poisson(λₒ * 10))                                                             # Number of recipients
tₒ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nₒ)             # Arrival dates
sampled_donors = rand(donors, nₒ)                                                       # Sampled recipients
new_donors = set_donor_arrival.(sampled_donors, tₒ)                                     # Adjust the arrival dates

# KidneyAllocation.get_arrival.(new_donors)


donor = new_donors[4]
arrival = donor.arrival
# eligible_index = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
# eligible_recipients = waiting_recipients[eligible_index]
# chosen_recipient_idx = allocate_one_donor(donor, eligible_recipients, fm, u)
# chosen_recipient = eligible_recipients[chosen_recipient_idx]
# @time allocate_one_donor(eligible_recipients, donor, fm, u)

eligible_mask = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_index = findall(eligible_mask)

ranked_indices = KidneyAllocation.rank_eligible_indices_by_score(donor,waiting_recipients, eligible_index )

@time KidneyAllocation.acceptance_probability(dm, waiting_recipients[ranked_indices], donor)

@time chosen_index = allocate_one_donor(donor, waiting_recipients[ranked_indices], dm)

chosen_recipient = waiting_recipients[ranked_indices[chosen_index]]

# Sanity checks
score.(donor, chosen_recipient)
KidneyAllocation.acceptance_probability(dm, chosen_recipient, donor)
KidneyAllocation.decide(dm, chosen_recipient, donor)


@time ind = allocate(new_donors, waiting_recipients, dm)


# Sanity checks
score(new_donors[100], waiting_recipients[ind[100]])
KidneyAllocation.acceptance_probability(dm, waiting_recipients[ind[100]], new_donors[100])
KidneyAllocation.decide(dm, waiting_recipients[ind[100]], new_donors[100])

@time ind = allocate(new_donors, waiting_recipients, dm, until = 1)

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