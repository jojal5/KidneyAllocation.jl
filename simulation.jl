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

initial = [recipient_by_CAN_ID[i] for i in active_can_id]


## Generate recipient arrivals for the next nyears

import KidneyAllocation.shift_recipient_timeline

nyears = 10

# Number of recipients
nᵣ = rand(Poisson(recipient_arrival_rate * nyears)) 
# Arrival dates                             
tᵣ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nᵣ)
# Sampled CAN_ID
sampled_can_id = rand(keys(recipient_by_CAN_ID), nᵣ)
# Sampled recipients with the adjusted timeline 
new_recipients = Vector{Recipient}(undef, nᵣ)
for (i,id) in enumerate(sampled_can_id)
    sampled_recipient = recipient_by_CAN_ID[id]
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

donor_by_don_id = KidneyAllocation.build_donor_registry(donor_filepath)

## Generate donor arrivals for the next nyears

# Number of recipients
nₒ = rand(Poisson(donor_arrival_rate * nyears)) 
# Arrival dates                             
tₒ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nₒ)
# Sampled DON_ID
sampled_don_id = rand(keys(donor_by_don_id), nₒ)

## Sampled donors 

kidney_by_don_id = KidneyAllocation.kidneys_given_by_donor(df_donors)

# Sampled donors with the adjusted arrival and the number of given kidneys
donors = Donor[]
for (i,id) in enumerate(sampled_don_id)
    sampled_donor = donor_by_don_id[id]
    for j = 1:kidney_by_don_id[id]
        push!(donors, KidneyAllocation.set_donor_arrival(sampled_donor, tₒ[i]))
    end
end


## Load decision model

@load "src/SyntheticData/TreeDecisionModel.jld2"

## Test the allocation for one donor

import KidneyAllocation: is_active, is_abo_compatible, allocate_one_donor

donor = donors[1]
arrival = donor.arrival

eligible_mask = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_index = findall(eligible_mask)

ranked_indices = KidneyAllocation.rank_eligible_indices_by_score(donor,waiting_recipients, eligible_index )

@time KidneyAllocation.acceptance_probability(dm, waiting_recipients[ranked_indices], donor)

@time chosen_index = allocate_one_donor(donor, waiting_recipients[ranked_indices], dm)

chosen_recipient = waiting_recipients[ranked_indices[chosen_index]]

# Sanity checks
KidneyAllocation.score.(donor, chosen_recipient)
KidneyAllocation.acceptance_probability(dm, chosen_recipient, donor)
KidneyAllocation.decide(dm, chosen_recipient, donor)

## Test the allocation for all donors

import KidneyAllocation.allocate

@time ind = allocate(donors, waiting_recipients, dm)


is_unallocated = trues(length(waiting_recipients))                 
allocated_recipient_index = zeros(Int64,length(donors))

    for donor_idx in eachindex(donors)
        println(donor_idx)
        donor = donors[donor_idx]

        allocated_recipient_index[donor_idx] = allocate_one_donor(donor, waiting_recipients, dm, is_unallocated)

        if allocated_recipient_index[donor_idx] != 0
            is_unallocated[allocated_recipient_index[donor_idx]] = false
        end
    end


donor_idx = 14
donor = donors[donor_idx]

eligible_indices = KidneyAllocation.get_eligible_recipient_indices(donor, waiting_recipients, is_unallocated)

ind = KidneyAllocation.is_abo_compatible.(donor, waiting_recipients)
count(ind)


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