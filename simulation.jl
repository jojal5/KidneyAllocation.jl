using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

## Estimation of the recipient arrival rate

import KidneyAllocation.load_recipient

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

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
tᵣ = KidneyAllocation.sample_days(origin, origin + Year(nyears), nᵣ)
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
tₒ = KidneyAllocation.sample_days(origin, origin + Year(nyears), nₒ)
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

donor = donors[100]
arrival = donor.arrival

eligible_mask = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_index = findall(eligible_mask)

ranked_indices = KidneyAllocation.rank_eligible_indices_by_score(donor,waiting_recipients, eligible_index )

@time p= KidneyAllocation.acceptance_probability(dm, waiting_recipients[ranked_indices], donor)

@time chosen_index = allocate_one_donor(donor, waiting_recipients[ranked_indices], dm)

chosen_recipient = waiting_recipients[ranked_indices[chosen_index]]

# Sanity checks
KidneyAllocation.score.(donor, chosen_recipient)
KidneyAllocation.acceptance_probability(dm, chosen_recipient, donor)
KidneyAllocation.decide(dm, chosen_recipient, donor)

## Test the allocation for all donors

import KidneyAllocation.allocate

@time ind = allocate(donors, waiting_recipients, dm)

# TODO: Bcp trop d'offres non acceptées. Changer le modèle de décision ou bien forcer les candidats à les accepter. 
count(ind .== 0)


# Sanity checks - does not work if ind[idx] == 0, i.e. if the kidney is not attributed
idx = 1000
KidneyAllocation.score(donors[idx], waiting_recipients[ind[idx]])
KidneyAllocation.acceptance_probability(dm, waiting_recipients[ind[idx]], donors[idx])
KidneyAllocation.decide(dm, waiting_recipients[ind[idx]], donors[idx])

@time ind = KidneyAllocation.allocate_until_next_offer(donors, waiting_recipients, dm, 100)


@time ind = KidneyAllocation.allocate_until_transplant(donors, waiting_recipients, dm, 100)



# TODO - verify the time before transplant and first offer using real recipients

pushfirst!(waiting_recipients, waiting_recipients[1]) # To be replaced by the target recipient

# 15863 - 0.093 years ≈ 34 days
r = Recipient(Date(1945,03,11),Date(1980,11,2),Date(2000,1,1), O,
    29, 29, 44, 44, 7, 7,
    0)
r = KidneyAllocation.shift_recipient_timeline(r, Date(2014,1,1))

waiting_recipients[1] = r

ind = KidneyAllocation.allocate_until_transplant(donors, waiting_recipients, dm, 1)
donors[ind].arrival - waiting_recipients[1].arrival 


# 15472 - 0.063 years ≈ 23 days
r = Recipient(Date(1931,09,17),Date(1999,10,14),Date(2000,1,1), O,
    1, 2, 35, 61, 103, 13,
    0)
r = shift_recipient_timeline(r, Date(2014,1,1))

waiting_recipients[1] = r
ind = KidneyAllocation.allocate_until_transplant(donors, waiting_recipients, dm, 1)
donors[ind].arrival - waiting_recipients[1].arrival


# 6072 - 891 days before TX 
r = recipient_by_CAN_ID[6072]
r = KidneyAllocation.shift_recipient_timeline(r, Date(2014,1,1))

waiting_recipients[1] = r
ind = KidneyAllocation.allocate_until_transplant(donors, waiting_recipients, dm, 1)
donors[ind].arrival - waiting_recipients[1].arrival



