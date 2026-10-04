using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, JLD2, Random

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"


## Estimation of the recipient arrival rate

df_recipient = KidneyAllocation.load_recipient(recipient_filepath)

G = groupby(df_recipient, :CAN_ID)

n = 0
for g in G
    filtered_df = KidneyAllocation.filter_outcomes(g)
    if filtered_df.CAN_LISTING_DT[1] ≥ DateTime(2011,12,31,23,59,59) && filtered_df.CAN_LISTING_DT[1] < DateTime(2020,1,1,0,0,0)
        n+=1
    end
end

recipient_arrival_rate = n/8

## Build recipient registry by CAN_ID

recipient_registry = KidneyAllocation.build_recipient_registry(recipient_filepath, cpra_filepath)

## Estimate the donor arrival rate

df_donors = KidneyAllocation.load_donor(donor_filepath)

G = groupby(df_donors, :DON_ID)

n = 0
for g in G
    r = first(g)
    if  r.DON_DEATH_TM ≥ Date(2011,12,31) && r.DON_DEATH_TM < Date(2020,1,1)
        n+=1
    end
end

donor_arrival_rate = n/6

## Build donor registry by DON_ID

donor_registry = KidneyAllocation.build_donor_registry(donor_filepath)

# ## Retrieve recipients for which waiting time has to be estimated

# # CPRA < 80
# # Registered between 2012 and 2020
# # Only for their first registration if past transplant

# recipient_ids = Int64[]

# for id in keys(recipient_registry)
#     r = recipient_registry[id]
#     if r.arrival ≥ Date(2012,1,1) && r.arrival < Date(2020,1,1)
#         if r.cpra ≤ 80
#             push!(recipient_ids, id)
#         end
#     end  
# end

## Dictionary of exit statuses and date of exit

exit_status = Dict{Int64, Tuple{String, DateTime}}()

for k in keys(recipient_registry)

    g = filter(row->row.CAN_ID == k, df_recipient)
    filtered_df = KidneyAllocation.filter_outcomes(g)
    sort!(filtered_df, :UPDATE_TM, rev= true)

    exit_status[k] = (filtered_df.OUTCOME[1], filtered_df.UPDATE_TM[1])

end


## Retrieve transplanted recipients for which waiting time has to be estimated

# CPRA < 80
# Registered between 2012 and 2020
# Transplanted before 2020
# Only for their first registration if past transplant

recipient_ids = Int64[]

for id in keys(recipient_registry)
    r = recipient_registry[id]
    if r.arrival ≥ Date(2012,1,1) && r.arrival < Date(2020,1,1)
        if r.cpra ≤ 80
            if first(exit_status[id]) == "TX"  && last(exit_status[id])< Date(2020,1,1)
                push!(recipient_ids, id)
            end
        end
    end  
end

## Selection of a recipient

i = 4
id = recipient_ids[i]

recipient = recipient_registry[id]

## Retrieve the outcome and the date of exit (if any)

exit_status[id]

obs_waiting_time = Date(last(exit_status[id])) - Date(recipient.arrival)

filter(row -> row.CAN_ID == id, df_recipient)

## Retrieve the waiting list when recipient arrived

# TODO: il faut retirer les patients qui ont été transplantés
initial_recipients = Recipient[]

for k in keys(recipient_registry)
    if KidneyAllocation.is_registered(recipient_registry[k], recipient.arrival)
        push!(initial_recipients, recipient_registry[k])
    end
end

ind = findfirst(initial_recipients .== recipient)
# Sanity check
initial_recipients[ind] == recipient

## Generate recipient arrivals for the next nyears

nyears = 5

# Number of new recipients
nᵣ = rand(Poisson(recipient_arrival_rate * nyears)) 
# Arrival dates                             
tᵣ = KidneyAllocation.sample_days(recipient.arrival, recipient.arrival + Year(nyears), nᵣ)
# Sampled CAN_ID
sampled_can_id = rand(keys(recipient_registry), nᵣ)
# Sampled recipients with the adjusted timeline 
new_recipients = Vector{Recipient}(undef, nᵣ)
for (i,id) in enumerate(sampled_can_id)
    sampled_recipient = recipient_registry[id]
    new_recipients[i] = KidneyAllocation.shift_recipient_timeline(sampled_recipient, tᵣ[i])
end

# Sanity check
# KidneyAllocation.get_arrival.(new_recipients) == tᵣ

waiting_recipients = vcat(initial_recipients, new_recipients)

# Verify the position of the considered recipient (Sanity check)
waiting_recipients[ind] == recipient


## Generate donor arrivals for the next nyears

# Number of recipients
nₒ = rand(Poisson(donor_arrival_rate * nyears)) 
# Arrival dates                             
tₒ = KidneyAllocation.sample_days(recipient.arrival, recipient.arrival + Year(nyears), nₒ)
# Sampled DON_ID
sampled_don_id = rand(keys(donor_registry), nₒ)

## Sampled donors 

kidney_by_don_id = KidneyAllocation.kidneys_given_by_donor(df_donors)

# Sampled donors with the adjusted arrival and the number of given kidneys
donors = Donor[]
for (i,id) in enumerate(sampled_don_id)
    sampled_donor = donor_registry[id]
    for j = 1:kidney_by_don_id[id]
        push!(donors, KidneyAllocation.set_donor_arrival(sampled_donor, tₒ[i]))
    end
end

## Load decision model

@load "src/SyntheticData/GLMDecisionModel.jld2"


## Allocate until first offer

# @time offer_ind = KidneyAllocation.allocate_until_next_offer(donors, waiting_recipients, dm, ind)

# donors[offer_ind].arrival - recipient.arrival


## Refactor

attributed_recipient_index = zeros(Int64, length(donors))
is_unallocated::AbstractVector{<:Bool}=trues(length(waiting_recipients))
# donor = donors[1]

for (i,donor) in enumerate(donors)

    eligible_indices = KidneyAllocation.get_eligible_recipient_indices(donor, waiting_recipients, is_unallocated)

    scored_indices = KidneyAllocation.rank_eligible_indices_by_score(donor, waiting_recipients, eligible_indices)

    # Sanity check
    KidneyAllocation.score.(donor, waiting_recipients[scored_indices])

    decision = KidneyAllocation.decide.(Ref(dm), waiting_recipients[scored_indices], donor)

    if any(decision)
        attributed_recipient_index[i] = scored_indices[findfirst(decision)]
        is_unallocated[attributed_recipient_index[i]] = false
    else
        attributed_recipient_index[i] = 0
    end

end

attributed_donor_index = findfirst(attributed_recipient_index .== ind)

estimated_waiting_time = donors[attributed_donor_index].arrival - recipient.arrival











## Test the allocation for one donor

import KidneyAllocation: is_active, is_abo_compatible, allocate_one_donor

donor = donors[100]
arrival = donor.arrival

eligible_mask = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_index = findall(eligible_mask)

ranked_indices = KidneyAllocation.rank_eligible_indices_by_score(donor, waiting_recipients, eligible_index )

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



