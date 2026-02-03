using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, load_recipient, is_active, is_expired, is_abo_compatible
import KidneyAllocation: load_donor, build_donor_registry
import KidneyAllocation: shift_recipient_timeline, set_donor_arrival
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision



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



ind_active = is_active.(recipients, Date(2013, 12, 31))

waiting_recipients = recipients[ind_active]

nᵣ = rand(Poisson(λᵣ * 10))
tᵣ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nᵣ)

sampled_recipients = rand(recipients, nᵣ)
new_recipients = KidneyAllocation.shift_recipient_timeline.(sampled_recipients, tᵣ)

append!(waiting_recipients, new_recipients)

ind_active = is_active.(waiting_recipients, Date(2013, 12, 31))

waiting_recipients[ind_active]


nₒ = rand(Poisson(λₒ * 10))
tₒ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nₒ)

sampled_donors = rand(donors, nₒ)
new_donors = set_donor_arrival.(sampled_donors, tₒ)

KidneyAllocation.get_arrival.(new_donors)



donor = new_donors[1]

@time for donor in new_donors

    filter!(x -> !is_expired(x, donor.arrival), waiting_recipients)

    id_recipient = eachindex(waiting_recipients)

    ind_active = is_active.(waiting_recipients, donor.arrival)
    ind_compatible = is_abo_compatible.(KidneyAllocation.get_bloodtype.(waiting_recipients), KidneyAllocation.get_bloodtype(donor))

    eligible_id = id_recipient[ind_active.&&ind_compatible]

    s = KidneyAllocation.score.(donor, waiting_recipients[eligible_id])

    ind = sortperm(s; rev=true)

    permute!(eligible_id, ind)

    accepting_recipient_id = 0

    for iter in eachindex(eligible_id)
        id = eligible_id[iter]
        decision = get_decision(donor, waiting_recipients[id], fm, u)
        if decision
            accepting_recipient_id = id
            break
        end
    end


    # Sanity checks
    # is_active(waiting_recipients[accepting_recipient_id ], donor.arrival)
    # is_abo_compatible.(KidneyAllocation.get_bloodtype.(waiting_recipients[accepting_recipient_id]), KidneyAllocation.get_bloodtype(donor))
    # KidneyAllocation.score.(donor, waiting_recipients[accepting_recipient_id])

    if accepting_recipient_id != 0
        deleteat!(waiting_recipients, accepting_recipient_id)
    end

end

KidneyAllocation.is_expired(waiting_recipients[1], Date(204, 1, 1))

waiting_recipients

is_active.(waiting_recipients, Date(2024, 1, 1))

waiting_recipients[is_active.(waiting_recipients, Date(2024, 1, 1))]

r = filter(x -> !is_expired(x, Date(2024, 1, 1)), waiting_recipients)

arr = KidneyAllocation.get_arrival.(r)

minimum(arr)





eligible_id = Int[]
empty!(eligible_id)
sizehint!(eligible_id, length(waiting_recipients))






import KidneyAllocation.score
import KidneyAllocation: years_between, fractionalyears_between

@inline blood_str(b::ABOGroup) =
    b == O ? "O" : b == A ? "A" : b == B ? "B" : "AB"

function get_decision(donor::Donor, recipient::Recipient,
    fm::StatsModels.TableRegressionModel, u::Real)::Bool
    is_abo_compatible(donor.blood, recipient.blood) || return false

    arrival = donor.arrival

    row = (DON_AGE=donor.age,
        KDRI=donor.kdri,
        CAN_AGE=years_between(recipient.birth, arrival),
        CAN_WAIT=fractionalyears_between(recipient.dialysis, arrival),
        CAN_BLOOD=blood_str(recipient.blood))

    return predict(fm, (row,))[1] > u
end


function allocate_one_donor!(
    waiting_recipients::Vector{Recipient},
    donor::Donor,
    fm,
    u;
    eligible_id::Vector{Int}=Int[],
    scores::Vector{Float64}=Float64[],
    perm::Vector{Int}=Int[])

    arrival = donor.arrival
    donor_bt = donor.blood

    empty!(eligible_id)
    empty!(scores)

    sizehint!(eligible_id, length(waiting_recipients))
    sizehint!(scores, length(waiting_recipients))

    @inbounds for i in eachindex(waiting_recipients)
        r = waiting_recipients[i]
        is_active(r, arrival) || continue
        is_abo_compatible(r.blood, donor_bt) || continue

        push!(eligible_id, i)
        push!(scores, score(donor, r))
    end

    isempty(eligible_id) && return 0

    # sortperm! with reusable buffer
    resize!(perm, length(scores))
    sortperm!(perm, scores; rev=true)

    accepting_id = 0
    @inbounds for k in perm
        id = eligible_id[k]
        if get_decision(donor, waiting_recipients[id], fm, u)
            accepting_id = id
            break
        end
    end

    return accepting_id
end

donor = new_donors[1]

@time allocate_one_donor(waiting_recipients, donor, fm, u)





# function allocate_one_donor(
#     waiting_recipients::Vector{Recipient},
#     donor::Donor,
#     fm,
#     u)::Int64

#     arrival = donor.arrival

#     eligible = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)

#     eligible_index = findall(eligible)

#     attribution_score = Vector{Float64}(undef, length(eligible_index))

#     for i in eachindex(eligible_index)
#         attribution_score[i] = score(donor, waiting_recipients[eligible_index[i]])
#     end

#     ind = sortperm(attribution_score, rev=true)

#     accepting_index = 0

#     for i in eachindex(ind)
#         if get_decision(donor, waiting_recipients[ind[i]], fm, u)
#             accepting_index = i
#             break
#         end
#     end

#     return accepting_index
# end

function allocate_one_donor(
    eligible_recipients::Vector{Recipient},
    donor::Donor,
    fm,
    u)::Int64

    arrival = donor.arrival

    attribution_score = score.(donor, eligible_recipients)

    ind = sortperm(attribution_score, rev=true)

    accepting_index = 0

    for i in eachindex(ind)
        if get_decision(donor, eligible_recipients[ind[i]], fm, u)
            accepting_index = ind[i]
            break
        end
    end

    return accepting_index
end


donor = new_donors[1]
arrival = donor.arrival
eligible_index = is_active.(waiting_recipients, arrival) .&& is_abo_compatible.(donor, waiting_recipients)
eligible_recipients = waiting_recipients[eligible_index]
allocate_one_donor(eligible_recipients, donor, fm, u)

@time allocate_one_donor(eligible_recipients, donor, fm, u)







function allocate(recipients::Vector{Recipient}, donors::Vector{Donor}, fm, u)

    non_allocated_recipients = trues(length(recipients))
    ind_allocated = Vector{Int64}(undef, length(donors))

    for i in eachindex(donors)

        arrival  = donors[i].arrival
        eligible_id = is_active.(recipients, arrival) .&& is_abo_compatible.(donors[i], recipients) .&& non_allocated_recipients
        eligible_recipients = recipients[eligible_id]

        ind = findall(eligible_id)

        id = allocate_one_donor(eligible_recipients , donors[i], fm, u)

        if id != 0
            ind_allocated[i] = ind[id]
            non_allocated_recipients[ind[id]] = false
        else
            ind_allocated[i] = 0
        end
    end

    return ind_allocated
end

function allocate_GPT(recipients::Vector{Recipient}, donors::Vector{Donor}, fm, u)

    is_unallocated = trues(length(recipients))                 # BitVector
    allocated_recipient_index = Vector{Int}(undef, length(donors))

    eligible_mask = falses(length(recipients))                 # reused BitVector
    eligible_indices = Int[]                                   # reused Vector{Int}

    for donor_idx in eachindex(donors)

        donor = donors[donor_idx]
        arrival = donor.arrival

        # Build eligibility mask with minimal temporaries
        eligible_mask .= is_unallocated                         # start from availability
        eligible_mask .&= is_active.(recipients, arrival)        # active at arrival
        eligible_mask .&= is_abo_compatible.(donor, recipients)  # ABO compatible

        # Reuse eligible_indices instead of allocating a new vector each time
        empty!(eligible_indices)
        append!(eligible_indices, findall(eligible_mask))

        if isempty(eligible_indices)
            allocated_recipient_index[donor_idx] = 0
            continue
        end

        # Find accepting recipient among eligible ones (no slicing recipients[eligible_mask])
        chosen_recipient = 0
        for k in eachindex(eligible_indices)
            r_idx = eligible_indices[k]
            if get_decision(donor, recipients[r_idx], fm, u)
                chosen_recipient = r_idx
                break
            end
        end

        if chosen_recipient != 0
            allocated_recipient_index[donor_idx] = chosen_recipient
            is_unallocated[chosen_recipient] = false
        else
            allocated_recipient_index[donor_idx] = 0
        end
    end

    return allocated_recipient_index
end


@time allocate(waiting_recipients, new_donors, fm, u)

@time allocate_GPT(waiting_recipients, new_donors, fm, u)



donors = new_donors
recipients = waiting_recipients

non_allocated_recipients = trues(length(recipients))
ind_allocated = Vector{Int64}(undef, length(donors))

i = 1
arrival  = donors[i].arrival
eligible_id = is_active.(recipients, arrival) .&& is_abo_compatible.(donors[i], recipients) .&& non_allocated_recipients
eligible_recipients = recipients[eligible_id]

@time ind = findall(eligible_id)

id = allocate_one_donor(eligible_recipients , donors[i], fm, u)

ind[id]



@time allocate_one_donor!(waiting_recipients, new_donors[1], fm, u)

function swapdeleteat!(v::Vector, i::Int)
    v[i] = v[end]
    pop!(v)
    return v
end





function allocate!(donors::Vector{Donor}, recipients::Vector{Recipient}, fm, u)

    eligible_id = Int[]
    scores = Float64[]
    perm = Int[]

    # Counting the non-attributed kidneys.
    n = 0

    for donor in donors

        # If you must drop expired each donor, still do it, but avoid extra work later
        filter!(r -> !is_expired(r, donor.arrival), recipients)

        id = allocate_one_donor!(recipients, donor, fm, u;
            eligible_id=eligible_id, scores=scores, perm=perm)
        if id != 0
            # deleteat!(waiting_recipients, id)  # see next section for faster option
            swapdeleteat!(recipients, id)
        else
            n += 1
        end

        # return id
    end
end



@time allocate!(new_donors, waiting_recipients, fm, u)




function allocate(donors::Vector{Donor}, recipients::Vector{Recipient}, fm, u)

    eligible_id = Int[]
    scores = Float64[]
    perm = Int[]



    for donor in donors

        # If you must drop expired each donor, still do it, but avoid extra work later
        filter!(r -> !is_expired(r, donor.arrival), recipients)

        id = allocate_one_donor!(recipients, donor, fm, u;
            eligible_id=eligible_id, scores=scores, perm=perm)
        if id != 0
            # deleteat!(waiting_recipients, id)  # see next section for faster option
            swapdeleteat!(recipients, id)
        else
            n += 1
        end

        # return id
    end
end