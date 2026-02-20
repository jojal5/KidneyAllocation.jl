
"""
    allocate_one_donor(donor, recipients, dm, is_unallocated) -> Int

Return the index of the first recipient who accepts the donor offer among
eligible and ranked candidates, or `0` if none accept.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    is_unallocated::BitVector=trues(length(recipients))
)
    eligible_indices = get_eligible_recipient_indices(donor, recipients, is_unallocated)
    ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

    return allocate_one_donor(donor, recipients, dm, ranked_indices)

end

"""
    allocate_one_donor(donor, recipients, dm, ranked_indices) -> Int

Return the index of the first accepting recipient in `ranked_indices`,
or `0` if none accept.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ranked_indices::AbstractVector{<:Int}
)
    accepted = decide(dm, recipients[ranked_indices], donor)

    if any(accepted)
        ind = findfirst(accepted)
        return ranked_indices[ind]
    else
        return 0
    end

end

"""
    allocate_one_donor(donor, recipients, dm, ranked_indices::Int) -> Int

Return `ranked_indices` if the corresponding recipient accepts the offer,
or `0` otherwise.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ranked_indices::Int
)
    
    return allocate_one_donor(donor, recipients, dm, [ranked_indices])
end

"""
    allocate(donors, recipients, dm) -> Vector{Int}

Allocate each donor in `donors` to at most one recipient in `recipients` using `dm`.

Returns a vector of allocated recipient indices (0 if unallocated).
"""
function allocate(donors::Vector{Donor}, recipients::Vector{Recipient}, dm::AbstractDecisionModel)

    is_unallocated = trues(length(recipients))                 
    allocated_recipient_index = zeros(Int64,length(donors))

    for donor_idx in eachindex(donors)

        donor = donors[donor_idx]

        allocated_recipient_index[donor_idx] = allocate_one_donor(donor, recipients, dm, is_unallocated)

        if allocated_recipient_index[donor_idx] != 0
            is_unallocated[allocated_recipient_index[donor_idx]] = false
        end
    end

    return allocated_recipient_index
end

"""
    allocate_until_transplant(donors, recipients, dm, ind) -> Int

Allocate donors sequentially using `dm` until recipient `ind` is allocated.
Return the donor index, or `0` if no transplant occurs.
"""
function allocate_until_transplant(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ind::Int,
)

    @assert 1 ≤ ind ≤ length(recipients) "Recipient index should be in 1 ≤ ind ≤ $(length(recipients)), got ind = $ind."
    
    is_unallocated = trues(length(recipients))

    for donor_idx in eachindex(donors)
        donor = donors[donor_idx]

        allocated_recipient_index = allocate_one_donor(donor, recipients, dm, is_unallocated)

        if allocated_recipient_index == ind
            return donor_idx
        end

        if allocated_recipient_index != 0
            is_unallocated[allocated_recipient_index] = false
        end
    end

    return 0
end

"""
    allocate_until_next_offer(donors, recipients, dm, ind) -> Int

Return the index of the first donor for which recipient `ind` would receive an offer, or `0` if none occurs.
"""
function allocate_until_next_offer(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ind::Int,
)

    @assert 1 ≤ ind ≤ length(recipients) "Recipient index should be in 1 ≤ ind ≤ $(length(recipients)), got ind = $ind."

    is_unallocated = trues(length(recipients))

    for donor_idx in eachindex(donors)
        donor = donors[donor_idx]

        eligible_indices = get_eligible_recipient_indices(donor, recipients, is_unallocated)
        ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

        allocated_recipient_index = allocate_one_donor(donor, recipients, dm, ranked_indices)

        # check whether ind would be offered before acceptance (or no acceptance)
        pos_ind = findfirst(==(ind), ranked_indices)
        if pos_ind !== nothing
            if allocated_recipient_index == 0
                return donor_idx
            else
                pos_alloc = findfirst(==(allocated_recipient_index), ranked_indices)
                if pos_alloc !== nothing && pos_ind ≤ pos_alloc
                    return donor_idx
                end
            end
        end

        # update availability after allocation
        if allocated_recipient_index != 0
            is_unallocated[allocated_recipient_index] = false
        end
    end

    return 0
end

"""
    get_eligible_recipient_indices(donor, recipients, is_unallocated) 

Return the indices of recipients eligible to receive an offer from `donor`.

A recipient is considered eligible if it:
- is currently unallocated,
- is active at the donor arrival date,
- is ABO-compatible with the donor, and
- is CPRA is lower a random value.

# Arguments
- `donor::Donor`: Donor being allocated.
- `recipients::Vector{Recipient}`: Current waiting list.
- `is_unallocated::BitVector`: Optional mask indicating which recipients are
  still available for allocation (default: all `true`).

# Returns
- `Vector{Int}`: Indices into `recipients` identifying eligible recipients.
"""
function get_eligible_recipient_indices(
    donor::Donor,
    recipients::Vector{Recipient},
    is_unallocated::AbstractVector{<:Bool} = trues(length(recipients)),
)

    arrival = donor.arrival

    eligible_mask = copy(is_unallocated)
    eligible_mask .&= is_active.(recipients, arrival)
    eligible_mask .&= is_abo_compatible.(donor, recipients)
    eligible_mask .&= sim_cpra_compatibility.(recipients)

    return findall(eligible_mask)
end


"""
    rank_eligible_indices_by_score(donor, recipients, eligible_indices)

Return a new vector of recipient indices ranked by decreasing attribution score
for `donor`.

The returned vector contains the same elements as `eligible_indices`, reordered
so that `score(donor, recipients[i])` is decreasing.

# Arguments
- `donor::Donor`: Donor being allocated.
- `recipients::Vector{Recipient}`: Current waiting list.
- `eligible_indices::AbstractVector{<:Integer}`: Indices into `recipients`
  identifying eligible recipients.
"""
function rank_eligible_indices_by_score(
    donor::Donor,
    recipients::Vector{Recipient},
    eligible_indices::AbstractVector{<:Integer},
)
    scores = score.(Ref(donor), recipients[eligible_indices])
    p = sortperm(scores; rev=true)
    return eligible_indices[p]
end

