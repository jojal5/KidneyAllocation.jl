
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    is_unallocated::AbstractVector{<:Bool}=trues(length(recipients));
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)
    ranked_indices = rank_eligible_recipient_indices(donor, recipients, is_unallocated)

    return allocate_one_donor(donor, recipients, dm, ranked_indices; mode=mode, rng=rng)
end


"""
    allocate_one_donor(donor, recipients, dm, ranked_indices;
                       mode=:random, rng=Random.default_rng()) -> Int

Return the index of the first recipient in `ranked_indices` who accepts
`donor`, or `0` if none accepts. `mode` and `rng` are forwarded to `decide`.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ranked_indices::AbstractVector{<:Int};
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)
    isempty(ranked_indices) && return 0

    accepted = decide(dm, recipients[ranked_indices], donor; mode=mode, rng=rng)

    first_accepted = findfirst(accepted)
    return isnothing(first_accepted) ? 0 : ranked_indices[first_accepted]
end


"""
    allocate_one_donor(donor, recipients, dm, ranked_index;
                       mode=:random, rng=Random.default_rng()) -> Int

Return `ranked_index` if the corresponding recipient accepts the donor offer,
or `0` otherwise.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ranked_index::Int;
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng()
)
    
    return allocate_one_donor(donor, recipients, dm, [ranked_index]; mode=mode, rng=rng)
end

"""
    allocate(donors, recipients, dm;
             mode=:random, rng=Random.default_rng()) -> Vector{Int}

Allocate each donor to at most one recipient using `dm`.

`mode` and `rng` are forwarded to `allocate_one_donor`. Return the allocated
recipient index for each donor, using `0` when no allocation is made.
"""
function allocate(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel;
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)
    is_unallocated = trues(length(recipients))
    allocated_recipient_indices = zeros(Int, length(donors))

    for donor_idx in eachindex(donors)
        donor = donors[donor_idx]

        recipient_idx = allocate_one_donor(donor, recipients, dm, is_unallocated; mode=mode, rng=rng)

        allocated_recipient_indices[donor_idx] = recipient_idx

        if recipient_idx != 0
            is_unallocated[recipient_idx] = false
        end
    end

    return allocated_recipient_indices
end

"""
    allocate_until_transplant(donors, recipients, dm, recipient_index;
                              mode=:random, rng=Random.default_rng()) -> Int

Allocate donors sequentially until `recipient_index` is allocated.

Return the corresponding donor index, or `0` if that recipient is never
allocated. `mode` and `rng` are forwarded to `allocate_one_donor`.
"""
function allocate_until_transplant(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    recipient_index::Int;
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)
    1 ≤ recipient_index ≤ length(recipients) ||
        throw(ArgumentError(
            "Recipient index must be in 1:$(length(recipients)); got $recipient_index",
        ))

    is_unallocated = trues(length(recipients))

    for donor_idx in eachindex(donors)
        recipient_idx = allocate_one_donor(donors[donor_idx], recipients, dm, is_unallocated; mode=mode, rng=rng)

        recipient_idx == recipient_index && return donor_idx

        if recipient_idx != 0
            is_unallocated[recipient_idx] = false
        end
    end

    return 0
end

"""
    allocate_until_next_offer(donors, recipients, dm, recipient_index;
                              mode=:random, rng=Random.default_rng()) -> Int

Return the index of the first donor for which `recipient_index` would receive
an offer, or `0` if that recipient is never offered a donor.
"""
function allocate_until_next_offer(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    recipient_index::Int;
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)
    1 ≤ recipient_index ≤ length(recipients) ||
        throw(ArgumentError(
            "Recipient index must be in 1:$(length(recipients)); got $recipient_index",
        ))

    is_unallocated = trues(length(recipients))

    for donor_idx in eachindex(donors)
        donor = donors[donor_idx]

        ranked_indices = rank_eligible_recipient_indices(donor, recipients, is_unallocated)

        isempty(ranked_indices) && continue


        allocated_recipient_index = allocate_one_donor(donor, recipients, dm, ranked_indices; mode=mode, rng=rng)

        target_position = findfirst(==(recipient_index), ranked_indices)
        acceptance_position = findfirst(
            ==(allocated_recipient_index),
            ranked_indices,
        )

        if !isnothing(target_position) &&
           (isnothing(acceptance_position) ||
            target_position ≤ acceptance_position)
            return donor_idx
        end

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
    is_unallocated::AbstractVector{<:Bool} = trues(length(recipients));
    rng::AbstractRNG=Random.default_rng(),
)

    arrival = donor.arrival

    eligible_mask = copy(is_unallocated)
    eligible_mask .&= is_active.(recipients, arrival)
    eligible_mask .&= is_abo_compatible.(donor, recipients)
    eligible_mask .&= sim_cpra_compatibility.(recipients, rng=rng)

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

"""
    rank_eligible_recipient_indices(donor, recipients, is_unallocated) -> Vector{Int}

Return the indices of recipients eligible for `donor`, ordered by allocation
score.
"""
function rank_eligible_recipient_indices(
    donor::Donor,
    recipients::Vector{Recipient},
    is_unallocated::AbstractVector{<:Bool}=trues(length(recipients));
    rng::AbstractRNG=Random.default_rng(),
)
    eligible_indices = get_eligible_recipient_indices(donor, recipients, is_unallocated; rng=rng)

    isempty(eligible_indices) && return Int[]

    return rank_eligible_indices_by_score(donor, recipients, eligible_indices)
end

"""
    get_recipient_offers(focal_recipient, donors, waitlist_recipients, dm,
                         is_unallocated; ...) -> Vector{Donor}

Return donors whose offers reach `focal_recipient` while competing recipients
are allocated sequentially. The focal recipient is forced to reject every
offer and therefore remains on the waiting list.
"""
function get_recipient_offers(
    focal_recipient::Recipient,
    donors::Vector{Donor},
    waitlist_recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    is_unallocated::BitVector=trues(length(waitlist_recipients));
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)::Vector{Donor}

    length(is_unallocated) == length(waitlist_recipients) ||
        throw(ArgumentError(
            "`is_unallocated` must have one entry per waiting-list recipient",
        ))

    candidate_recipients = [waitlist_recipients; focal_recipient]
    focal_index = length(candidate_recipients)
    offered_donors = Donor[]

    for donor in donors
        # The focal recipient is always still unallocated.
        candidate_is_unallocated = BitVector([is_unallocated; true])

        ranked_indices = rank_eligible_recipient_indices(
            donor,
            candidate_recipients,
            candidate_is_unallocated;
            rng=rng,
        )
        isempty(ranked_indices) && continue

        focal_position = findfirst(==(focal_index), ranked_indices)

        if focal_position === nothing
            allocated_index = allocate_one_donor(
                donor,
                candidate_recipients,
                dm,
                ranked_indices;
                mode=mode,
                rng=rng,
            )

            allocated_index != 0 &&
                (is_unallocated[allocated_index] = false)

            continue
        end

        # A higher-ranked recipient accepting prevents an offer to the focal one.
        if focal_position > 1
            higher_priority_indices = ranked_indices[1:(focal_position - 1)]

            allocated_index = allocate_one_donor(
                donor,
                candidate_recipients,
                dm,
                higher_priority_indices;
                mode=mode,
                rng=rng,
            )

            if allocated_index != 0
                is_unallocated[allocated_index] = false
                continue
            end
        end

        # The offer reaches the focal recipient, who counterfactually refuses it.
        push!(offered_donors, donor)

        # The offer can then proceed to lower-ranked recipients.
        if focal_position < length(ranked_indices)
            lower_priority_indices = ranked_indices[(focal_position + 1):end]

            allocated_index = allocate_one_donor(
                donor,
                candidate_recipients,
                dm,
                lower_priority_indices;
                mode=mode,
                rng=rng,
            )

            allocated_index != 0 &&
                (is_unallocated[allocated_index] = false)
        end
    end

    return offered_donors
end



"""
    allocate_until_transplant(focal_recipient, donors, waitlist_recipients, dm,
                              is_unallocated; ...) -> Union{Donor,Nothing}

Process donors sequentially until `focal_recipient` accepts an offer. Return
the transplanting donor, or `nothing` if no offer is accepted.

`is_unallocated` is updated in place when a waiting-list recipient accepts an
offer.
"""
function allocate_until_transplant(
    focal_recipient::Recipient,
    donors::Vector{Donor},
    waitlist_recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    is_unallocated::BitVector=trues(length(waitlist_recipients));
    mode::Symbol=:random,
    rng::AbstractRNG=Random.default_rng(),
)::Union{Donor, Nothing}

    length(is_unallocated) == length(waitlist_recipients) ||
        throw(ArgumentError("`is_unallocated` must have one entry per waiting-list recipient"))

    candidate_recipients = [waitlist_recipients; focal_recipient]
    focal_index = length(candidate_recipients)

    for donor in donors
        # The focal recipient is always still unallocated.
        candidate_is_unallocated = BitVector([is_unallocated; true])

        ranked_indices = rank_eligible_recipient_indices(donor, candidate_recipients, candidate_is_unallocated; rng=rng)

        isempty(ranked_indices) && continue

        focal_position = findfirst(==(focal_index), ranked_indices)

        if focal_position === nothing
            allocated_index = allocate_one_donor(donor, candidate_recipients, dm, ranked_indices; mode=mode, rng=rng)

            allocated_index != 0 &&
                (is_unallocated[allocated_index] = false)

            continue
        end

        # A higher-ranked recipient accepting prevents an offer to the focal one.
        if focal_position > 1
            higher_priority_indices = ranked_indices[1:(focal_position - 1)]

            allocated_index = allocate_one_donor(donor, candidate_recipients, dm, higher_priority_indices; mode=mode, rng=rng)

            if allocated_index != 0
                is_unallocated[allocated_index] = false
                continue
            end
        end

        # The offer reaches the focal recipient. If they accept it, return the donor; otherwise, continue allocation.
        if decide(dm, focal_recipient, donor; mode=mode, rng=rng)
            return donor
        end  

        # The offer can then proceed to lower-ranked recipients.
        if focal_position < length(ranked_indices)
            lower_priority_indices = ranked_indices[(focal_position + 1):end]

            allocated_index = allocate_one_donor(donor, candidate_recipients, dm, lower_priority_indices; mode=mode, rng=rng)

            allocated_index != 0 &&
                (is_unallocated[allocated_index] = false)
        end
    end

    return nothing
end