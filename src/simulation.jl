"""
    allocate_one_donor(
        eligible_recipients::Vector{Recipient},
        donor::Donor,
        fm,
        u
    ) -> Int64

Determine which eligible recipient (if any) accepts an offer from a given donor.

This function ranks all `eligible_recipients` by decreasing allocation
`score(donor, recipient)`, then evaluates recipients in that order using
`get_decision`. The first recipient who accepts the offer is selected.

If no recipient accepts, the function returns `0`.

# Arguments
- `eligible_recipients::Vector{Recipient}`: Vector of recipients who are
  eligible for allocation (active, compatible, and unallocated).
- `donor::Donor`: The donor whose organ is being allocated.
- `fm`: Fitted statistical model used by `get_decision`.
- `u`: Acceptance threshold or decision parameter.

# Returns
- `Int64`: Index of the accepting recipient in `eligible_recipients`,
  or `0` if no recipient accepts.

# Algorithm
1. Compute allocation scores for all eligible recipients.
2. Sort recipients in decreasing order of score.
3. Evaluate `get_decision` sequentially in that order.
4. Return the index of the first accepting recipient.

# Performance Notes
- This implementation allocates temporary arrays for scores and sorting.
- Time complexity is `O(n log n)` in the number of eligible recipients.
- For large-scale simulations, consider using a buffer-reusing or
  partial-sorting variant to reduce allocations and runtime.

# See also
- [`score`](@ref)
- [`get_decision`](@ref)
"""
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

"""
    allocate(
        recipients::Vector{Recipient},
        donors::Vector{Donor},
        fm,
        u;
        until::Int64 = -9999
    ) -> Vector{Int64}

Simulate a sequential kidney allocation process over a stream of donors, with
optional early stopping.

For each donor (processed in the order provided by `donors`), the function:

1. Identifies eligible recipients among `recipients` who are
   (i) still unallocated, (ii) active at the donor's `arrival`, and
   (iii) ABO-compatible with the donor.
2. Extracts the subvector of eligible recipients and calls
   [`allocate_one_donor`](@ref) to select an accepting recipient using a
   score-ranked offer sequence.
3. Maps the selected index back to the original `recipients` vector and marks
   the recipient as allocated.

If no eligible recipient exists, or if none accepts, the donor is left
unallocated and the output for that donor is `0`.

If the keyword argument `until` is provided and a recipient with index `until`
is allocated, the allocation process terminates immediately.

# Arguments
- `recipients::Vector{Recipient}`: Candidate recipients. This vector is not
  modified; allocation status is tracked internally.
- `donors::Vector{Donor}`: Donor stream processed sequentially.
- `fm`: Fitted statistical model used by `get_decision` inside
  `allocate_one_donor`.
- `u`: Decision parameter (e.g., acceptance threshold) passed through to
  `allocate_one_donor` and `get_decision`.

# Keyword Arguments
- `until::Int64 = -9999`: Optional stopping index. If a recipient with this
  index (in `recipients`) is allocated, the simulation stops immediately.
  By default (`-9999`), no early stopping is applied.

# Returns
- `Vector{Int64}`: For each donor `donors[i]`, returns the index (in the original
  `recipients` vector) of the allocated recipient, or `0` if no allocation
  occurred. Entries after early stopping (if any) remain `0`.

# Behavior and Assumptions
- Each recipient can be allocated at most once.
- Eligibility is evaluated at each donor's `arrival` time using `is_active` and
  `is_abo_compatible`.
- Within each donor, recipient selection is delegated to
  `allocate_one_donor`, which ranks eligible recipients by decreasing
  `score(donor, recipient)` and applies `get_decision`.

# Index Mapping
`allocate_one_donor` returns an index relative to the eligible-recipient
subvector `recipients[eligible_mask]`. This function maps that index back to
the original `recipients` vector using `eligible_indices`.

# Performance Notes
- Reuses `eligible_mask` and `eligible_indices` to reduce per-donor allocations.
- The slicing operation `recipients[eligible_mask]` allocates a new vector of
  eligible recipients at each iteration.
- Broadcasting in `is_active.(...)` and `is_abo_compatible.(...)`, and calls to
  `findall`, also generate temporary allocations.
- For large-scale simulations, consider index-based variants to eliminate
  slicing and broadcasting.

# See also
- [`allocate_one_donor`](@ref)
- [`score`](@ref)
- [`get_decision`](@ref)
- [`is_active`](@ref)
- [`is_abo_compatible`](@ref)
"""
function allocate(recipients::Vector{Recipient}, donors::Vector{Donor}, fm, u; until::Int64=-9999)

    is_unallocated = trues(length(recipients))                 
    allocated_recipient_index = zeros(Int64,length(donors))

    eligible_mask = falses(length(recipients))
    eligible_indices = Int64[]

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

        chosen_recipient = allocate_one_donor(recipients[eligible_mask], donor, fm, u)

        if chosen_recipient != 0
            is_unallocated[eligible_indices[chosen_recipient]] = false
            allocated_recipient_index[donor_idx] = eligible_indices[chosen_recipient]
            if eligible_indices[chosen_recipient] == until
                break
            end
        else
            allocated_recipient_index[donor_idx] = 0
        end

    end

    return allocated_recipient_index
end