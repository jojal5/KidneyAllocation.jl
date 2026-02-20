using Pkg
Pkg.activate(".")

using Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: get_eligible_recipient_indices, rank_eligible_indices_by_score, allocate_one_donor

"""
    allocate_until_next_offer(donors, recipients, dm, ind) -> Int

Return the index of the first donor for which recipient `ind` would receive an
offer under the ranked allocation process, or `0` if none occurs.
"""
function allocate_until_next_offer(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    ind::Int,
)
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


import KidneyAllocation: allocate_one_donor, allocate

@testset "allocate_until_next_offer" begin

    @import KidneyAllocation.allocate_until_next_offer
    @load "data/tree_decision_model.jld2"

    recipients = [
        Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
        Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), A, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1979, 1, 1), Date(1996, 1, 1), Date(1997, 1, 1), AB, 68, 203, 39, 77, 15, 17, 0),
    ]

    donors = [
        Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.6),
        Donor(Date(2001, 1, 2), 55, O, 2, 33, 37, 53, 4, 11, 1.2),
        Donor(Date(2001, 1, 3), 55, O, 2, 33, 37, 53, 4, 11, 1.2),
        Donor(Date(2001, 1, 4), 55, AB, 2, 33, 37, 53, 4, 11, 1.2),
    ]

    # Recipient 1 refuses the firts offer
    @test allocate_until_next_offer(donors, recipients, dm, 1) == 1

    # Recipient 2 refuses the first offer
    @test allocate_until_next_offer(donors, recipients, dm, 2) == 1

    # Recipient 3 has not been offered a donor
    @test allocate_until_next_offer(donors, recipients, dm, 3) == 0

    # Recipient 4 has been offered donor 4
    @test allocate_until_next_offer(donors, recipients, dm, 4) == 4

end
