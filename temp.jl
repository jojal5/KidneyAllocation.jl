using Pkg
Pkg.activate(".")

using Dates, JLD2, Random, Test

using KidneyAllocation


function allocate_until_next_offer(donors::Vector{Donor}, recipients::Vector{Recipient}, dm::AbstractDecisionModel, ind::Int)
    
    is_unallocated = trues(length(recipients))                 
    
    for donor_idx in eachindex(donors)

        donor = donors[donor_idx]

        eligible_indices = get_eligible_recipient_indices(donor, recipients, is_unallocated)
        ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

        allocated_recipient_index[donor_idx] = allocate_one_donor(donor, recipients, dm, ranked_indices)

        if ind ∈ ranked_indices end


        if allocated_recipient_index[donor_idx] != 0
            is_unallocated[allocated_recipient_index[donor_idx]] = false
        end

        if allocated_recipient_index[donor_idx] == until
            break
        end

    end

    return allocated_recipient_index
end




"""
    comes_before(v::Vector{Int}, a::Int, b::Int) -> Bool

Return `true` if `a` appears before `b` in `v`, and `false` otherwise.
"""
function comes_before(v::AbstractVector{<:Int}, a::Int, b::Int)
    @assert a ≠ b "`a` should be different than `b, got a = b = $a"

    found_a = false
    for x in v
        if x == a
            found_a = true
        elseif x == b
            return found_a
        end
    end
    return false
end

using Test

@testset "comes_before" begin
    import KidneyAllocation.comes_before

    v = [1, 2, 3, 4, 5]

    @test_throws AssertionError comes_before(v, 4, 4)
    @test comes_before(v, 2, 4) == true
    @test comes_before(v, 4, 2) == false
    @test comes_before(v, -2, 4) == false
    @test comes_before(v, 2, 0) == false

end
