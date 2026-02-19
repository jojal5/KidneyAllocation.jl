using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2, Random, Test

using KidneyAllocation

recipients = [
        Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981,1,1), Date(1997,1,1), Date(2000,6,1), A, 69, 2403, 7, 35, 4, 103, 10),
        Recipient(Date(1963,1,1), Date(1998,1,1), Date(2001,5,1), B, 25, 68, 67, 5102, 11, 16, 20),
    ]


"""
    sim_cpra_compatibility(recipient) -> Bool

Return `true` if the recipient is compatible with a donor, based on CPRA.

Compatibility is simulated as a Bernoulli trial with probability
`1 - recipient.cpra/100`.
"""
function sim_cpra_compatibility(recipient::Recipient)
    return rand() ≥ recipient.cpra / 100
end

@testset "sim_cpra_compatibility" begin

    import KidneyAllocation.sim_cpra_compatibility

    recipients = [
        Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981,1,1), Date(1997,1,1), Date(2000,6,1), A, 69, 2403, 7, 35, 4, 103, 100),
    ]

    @test sim_cpra_compatibility(recipients[1]) == true
    @test sim_cpra_compatibility(recipients[2]) == false

end



@time sim_cpra_compatibility.(recipients)
