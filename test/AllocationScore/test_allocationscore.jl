# Note: All donor–recipient combinations in these tests are synthetic and do not correspond to real individuals, to maintain confidentiality.

@testset "allocationscore.jl" begin

    @testset "score()" begin
        
        d = Donor(DateTime(2023, 07, 30, 12, 0, 0),
            33, A, 2, 23, 44, 51, 11, 12, 1.)

        r = Recipient(Date(1968, 9, 14),
            Date(2016, 12, 22),
            DateTime(2023, 5, 27, 9, 35, 0),
            A, 3, 24, 18, 35, 17, 11, 46)

        @test KidneyAllocation.score(d, r) ≈ 12.83 rtol=.01 

    end

end