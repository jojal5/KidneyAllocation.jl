@testset "decision.jl" begin


    @testset "get_eligible_recipient_indice" begin

        import KidneyAllocation.get_eligible_recipient_indices

        # Registry (tiny and fakes)
        recipients = [
            Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2000, 7, 1), O, 25, 68, 67, 5102, 11, 16, 100),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), O, 25, 68, 67, 5102, 11, 16, 0),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), AB, 25, 68, 67, 5102, 11, 16, 0),
        ]
        donors = [
            Donor(Date(2000, 1, 1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5),
            Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.6),
            Donor(Date(2002, 1, 1), 55, A, 2, 33, 37, 53, 4, 11, 1.6),
        ]

        @test get_eligible_recipient_indices(donors[1], recipients) == [1]
        @test get_eligible_recipient_indices(donors[2], recipients) == [1; 2]
        @test get_eligible_recipient_indices(donors[2], recipients, [false, true, true, true, true]) == [2]
        @test get_eligible_recipient_indices(donors[3], recipients) == [5]

    end

    @testset "rank_eligible_indices_by_score" begin

        import KidneyAllocation: rank_eligible_indices_by_score

        # Registry (tiny and fakes)
        recipients = [
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2000, 7, 1), A, 25, 68, 67, 5102, 11, 16, 0),
            Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
            Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2000, 7, 1), O, 25, 68, 67, 5102, 11, 16, 100),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), O, 25, 68, 67, 5102, 11, 16, 0),
            Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), AB, 25, 68, 67, 5102, 11, 16, 0),
        ]

        donor = Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.6)

        eligible_indices = [2; 3]
        ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

        @test ranked_indices == [3, 2]

    end


end