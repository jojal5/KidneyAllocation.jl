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

    @load "data/tree_decision_model.jld2"

    @testset "allocate_one_donor()" begin
        import KidneyAllocation: allocate_one_donor

        @testset "Eligible recipient with the higher score accept the offer" begin

            # Registry (tiny and fakes)
            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            ]

            donor = Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.2)

            is_unallocated = trues(length(recipients))

            # Les deux candidats sont éligibles, les deux accepteraient l'offre et le deuxième a le score le plus élevé. C'est donc au 2e que l'offre sera attribuée.
            @test allocate_one_donor(donor, recipients, dm, is_unallocated) == 2

            # Lorsque le 2e candidat a reçu une offre, il ne fait plus partie de la compétition et le premier devrait recevoir l'offre similaire.
            is_unallocated[2] = false
            @test allocate_one_donor(donor, recipients, dm, is_unallocated) == 1

        end

        @testset "Eligible recipient with the higher score refuses the offer, the second accept" begin

            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1939, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 69, 203, 39, 77, 15, 17, 0),
            ]

            donor = Donor(Date(2001, 1, 1), 65, O, 69, 2403, 7, 35, 4, 103, 1.5)

            # Les deux candidats sont éligibles, le 2e accepterait l'offre, mais le 1er a le score le plus élevé. C'est donc au 2e que l'offre sera attribuée.
            @test allocate_one_donor(donor, recipients, dm) == 2
        end

        @testset "Eligible recipients refuse the offer" begin

            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1949, 1, 1), Date(1997, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            ]

            donor = Donor(Date(2001, 1, 1), 45, O, 2, 33, 37, 53, 4, 11, 1.6)

            # Les deux candidats sont éligibles, les deux refuseraient l'offre. L'offre n'est donc pas attribuée
            @test allocate_one_donor(donor, recipients, dm) == 0
        end
    end

    @testset "allocate" begin

        import KidneyAllocation.allocate

        @testset "allocation without refusal" begin

            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            ]

            donors = [
                Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.2),
                Donor(Date(2001, 1, 2), 55, O, 2, 33, 37, 53, 4, 11, 1.2)
            ]

            # Both donors attributed to the recipient with the highest score
            @test allocate(donors, recipients, dm) == [2, 1]

        end

        @testset "allocation begining with a refusal" begin

            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            ]

            donors = [
                Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.6),
                Donor(Date(2001, 1, 2), 55, O, 2, 33, 37, 53, 4, 11, 1.2),
                Donor(Date(2001, 1, 3), 55, O, 2, 33, 37, 53, 4, 11, 1.2)
            ]

            # First donor not accepted, both following donors attributed to the recipient with the highest score
            @test allocate(donors, recipients, dm) == [0, 2, 1]

        end

    end

    @testset "allocation_until_transplant" begin

        import KidneyAllocation.allocate_until_transplant
        @load "data/tree_decision_model.jld2"

            recipients = [
                Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), O, 69, 2403, 7, 35, 4, 103, 0),
                Recipient(Date(1969, 1, 1), Date(1995, 1, 1), Date(1996, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
            ]

            donors = [
                Donor(Date(2001, 1, 1), 55, O, 2, 33, 37, 53, 4, 11, 1.6),
                Donor(Date(2001, 1, 2), 55, O, 2, 33, 37, 53, 4, 11, 1.2),
                Donor(Date(2001, 1, 3), 55, O, 2, 33, 37, 53, 4, 11, 1.2)
            ]

            # First donor not accepted, early stop when recipient 2 accept the offer
            @test allocate_until_transplant(donors, recipients, dm, 2) == 2

    end

    @testset "allocate_until_next_offer" begin

        import KidneyAllocation.allocate_until_next_offer
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


end