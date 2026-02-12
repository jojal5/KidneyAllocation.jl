@testset "simulation.jl" begin
    rng = MersenneTwister(42)
    origin = Date(2000,1,1)

    # Registry (tiny and fakes)
    recipients = [
        Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981,1,1), Date(1997,1,1), Date(2000,6,1), A, 69, 2403, 7, 35, 4, 103, 10),
        Recipient(Date(1963,1,1), Date(1998,1,1), Date(2001,5,1), B, 25, 68, 67, 5102, 11, 16, 20),
    ]
    donors = [
        Donor(Date(2000,1,1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5),
        Donor(Date(2001,1,1), 55, A, 2, 33, 37, 53, 4, 11, 1.6),
    ]

    fm = nothing   # mock; allocation ignores it
    u  = 0.5

    @testset "generate_arrivals" begin
        import KidneyAllocation.generate_arrivals
        idx, dates = generate_arrivals(eachindex(recipients), 2.0; origin=origin, nyears=1, rng=rng)
        @test all(i -> i in eachindex(recipients), idx)
        @test all(d -> origin <= d <= origin + Year(1), dates)
    end

    @testset "reconstruct_recipients" begin
        import KidneyAllocation: get_HLA, reconstruct_recipients
        idx = [1,3]
        dates = [Date(2000,1,10), Date(2000,2,20)]
        rec = reconstruct_recipients(recipients, idx, dates)
        @test length(rec) == length(idx)

        for i in eachindex(idx)
            shift = dates[i]-recipients[idx[i]].arrival
            @test rec[i].arrival == dates[i]
            @test rec[i].birth == recipients[idx[i]].birth + shift
            @test rec[i].dialysis == recipients[idx[i]].dialysis + shift

            @test rec[i].blood == recipients[idx[i]].blood
            @test get_HLA(rec[i]) == get_HLA(recipients[idx[i]])
            @test rec[i].cpra == recipients[idx[i]].cpra
        end
    end

    @testset "reconstruct_donors" begin
        import KidneyAllocation: get_HLA, reconstruct_donors
        idx = [1]
        dates = [Date(2000,3,3)]

        don = reconstruct_donors(donors, idx, dates)
        @test length(don) == length(idx)

        for i in eachindex(idx)
            @test don[i].arrival == dates[i]
            @test don[i].age == donors[idx[i]].age
            @test don[i].blood == donors[idx[i]].blood
            @test get_HLA(don[i]) == get_HLA(donors[idx[i]])
            @test don[i].kdri == donors[idx[i]].kdri
        end
    end

# TODO: test when DecisionModel is created
    # @testset "simulate_initial_state_indexed" begin
    #     rng2 = MersenneTwister(123)

    #     final_idx, shifted_dates = simulate_initial_state_indexed(
    #         recipients, donors, fm, u;
    #         start_date = origin,
    #         nyears = 1,
    #         donor_rate = 1.0,
    #         recipient_rate = 1.0,
    #         origin_date = origin,
    #         rng = rng2,
    #     )

    #     @test length(final_idx) == length(shifted_dates)
    #     @test all(i -> i in eachindex(recipients), final_idx)
    #     # shifted dates are Dates; no strict bounds needed, but they should be Dates
    #     @test all(d -> d isa Date, shifted_dates)
    # end
end