@testset "simulation.jl" begin
    rng = MersenneTwister(42)
    origin = Date(2000,1,1)

    # Registry (tiny)
    recipients = [
        Recipient(Date(1999,1,1), Date(1980,1,1), Date(1998,1,1), O),
        Recipient(Date(2001,1,1), Date(1975,1,1), Date(2000,6,1), A),
        Recipient(Date(2002,1,1), Date(1990,1,1), Date(2001,5,1), B),
    ]
    donors = [
        Donor(Date(2000,2,1), 40.0, 1.2, O),
        Donor(Date(2003,1,1), 55.0, 1.8, A),
    ]

    fm = nothing   # mock; allocation ignores it
    u  = 0.5

    @testset "generate_arrivals" begin
        idx, dates = generate_arrivals(eachindex(recipients), 2.0; origin=origin, nyears=1, rng=rng)
        @test all(i -> i in eachindex(recipients), idx)
        @test all(d -> origin <= d <= origin + Year(1), dates)
    end

    @testset "reconstruct_recipients / donors" begin
        idx = [1,3]
        dates = [Date(2000,1,10), Date(2000,2,20)]
        rec = reconstruct_recipients(recipients, idx, dates)
        @test length(rec) == 2
        @test rec[1] == recipients[idx[1]]
        @test rec[2] == recipients[idx[2]]

        didx = [2]
        ddates = [Date(2000,3,3)]
        don = reconstruct_donors(donors, didx, ddates)
        @test don[1] == donors[didx[1]]

        @test_throws ArgumentError reconstruct_recipients(recipients, [1,2], [Date(2000,1,1)])
        @test_throws ArgumentError reconstruct_donors(donors, [1,2], [Date(2000,1,1)])
    end

    @testset "simulate_initial_state_indexed" begin
        rng2 = MersenneTwister(123)

        final_idx, shifted_dates = simulate_initial_state_indexed(
            recipients, donors, fm, u;
            start_date = origin,
            nyears = 1,
            donor_rate = 1.0,
            recipient_rate = 1.0,
            origin_date = origin,
            rng = rng2,
        )

        @test length(final_idx) == length(shifted_dates)
        @test all(i -> i in eachindex(recipients), final_idx)
        # shifted dates are Dates; no strict bounds needed, but they should be Dates
        @test all(d -> d isa Date, shifted_dates)
    end
end