@testset "decision.jl" begin
    @testset "offered_recipients()" begin

        import KidneyAllocation.offered_recipients

        # Missing columns
        df = DataFrame(STATUS="TX", DON_CAN_SCORE=30)
        @test_throws AssertionError offered_recipients(df)

        df = DataFrame(DON_ID=1, DON_CAN_SCORE=30)
        @test_throws AssertionError offered_recipients(df)

        df = DataFrame(DON_ID=1, STATUS="TX")
        @test_throws AssertionError offered_recipients(df)

        # Not a single donor
        df = DataFrame(DON_ID=[1, 2], STATUS="TX", DON_CAN_SCORE=[30, 29])
        @test_throws AssertionError offered_recipients(df)

        # Single row transplanted
        df = DataFrame(DON_ID=1, STATUS="TX", DON_CAN_SCORE=30)
        df_offered = offered_recipients(df)
        @test df_offered.DON_CAN_SCORE == [30]

        # The first offer is refused, but the second is accepted (not sorted)
        df = DataFrame(DON_ID=[1, 1], STATUS=["TX", missing], DON_CAN_SCORE=[29, 30])
        df_offered = offered_recipients(df)
        @test df_offered.DON_CAN_SCORE == [30, 29]

        # All the offers are refused
        df = DataFrame(DON_ID=[1, 1], STATUS=[missing, missing], DON_CAN_SCORE=[30, 29])
        df_offered = offered_recipients(df)
        @test df_offered.DON_CAN_SCORE == [30, 29]

        # The fourth offers is the last accepted.
        df = DataFrame(DON_ID=1, STATUS=[missing, "TX", missing, "TX", missing], DON_CAN_SCORE=[30, 29, 28, 27, 26])
        df_offered = offered_recipients(df)
        @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27]

        # Three offers are wrongly marked as aTX.
        df = DataFrame(DON_ID=1, STATUS=[missing, "TX", missing, "TX", "TX"], DON_CAN_SCORE=[30, 29, 28, 27, 26])
        df_offered = offered_recipients(df)
        @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27]
    end
end