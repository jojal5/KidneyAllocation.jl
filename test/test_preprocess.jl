
@testset "preprocess.jl" begin

    @testset "fill_hla_pairs" begin
        import KidneyAllocation.fill_hla_pairs!

        df = DataFrame(DON_A1=[2, 2, 2, missing], DON_A2=[3, missing, 3, 3], DON_B1=[missing, 4, 4, 4], DON_B2=[5, 5, missing, missing], DON_DR1=[missing, 6, 6, 6], DON_DR2=[7, 7, 7, missing])

        fill_hla_pairs!(df, "DON")

        @test df.DON_A1 == [2, 2, 2, 3]
        @test df.DON_A2 == [3, 2, 3, 3]
        @test df.DON_B1 == [5, 4, 4, 4]
        @test df.DON_B2 == [5, 5, 4, 4]
        @test df.DON_DR1 == [7, 6, 6, 6]
        @test df.DON_DR2 == [7, 7, 7, 6]
    end

    @testset "harmonize_col!()" begin

        import KidneyAllocation.harmonize_col!

        df = DataFrame(CAN_ID=1, CAN_DIAL_DT=Date(2000, 1, 1), UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)])
        append!(df, DataFrame(CAN_ID=2, CAN_DIAL_DT=Date(2000, 1, 2), UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)]))

        harmonize_col!(df, col=:CAN_DIAL_DT)
        G = groupby(df, :CAN_ID)
        @test all(G[1].CAN_DIAL_DT .== Date(2000, 1, 1))
        @test all(G[2].CAN_DIAL_DT .== Date(2000, 1, 2))

        df = DataFrame(CAN_ID=1, CAN_DIAL_DT=[Date(2000, 1, 1), missing, missing], UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)])
        append!(df, DataFrame(CAN_ID=2, CAN_DIAL_DT=[missing, Date(2000, 1, 2), Date(2000, 1, 2)], UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)]))

        harmonize_col!(df, col=:CAN_DIAL_DT)
        G = groupby(df, :CAN_ID)
        @test all(G[1].CAN_DIAL_DT .== Date(2000, 1, 1))
        @test all(G[2].CAN_DIAL_DT .== Date(2000, 1, 2))

        df = DataFrame(CAN_ID=1, CAN_DIAL_DT=[Date(2000, 1, 1), Date(2100, 1, 1), missing], UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)])
        append!(df, DataFrame(CAN_ID=2, CAN_DIAL_DT=Date(2000, 1, 2), UPDATE_TM=[Date(2001, 1, 1), Date(2002, 1, 1), Date(2003, 1, 1)]))

        @test_throws ArgumentError harmonize_col!(df, col=:CAN_DIAL_DT)

    end

    @testset "parse_hla_int" begin
        import KidneyAllocation.parse_hla_int

        @test parse_hla_int("24") == 24
        @test parse_hla_int("24L") == 24
        @test parse_hla_int("24Low") == 24

        @test ismissing(parse_hla_int(missing))
        @test_throws ArgumentError parse_hla_int("invalid")

    end

    @testset "coalesce_listing!" begin
        import KidneyAllocation.coalesce_listing!

        df = DataFrame(CAN_ID=[1, 2, 3], CAN_DIAL_DT=[Date(2000, 1, 1), Date(2001, 1, 1), Date(2002, 1, 1)], CAN_LISTING_DT=[Date(2001, 1, 1), Date(1990, 1, 1), missing])
        coalesce_listing!(df)

        @test df.CAN_LISTING_DT[1] == Date(2001, 1, 1)
        @test df.CAN_LISTING_DT[2] == Date(1990, 1, 1)
        @test df.CAN_LISTING_DT[3] == Date(2002, 1, 1)

        allowmissing!(df, 2)
        df.CAN_DIAL_DT[1] = missing

        @test_throws ArgumentError coalesce_listing!(df)

    end

    df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)
    G = groupby(df, :CAN_ID)

    @testset "check_df_at_most_one_exit_outcome()" begin
        import KidneyAllocation.check_df_at_most_one_exit_outcome

        @test check_df_at_most_one_exit_outcome(G[1]) === nothing
        @test_throws ArgumentError check_df_at_most_one_exit_outcome(G[2])
        @test check_df_at_most_one_exit_outcome(G[3]) === nothing
        @test_throws ArgumentError check_df_at_most_one_exit_outcome(G[4])
        @test check_df_at_most_one_exit_outcome(G[5]) === nothing
        @test check_df_at_most_one_exit_outcome(G[6]) === nothing

    end

    @testset "filter_dialysis_episode()" begin
        import KidneyAllocation.filter_dialysis_episode

        filtered_df = filter_dialysis_episode(G[1])
        @test nrow(filtered_df) == 7

        filtered_df = filter_dialysis_episode(G[2])
        @test nrow(filtered_df) == 2

        filtered_df = filter_dialysis_episode(G[3])
        @test nrow(filtered_df) == 12

        filtered_df = filter_dialysis_episode(G[4])
        @test nrow(filtered_df) == 2

        filtered_df = filter_dialysis_episode(G[5])
        @test nrow(filtered_df) == 2

        filtered_df = filter_dialysis_episode(G[6])
        @test nrow(filtered_df) == 1

    end

    @testset "first_exit_after_dialysis()" begin
        import KidneyAllocation.first_exit_after_dialysis

        @test first_exit_after_dialysis(G[1]) == DateTime(2012,4,24,14,33,41)
        @test first_exit_after_dialysis(G[2]) == DateTime(2004,1,14,22,55,12)
        @test first_exit_after_dialysis(G[3]) == DateTime(2012, 07, 01, 17, 0, 10)
        @test first_exit_after_dialysis(G[4]) == DateTime(2004, 04, 20, 00, 30, 14)
        @test first_exit_after_dialysis(G[5]) === nothing
        @test first_exit_after_dialysis(G[6]) === nothing

    end

    @testset "get_exit_date()" begin
        import KidneyAllocation.get_exit_date
        @test get_exit_date(df[1, :]) == DateTime(2012,4,24,14,33,41) # TX
        @test get_exit_date(df[2, :]) === nothing # 1
        @test get_exit_date(df[4, :]) === nothing #0
        @test get_exit_date(df[8, :]) == DateTime(2012,12,3,15,49,00) # TX VIVANT
        @test get_exit_date(df[29, :]) == DateTime(2013,12,4,11,26,25) # X
    end

    @testset "get_expiration_date()" begin
        import KidneyAllocation.get_expiration_date

        @test get_expiration_date(first(G[1])) === nothing
        @test get_expiration_date(first(G[2])) === DateTime(2012, 12, 3, 15, 49, 00)
        @test get_expiration_date(first(G[3])) === nothing
        @test get_expiration_date(first(G[4])) == DateTime(2013, 12, 4, 11, 26, 25)
        @test get_expiration_date(first(G[5])) === nothing
        @test get_expiration_date(first(G[6])) === nothing

    end

    @testset "recipient_active_waiting_proportion" begin
        import KidneyAllocation.recipient_active_waiting_proportion

        @test recipient_active_waiting_proportion(G[1]) ≈ 0.8855 atol=1e-4
        @test_throws ArgumentError recipient_active_waiting_proportion(G[2])
        @test recipient_active_waiting_proportion(G[3]) ≈ 0.7929 atol=1e-4
        @test_throws ArgumentError recipient_active_waiting_proportion(G[4])
        @test recipient_active_waiting_proportion(G[5]) ≈ 1. atol=1e-4
        @test recipient_active_waiting_proportion(G[6]) ≈ 1. atol=1e-4
    end

    @testset "recipient_from_row()" begin
        import KidneyAllocation.recipient_from_row

        df = DataFrame(CAN_ID=1, CAN_BTH_DT=Date(1970, 1, 1),
        CAN_DIAL_DT=Date(1999, 1, 1),
        CAN_LISTING_DT=Date(2000, 1, 1),
        CAN_BLOOD="O",
        CAN_A1=3, CAN_A2=3, CAN_B1=7, CAN_B2=8, CAN_DR1=7, CAN_DR2=8,
        CPRA = 0,
        OUTCOME = "TX",
        active_waiting_proportion = 1.0)

        r = first(df)

        recipient = recipient_from_row(r)

        @test recipient.birth == Date(1970, 1, 1)
        @test recipient.dialysis == Date(1999, 1, 1)
        @test recipient.arrival == Date(2000, 1, 1)
        @test recipient.blood == O
        @test get_HLA(recipient) == (3, 3, 7, 8, 7, 8)
        @test recipient.cpra == 0
        @test isnothing(recipient.expiration_date)

    end

    @testset "donor_from_row()" begin

        import KidneyAllocation: donor_from_row, evaluate_kdri, parse_abo, creatinine_mgdl, get_HLA

        df = DataFrame(DON_ID=1, DON_DEATH_TM=Date(2000, 1, 1), DON_AGE=60, DON_BLOOD="O", HEIGHT=1.8, WEIGHT=60., HYPERTENSION=1, DIABETES=0, DEATH=4, CREATININE=8., DCD=0,
            DON_A1=3, DON_A2=3, DON_B1=7, DON_B2=8, DON_DR1=7, DON_DR2=8)

        r = first(df)

        d = donor_from_row(r)

        @test d.arrival == Date(2000, 1, 1)
        @test d.age == 60
        @test d.blood == O
        @test get_HLA(d) == (3, 3, 7, 8, 7, 8)
        @test d.kdri ≈ 3.7152 atol = 1e-4

    end

    df = CSV.read("data/offers.csv", DataFrame)
    G = groupby(df, :DON_ID)

    @testset "filter_offers_until_last_acceptance()" begin
        import KidneyAllocation.filter_offers_until_last_acceptance

        filtered_df = filter_offers_until_last_acceptance(G[1])
        @test nrow(filtered_df) == 13
        @test uppercase(strip(last(filtered_df.DECISION))) == "ACCEPTATION"

        filtered_df = filter_offers_until_last_acceptance(G[2])
        @test nrow(filtered_df) == 42

        filtered_df = filter_offers_until_last_acceptance(G[3])
        @test nrow(filtered_df) == 3
        @test uppercase(strip(last(filtered_df.DECISION))) == "ACCEPTATION"

        filtered_df = filter_offers_until_last_acceptance(G[4])
        @test nrow(filtered_df) == 1
        @test uppercase(strip(last(filtered_df.DECISION))) == "ACCEPTATION"

    end

    @testset "infer_number_of_kidneys_offered()" begin
        import KidneyAllocation.infer_number_of_kidneys_offered
        @test infer_number_of_kidneys_offered(G[1]) == 2
        @test infer_number_of_kidneys_offered(G[2]) == 2 # No kidney has been accepted but we suppos that 2 were offered
        @test infer_number_of_kidneys_offered(G[3]) == 2 # 3 recipients have accepted, but only 2 kidneys were transplanted
        @test infer_number_of_kidneys_offered(G[4]) == 1 # Only one kidney offered and it is accepted
    end



    



    

    


    

    

    # @testset "kidneys_given_by_donor()" begin

    #     import KidneyAllocation.kidneys_given_by_donor

    #     # Missing columns
    #     df = DataFrame(STATUS="TX")
    #     @test_throws AssertionError kidneys_given_by_donor(df)
    #     df = DataFrame(DON_ID=1)
    #     @test_throws AssertionError kidneys_given_by_donor(df)

    #     df = DataFrame(DON_ID=1, STATUS=[missing, missing, "TX", "TX"])
    #     append!(df, DataFrame(DON_ID=2, STATUS=[missing, "TX", missing]))
    #     append!(df, DataFrame(DON_ID=3, STATUS=[missing, missing, missing]))
    #     d = kidneys_given_by_donor(df)

    #     @test d[1] == 2
    #     @test d[2] == 1
    #     @test d[3] == 0

    # end

end

# @testset "offered_recipients()" begin

#     import KidneyAllocation.filter_offers_until_last_acceptance

#     # Single row transplanted
#     df = DataFrame(DON_ID=1, STATUS="TX", DON_CAN_SCORE=30)
#     df_offered = filter_offers_until_last_acceptance(df)
#     @test df_offered.DON_CAN_SCORE == [30]

#     # The first offer is refused, but the second is accepted (not sorted)
#     df = DataFrame(DON_ID=[1, 1], STATUS=["TX", missing], DON_CAN_SCORE=[29, 30])
#     df_offered = filter_offers_until_last_acceptance(df)
#     @test df_offered.DON_CAN_SCORE == [30, 29]

#     # All the offers are refused
#     df = DataFrame(DON_ID=[1, 1], STATUS=[missing, missing], DON_CAN_SCORE=[30, 29])
#     df_offered = filter_offers_until_last_acceptance(df)
#     @test isempty(df_offered)

#     # The fourth offers is the last accepted.
#     df = DataFrame(DON_ID=1, STATUS=[missing, "TX", missing, "TX", missing], DON_CAN_SCORE=[30, 29, 28, 27, 26])
#     df_offered = filter_offers_until_last_acceptance(df)
#     @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27]

#     # Three offers are marked as TX (we keep all of them)
#     df = DataFrame(DON_ID=1, STATUS=[missing, "TX", missing, "TX", "TX"], DON_CAN_SCORE=[30, 29, 28, 27, 26])
#     df_offered = filter_offers_until_last_acceptance(df)
#     @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27, 26]
# end

