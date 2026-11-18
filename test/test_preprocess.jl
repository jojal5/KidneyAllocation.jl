
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


    @testset "filter_first_waitlist_episode" begin
        import KidneyAllocation.filter_first_waitlist_episode

        @test nrow(filter_first_waitlist_episode(G[1])) == nrow(G[1])
        @test nrow(filter_first_waitlist_episode(G[2])) == 2
        @test nrow(filter_first_waitlist_episode(G[3])) == nrow(G[3])
        @test nrow(filter_first_waitlist_episode(G[4])) == 2
        @test nrow(filter_first_waitlist_episode(G[5])) == nrow(G[5])
        @test nrow(filter_first_waitlist_episode(G[6])) == nrow(G[6])
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

end

