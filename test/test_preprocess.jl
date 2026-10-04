
@testset "preprocess.jl" begin

    @testset "parse_hla_int" begin
        import KidneyAllocation.parse_hla_int

        @test parse_hla_int("24") == 24
        @test parse_hla_int("24L") == 24
        @test parse_hla_int("24Low") == 24

        @test ismissing(parse_hla_int(missing))
        @test_throws ArgumentError parse_hla_int("invalid")

    end

    @testset "check dataframe columns" begin

        import KidneyAllocation: check_df_columns, check_df_column_constant, check_df_at_most_one_exit_outcome
        df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)
        g = groupby(df, :CAN_ID)

        @testset "check_df_columns" begin
            @test check_df_columns(df, :CAN_ID, :UPDATE_TM, :OUTCOME) === nothing

            df2 = similar(df, 0)
            @test check_df_columns(df2, :CAN_ID, :UPDATE_TM, :OUTCOME) === nothing

            @test_throws ArgumentError check_df_columns(df, :INEXISTENT)

            df3 = df[1:2, :]
            allowmissing!(df3, 5)
            df3[1, 5] = missing

            @test_throws ArgumentError check_df_columns(df3, :CAN_ID, :UPDATE_TM, :OUTCOME)
        end

        @testset "check_df_column_constant" begin

            @test_throws ArgumentError check_df_column_constant(df, :CAN_ID)

            @test check_df_column_constant(g[1], :CAN_ID) === nothing

        end

        @testset "check_df_at_most_one_exit_outcome" begin
            @test check_df_at_most_one_exit_outcome(g[1]) === nothing
            @test_throws ArgumentError check_df_at_most_one_exit_outcome(g[2])
            @test check_df_at_most_one_exit_outcome(g[3]) === nothing
            @test_throws ArgumentError check_df_at_most_one_exit_outcome(g[4])
            @test check_df_at_most_one_exit_outcome(g[5]) === nothing
            @test check_df_at_most_one_exit_outcome(g[6]) === nothing
        end
    end

    @testset "filter_outcomes" begin

        import KidneyAllocation.filter_outcomes

        df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)

        g = groupby(df, :CAN_ID)

        is_exit_outcome(outcome::AbstractString) =
            uppercase(strip(outcome)) ∈ ("X", "TX VIVANT", "DCD", "TX")

        n_exit_outcome = [1, 1, 1, 1, 0, 0]
        first_exit_date = [DateTime(2012, 4, 24, 14, 33, 41),
            DateTime(2004, 1, 14, 22, 55, 12),
            DateTime(2012, 7, 1, 17, 0, 10),
            DateTime(2004, 4, 20, 0, 30, 14),
            nothing,
            nothing]

        for i in 1:length(g)
            filtered_df = filter_outcomes(g[i])
            exit_rows = is_exit_outcome.(filtered_df.OUTCOME)
            @test count(exit_rows) == n_exit_outcome[i]
            if !isnothing(first_exit_date[i])
                @test filtered_df[exit_rows, :UPDATE_TM][1] == first_exit_date[i]
            else
                @test all(.!(exit_rows))
            end
        end

    end

    @testset "get_last_update()" begin
        import KidneyAllocation.get_last_update

        df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)
        G = groupby(df, :CAN_ID)

        (date, status) = get_last_update(G[1])
        @test date == DateTime(2012, 04, 24, 14, 33, 41)
        @test status == "TX"

        (date, status) = get_last_update(G[2])
        @test date == DateTime(2004, 1, 14, 22, 55, 12)
        @test status == "TX"

        (date, status) = get_last_update(G[3])
        @test date == DateTime(2012, 7, 1, 17, 00, 10)
        @test status == "TX"

        (date, status) = get_last_update(G[4])
        @test date == DateTime(2004, 4, 20, 0, 30, 14)
        @test status == "X"

        (date, status) = get_last_update(G[5])
        @test date == DateTime(2023, 5, 24, 13, 57, 4)
        @test status == "0"

        (date, status) = get_last_update(G[6])
        @test date == DateTime(2023, 5, 18, 10, 10, 51)
        @test status == "1"
    end

    @testset "get_expiration_date()" begin
        import KidneyAllocation.get_expiration_date

        df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)
        G = groupby(df, :CAN_ID)

        date = get_expiration_date(G[1])
        @test date === nothing

        date = get_expiration_date(G[4])
        @test date == DateTime(2004, 4, 20, 0, 30, 14)

        date = get_expiration_date(G[5])
        @test date === nothing
    end

    @testset "get_exit_date()" begin
        import KidneyAllocation.get_exit_date

        df = CSV.read("data/unfiltered_outcomes.csv", DataFrame)
        G = groupby(df, :CAN_ID)

        date = get_exit_date(G[1])
        @test date === DateTime(2012,4,24,14,33,41)

        date = get_exit_date(G[4])
        @test date == DateTime(2004, 4, 20, 0, 30, 14)

        date = get_exit_date(G[5])
        @test date === nothing
    end

    @testset "recipient_active_waiting_proportion" begin

        import KidneyAllocation.recipient_active_waiting_proportion

        df = CSV.read("data/filtered_outcomes.csv", DataFrame)
        g = groupby(df, :CAN_ID)

        @test recipient_active_waiting_proportion(g[1]) ≈ 0.8855 atol=1e-4
        @test recipient_active_waiting_proportion(g[2]) ≈ 1. atol=1e-4
        @test recipient_active_waiting_proportion(g[3]) ≈ 0.7929 atol=1e-4
        @test recipient_active_waiting_proportion(g[4]) ≈ 1. atol=1e-4
        @test recipient_active_waiting_proportion(g[5]) ≈ 1. atol=1e-4
        @test recipient_active_waiting_proportion(g[6]) ≈ 1. atol=1e-4
    end

    @testset "build_last_cpra_registry" begin
        import KidneyAllocation.build_last_cpra_registry

        cpra_filepath = "data/candidate_cpra.csv"
        d = build_last_cpra_registry(cpra_filepath)

        @test haskey(d, 3)
        @test d[3] == 37
        @test haskey(d, 14)
        @test d[14] == 12

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

    @testset "recipient_from_row()" begin

        import KidneyAllocation.recipient_from_row

        df = DataFrame(CAN_ID=1, CAN_BTH_DT=Date(1970, 1, 1), CAN_DIAL_DT=Date(1999, 1, 1), CAN_LISTING_DT=Date(2000, 1, 1), CAN_BLOOD="O",
            CAN_A1=3, CAN_A2=3, CAN_B1=7, CAN_B2=8, CAN_DR1=7, CAN_DR2=8)

        r = first(df)

        r = recipient_from_row(r)

        @test r.birth == Date(1970, 1, 1)
        @test r.dialysis == Date(1999, 1, 1)
        @test r.arrival == Date(2000, 1, 1)
        @test r.blood == O
        @test get_HLA(r) == (3, 3, 7, 8, 7, 8)
        @test r.cpra == 0
        @test isnothing(r.expiration_date)

    end


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

    # @testset "recipient_arrival_departure" begin

    #     import KidneyAllocation.recipient_arrival_departure

    #     @testset "recipient permanently removed" begin
    #         df = DataFrame(CAN_ID=2311, CAN_LISTING_DT=Date(2009, 9, 17),
    #             OUTCOME=["X", "0", "1", "1"], UPDATE_TM=[Date(2013, 1, 10), Date(2012, 8, 3), Date(2012, 2, 29), Date(2011, 2, 1)])

    #         arrival, departure = recipient_arrival_departure(df)

    #         @test arrival == Date(2009, 9, 17)
    #         @test departure == Date(2013, 1, 10)
    #     end


    #     @testset "transplanted recipient" begin
    #         df = DataFrame(CAN_ID=5695, CAN_LISTING_DT=Date(2017, 6, 19),
    #             OUTCOME=["TX", "1"], UPDATE_TM=[Date(2017, 9, 14), Date(2017, 7, 14)])

    #         arrival, departure = recipient_arrival_departure(df)

    #         @test arrival == Date(2017, 6, 19)
    #         @test departure == Date(2017, 9, 14)
    #     end

    #     @testset "still wainting recipient" begin
    #         df = DataFrame(CAN_ID=18725, CAN_LISTING_DT=Date(2021, 6, 14),
    #             OUTCOME=["1"], UPDATE_TM=[Date(2021, 9, 22)])

    #         arrival, departure = recipient_arrival_departure(df)

    #         @test arrival == Date(2021, 6, 14)
    #         @test departure == Date(2100, 1, 1)
    #     end

    # end


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

    @testset "kidneys_given_by_donor()" begin

        import KidneyAllocation.kidneys_given_by_donor

        # Missing columns
        df = DataFrame(STATUS="TX")
        @test_throws AssertionError kidneys_given_by_donor(df)
        df = DataFrame(DON_ID=1)
        @test_throws AssertionError kidneys_given_by_donor(df)

        df = DataFrame(DON_ID=1, STATUS=[missing, missing, "TX", "TX"])
        append!(df, DataFrame(DON_ID=2, STATUS=[missing, "TX", missing]))
        append!(df, DataFrame(DON_ID=3, STATUS=[missing, missing, missing]))
        d = kidneys_given_by_donor(df)

        @test d[1] == 2
        @test d[2] == 1
        @test d[3] == 0

    end

end

