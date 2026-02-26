
@testset "preprocess.jl" begin

    @testset "parse_hla_int" begin
        import KidneyAllocation.parse_hla_int

        @test parse_hla_int("24") == 24
        @test parse_hla_int("24L") == 24
        @test parse_hla_int("24Low") == 24

        @test ismissing(parse_hla_int(missing))
        @test_throws ArgumentError parse_hla_int("invalid")

    end

    @testset "infer_recipient_expiration_date" begin
        import KidneyAllocation.infer_recipient_expiration_date

        df = CSV.read("data/expiration_date.csv", DataFrame)

        g = groupby(df, :CAN_ID)

        @test infer_recipient_expiration_date(g[1]) === nothing
        @test infer_recipient_expiration_date(g[2]) == Date(2012, 2, 22)
        @test infer_recipient_expiration_date(g[3]) === nothing
        @test infer_recipient_expiration_date(g[4]) == Date(2008, 3, 6)

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

    @testset "recipient_arrival_departure" begin

        import KidneyAllocation.recipient_arrival_departure

        @testset "recipient permanently removed" begin
            df = DataFrame(CAN_ID=2311, CAN_LISTING_DT=Date(2009, 9, 17),
                OUTCOME=["X", "0", "1", "1"], UPDATE_TM=[Date(2013, 1, 10), Date(2012, 8, 3), Date(2012, 2, 29), Date(2011, 2, 1)])

            arrival, departure = recipient_arrival_departure(df)

            @test arrival == Date(2009, 9, 17)
            @test departure == Date(2013, 1, 10)
        end


        @testset "transplanted recipient" begin
            df = DataFrame(CAN_ID=5695, CAN_LISTING_DT=Date(2017, 6, 19),
                OUTCOME=["TX", "1"], UPDATE_TM=[Date(2017, 9, 14), Date(2017, 7, 14)])

            arrival, departure = recipient_arrival_departure(df)

            @test arrival == Date(2017, 6, 19)
            @test departure == Date(2017, 9, 14)
        end

        @testset "still wainting recipient" begin
            df = DataFrame(CAN_ID=18725, CAN_LISTING_DT=Date(2021, 6, 14),
                OUTCOME=["1"], UPDATE_TM=[Date(2021, 9, 22)])

            arrival, departure = recipient_arrival_departure(df)

            @test arrival == Date(2021, 6, 14)
            @test departure == Date(2100, 1, 1)
        end

        @testset "several recipients" begin
            df = DataFrame(CAN_ID=[5695, 1000], CAN_LISTING_DT=Date(2017, 6, 19),
                OUTCOME=["TX", "1"], UPDATE_TM=[Date(2017, 9, 14), Date(2017, 7, 14)])

            @test_throws AssertionError recipient_arrival_departure(df)
        end
    end


end