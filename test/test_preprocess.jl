
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


end