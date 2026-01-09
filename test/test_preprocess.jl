
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
        @test infer_recipient_expiration_date(g[2]) == Date(2012,2,22)
        @test infer_recipient_expiration_date(g[3]) === nothing
        @test infer_recipient_expiration_date(g[4]) == Date(2008,3,6)
        
    end


end