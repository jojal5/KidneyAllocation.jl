
@testset "preprocess.jl" begin
    
    @testset "parse_hla_int" begin
        import KidneyAllocation.parse_hla_int

        @test parse_hla_int("24") == 24
        @test parse_hla_int("24L") == 24
        @test parse_hla_int("24Low") == 24
        
    end


end