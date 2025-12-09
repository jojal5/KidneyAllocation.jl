@testset "Quality.jl" begin
    
    @testset "calculate_epts" begin
        import KidneyAllocation.calculate_epts

        r = Recipient(Date(1967, 9, 14),
            Date(2017, 3, 14),
            DateTime(2022, 7, 21, 12, 30, 0),
            A,
            11, 24, 27, 52, 17, 11,
            0
        )
        current_date = Date(2025,12,8)

        @test calculate_epts(r, false, false, current_date) ≈ 99.56 atol=.1
        @test calculate_epts(r, true, false, current_date) ≈ 99.73 atol=.1
        @test calculate_epts(r, true, true, current_date) ≈ 99.85 atol=.1

    end

    @testset "calculate_kdri" begin
        # Source: A Guide to Calculating and Interpreting the Kidney Donor Profile Index, April 21, 2025.
        @test KidneyAllocation.evaluate_kdri(52, 183, 81, true, false, true, 1.7, true) ≈ 1.71235565748184
    end


end