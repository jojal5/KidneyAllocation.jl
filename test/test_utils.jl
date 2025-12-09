@testset "utils.jl" begin

    import KidneyAllocation: _dt, days_between, fractionalyears_between, years_between 

    @testset "_dt helper conversion" begin

        # Date input
        d = Date(2024, 2, 12)
        dt = _dt(d)

        @test dt isa DateTime
        @test year(dt) == 2024
        @test month(dt) == 2
        @test day(dt) == 12
        @test hour(dt) == 0
        @test minute(dt) == 0
        @test second(dt) == 0

        # DateTime input
        dt2 = DateTime(2021, 7, 8, 13, 45, 10)
        @test _dt(dt2) === dt2

    end

    @testset "days_between()" begin
        @test days_between(Date(2020,5,10), DateTime(2020,5,15)) == 5
        @test days_between(DateTime(2020,5,15), Date(2020,5,10)) == -5
        @test days_between(Date(2020,1,1), Date(2020,1,1)) == 0

        @test days_between(Date(2020,5,10), DateTime(2021,5,15)) == 370
        @test days_between(DateTime(2021,5,15), Date(2020,5,10)) == -370
    end

    @testset "fractionalyears_between()" begin
        @test fractionalyears_between(Date(2000,1,1), Date(2000,1,2)) ≈ (length(Date(2000,1,1):Date(2000,1,2))-1)/365.25
        @test fractionalyears_between(Date(2000,1,1), Date(2001,7,1)) ≈ (length(Date(2000,1,1):Date(2001,7,1))-1)/365.25
        @test fractionalyears_between(Date(2001,1,2), Date(2000,1,1)) ≈ -(length(Date(2001,1,2):-Day(1):Date(2000,1,1))-1)/365.25
    end

    @testset "years_between()" begin
        @test years_between(Date(2000, 5, 10), DateTime(2020, 5, 9)) == 19
        @test years_between(DateTime(2020, 5, 9), Date(2000, 5, 10)) == -19
        @test years_between(Date(2000, 5, 10), Date(2020, 5, 10)) == 20
        @test years_between(DateTime(2020, 5, 10), Date(2000, 5, 10)) == -20

        # Same year
        @test years_between(Date(2020, 12, 31), DateTime(2021, 1, 1)) == 0
        @test years_between(DateTime(2021, 1, 1), Date(2020, 12, 31)) == 0
    end

    @testset "creatinine_mgdl()" begin

        import KidneyAllocation.creatinine_mgdl
       
        @test creatinine_mgdl(88.4) ≈ 1.0 
        
        # very high values (warning expected)
        @test creatinine_mgdl(1000) ≈ (1000/88.4)

        @test_throws AssertionError creatinine_mgdl(-5)

    end

end