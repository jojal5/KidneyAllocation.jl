
@testset "ABOGroup.jl" begin

    @testset "ABOGroup Enum Constructor" begin

        # --- Test that enum values exist and have correct type ---
        @test O isa ABOGroup
        @test A isa ABOGroup
        @test B isa ABOGroup
        @test AB isa ABOGroup

        # --- Test that constructor from Int works as expected ---
        @test ABOGroup(0) == O
        @test ABOGroup(1) == A
        @test ABOGroup(2) == B
        @test ABOGroup(3) == AB

        # --- Test that Int conversions back to Int work ---
        @test Int(O) == 0
        @test Int(A) == 1
        @test Int(B) == 2
        @test Int(AB) == 3

        @test_throws ArgumentError ABOGroup(4)
        @test_throws ArgumentError ABOGroup(-1)

    end

    @testset "ABO Compatibility" begin

        import KidneyAllocation.is_abo_compatible

        # O donations
        @test is_abo_compatible(O, O) == true
        @test is_abo_compatible(O, A) == false
        @test is_abo_compatible(O, B) == false
        @test is_abo_compatible(O, AB) == false

        # A donations
        @test is_abo_compatible(A, A) == true
        @test is_abo_compatible(A, AB) == true
        @test is_abo_compatible(A, O) == false
        @test is_abo_compatible(A, B) == false

        # B donations
        @test is_abo_compatible(B, B) == true
        @test is_abo_compatible(B, AB) == true
        @test is_abo_compatible(B, O) == false
        @test is_abo_compatible(B, A) == false

        # AB donations
        @test is_abo_compatible(AB, AB) == true
        @test is_abo_compatible(AB, O) == false
        @test is_abo_compatible(AB, A) == false
        @test is_abo_compatible(AB, B) == false

    end

end