@testset "TransplantEntity.jl" begin

    include("TransplantEntity/test_donor.jl")
    include("TransplantEntity/test_recipient.jl")

    @testset "is_hetero on Recipient and Donor" begin

        import KidneyAllocation.is_hetero

        birth = Date(1980, 1, 1)
        arrival = Date(2024, 1, 1)
        dialysis = Date(2015, 1, 1)

        # Heterozygous at DR: 1 and 4
        dr_het1 = HLA(1)
        dr_het2 = HLA(4)
        # Homozygous at DR: 1 and 1
        dr_hom1 = HLA(1)
        dr_hom2 = HLA(1)

        a1 = HLA(24)
        a2 = HLA(26)
        b1 = HLA(44)
        b2 = HLA(51)

        cpra = 0
        blood = A
        kdri = 1.2
        age = 40

        # Recipient heterozygous at DR
        r_het = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr_het1, dr_het2,
            cpra)

        # Recipient homozygous at DR
        r_hom = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr_hom1, dr_hom2,
            cpra)

        @test is_hetero(r_het) == true
        @test is_hetero(r_hom) == false

        # Donor heterozygous at DR
        d_het = Donor(arrival, age, blood,
            a1, a2, b1, b2,
            dr_het1, dr_het2,
            kdri)

        # Donor homozygous at DR
        d_hom = Donor(arrival, age, blood,
            a1, a2, b1, b2,
            dr_hom1, dr_hom2,
            kdri)

        @test is_hetero(d_het) == true
        @test is_hetero(d_hom) == false
    end

    @testset "property extraction" begin

        import KidneyAllocation: get_arrival, get_HLA, get_HLA_A, get_HLA_B, get_HLA_DR, get_bloodtype

        arrival = Date(2025, 1, 1)

        dr1 = HLA(1)
        dr2 = HLA(4)
        a1 = HLA(24)
        a2 = HLA(26)
        b1 = HLA(44)
        b2 = HLA(51)

        age = 45
        kdri = 1.5
        blood = A

        d = Donor(arrival, age, blood,
            a1, a2, b1, b2,
            dr1, dr2, kdri)

        @test get_arrival(d) == arrival
        @test get_HLA(d) == (a1, a2, b1, b2, dr1, dr2)
        @test get_HLA_A(d) == (a1, a2)
        @test get_HLA_B(d) == (b1, b2)
        @test get_HLA_DR(d) == (dr1, dr2)
        @test get_bloodtype(d) == A

        birth = Date(1980, 1, 1)
        dialysis = Date(2020, 1, 1)
        cpra = 0

        r = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra)

        @test get_arrival(r) == arrival
        @test get_HLA(r) == (a1, a2, b1, b2, dr1, dr2)
        @test get_bloodtype(r) == A

    end

    @testset "mismatch_count" begin
        import KidneyAllocation: get_HLA_A, get_HLA_B, get_HLA_DR, mismatch_locus, mismatch_count
        r = Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 73, 77, 15, 17, 0)
        d = Donor(Date(2000, 1, 1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5)

        @test mismatch_locus(get_HLA_A(r), get_HLA_A(d)) == 2
        @test mismatch_locus(get_HLA_B(r), get_HLA_B(d)) == 0
        @test mismatch_locus(get_HLA_DR(r), get_HLA_DR(d)) == 1

        @test mismatch_count(r, d) == 3
    end

end
