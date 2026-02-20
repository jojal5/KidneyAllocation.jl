@testset "Donor constructor" begin
    # A valid arrival time
    arrival = DateTime(2025, 1, 1, 12, 0, 0)

    # Choose valid HLA alleles from your allowed sets
    dr1 = HLA(1)
    dr2 = HLA(4)
    a1 = HLA(24)
    a2 = HLA(26)
    b1 = HLA(44)
    b2 = HLA(51)

    age = 45
    kdri = 1.5
    blood = A

    # --- Valid construction ---
    d = Donor(arrival, age, blood,
        a1, a2, b1, b2,
        dr1, dr2, kdri)

    @test d isa Donor
    @test d.arrival == Date(arrival)
    @test d.age == age
    @test d.blood == blood
    @test d.a1 == a1
    @test d.a2 == a2
    @test d.b1 == b1
    @test d.b2 == b2
    @test d.dr1 == dr1
    @test d.dr2 == dr2
    @test d.kdri == kdri

    # --- Invalid DR alleles ---
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, b1, b2,
        HLA(99), dr2, kdri)
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, b1, b2,
        dr1, HLA(99), kdri)

    # --- Invalid A alleles ---
    @test_throws ArgumentError Donor(arrival, age, blood,
        HLA(0), a2, b1, b2,
        dr1, dr2, kdri)
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, HLA(0), b1, b2,
        dr1, dr2, kdri)

    # --- Invalid B alleles ---
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, HLA(0), b2,
        dr1, dr2, kdri)
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, b1, HLA(0),
        dr1, dr2, kdri)

    # --- Invalid age ---
    @test_throws ArgumentError Donor(arrival, -1, blood,
        a1, a2, b1, b2,
        dr1, dr2, kdri)

    # --- Invalid KDRI ---
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, b1, b2,
        dr1, dr2, 0.0)
    @test_throws ArgumentError Donor(arrival, age, blood,
        a1, a2, b1, b2,
        dr1, dr2, -1.0)
end


@testset "set_donor_arrival" begin

    import KidneyAllocation: set_donor_arrival
    # Construct a reference donor
    arrival = Date(2020, 1, 1)
    donor = Donor(arrival,
        45,
        A,
        HLA(24), HLA(26),
        HLA(44), HLA(52),
        HLA(7), HLA(15),
        1.23)

    new_arrival = Date(2025, 1, 1)

    donor2 = set_donor_arrival(donor, new_arrival)

    @test donor2.arrival == new_arrival

    # All other fields must remain unchanged
    @test donor2.age == donor.age
    @test donor2.blood == donor.blood

    @test donor2.a1 == donor.a1
    @test donor2.a2 == donor.a2
    @test donor2.b1 == donor.b1
    @test donor2.b2 == donor.b2
    @test donor2.dr1 == donor.dr1
    @test donor2.dr2 == donor.dr2

    @test donor2.kdri == donor.kdri
end