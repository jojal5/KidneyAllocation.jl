@testset "TransplantEntity.jl" begin

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

    @testset "Recipient constructor" begin
        birth = DateTime(1980, 1, 1)
        arrival = DateTime(2024, 1, 1)
        dialysis = DateTime(2015, 1, 1)

        # Valid HLA alleles taken from valid sets:
        dr1 = HLA(1)
        dr2 = HLA(4)
        a1 = HLA(24)
        a2 = HLA(26)
        b1 = HLA(44)
        b2 = HLA(51)

        cpra = 30
        blood = A  # ABOGroup.A

        # --- Valid construction with default expiration_date (nothing) ---
        r = Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            cpra)

        @test r isa Recipient
        @test r.birth == birth
        @test r.dialysis == dialysis
        @test r.arrival == arrival
        @test r.dr1 == dr1
        @test r.dr2 == dr2
        @test r.a1 == a1
        @test r.a2 == a2
        @test r.b1 == b1
        @test r.b2 == b2
        @test r.cpra == cpra
        @test r.blood == blood
        @test r.expiration_date === nothing

        # --- Valid construction with explicit expiration_date ---
        exp_date = DateTime(2026, 1, 1)
        r2 = Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            cpra; expiration_date=exp_date)

        @test r2.expiration_date == exp_date

        # --- CPRA boundary values: 0 and 100 should be accepted ---
        r0 = Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            0)
        @test r0.cpra == 0

        r100 = Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            100)
        @test r100.cpra == 100

        # --- Invalid CPRA values ---
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            -1)
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            120)

        # --- Invalid DR alleles ---
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            HLA(99), dr2,
            cpra)
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, HLA(99),
            cpra)

        # --- Invalid A alleles ---
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            HLA(0), a2, b1, b2,
            dr1, dr2,
            cpra)
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, HLA(0), b1, b2,
            dr1, dr2,
            cpra)

        # --- Invalid B alleles ---
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, HLA(0), b2,
            dr1, dr2,
            cpra)
        @test_throws ArgumentError Recipient(birth, dialysis, arrival,
            blood,
            a1, a2, b1, HLA(0),
            dr1, dr2,
            cpra)
    end

    @testset "Recipient expiration and activity helpers" begin

        import KidneyAllocation: has_expiration, is_expired, is_active

        birth = Date(1980, 1, 1)
        arrival = Date(2024, 1, 1)
        dialysis = Date(2015, 1, 1)

        # Valid HLA alleles from your allowed sets
        dr1 = HLA(1)
        dr2 = HLA(4)
        a1 = HLA(24)
        a2 = HLA(26)
        b1 = HLA(44)
        b2 = HLA(51)

        cpra = 30
        blood = A   # ABOGroup.A

        # --- Recipient without expiration date (expiration_date = nothing) ---
        r_noexp = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra)

        @test has_expiration(r_noexp) == false

        t_before = Date(2023, 12, 31)
        t_at_arrival = arrival
        t_after = Date(2025, 1, 1)

        # No expiration: never expired
        @test is_expired(r_noexp, t_before) == false
        @test is_expired(r_noexp, t_at_arrival) == false
        @test is_expired(r_noexp, t_after) == false

        # Activity for no-expiration recipient
        @test is_active(r_noexp, t_before) == false     # before arrival
        @test is_active(r_noexp, t_at_arrival) == true     # exactly at arrival
        @test is_active(r_noexp, t_after) == true      # after arrival, no expiration

        # --- Recipient with explicit expiration date ---
        exp_date = Date(2025, 1, 1)
        r_exp = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra; expiration_date=exp_date)

        @test has_expiration(r_exp) == true

        t_before_exp = Date(2024, 6, 1)
        t_at_exp = exp_date
        t_after_exp = Date(2025, 6, 1)

        # Expired logic: uses strict < t
        @test is_expired(r_exp, t_before_exp) == false
        @test is_expired(r_exp, t_at_exp) == false   # equal → NOT expired
        @test is_expired(r_exp, t_after_exp) == true    # after expiration → expired

        # Activity with expiration:
        # - before arrival → inactive
        # - between arrival and exp_date inclusive → active
        # - after expiration → inactive
        @test is_active(r_exp, t_before) == false
        @test is_active(r_exp, arrival) == true
        @test is_active(r_exp, t_before_exp) == true
        @test is_active(r_exp, t_at_exp) == true
        @test is_active(r_exp, t_after_exp) == false
    end

    @testset "shift_recipient_timeline" begin
        import KidneyAllocation: days_between, shift_recipient_timeline

        birth = Date(1980, 1, 1)
        dialysis = Date(2010, 1, 1)
        arrival = Date(2020, 1, 1)
        expiry = Date(2025, 1, 1)

        r = Recipient(
            birth,
            dialysis,
            arrival,
            A,
            24, 26,
            44, 51,
            1, 4,
            80;
            expiration_date=expiry
        )

        # --- Forward shift ---
        new_arrival = Date(2022, 1, 1)
        r2 = shift_recipient_timeline(r, new_arrival)

        shift_days = days_between(arrival, new_arrival)

        @test r2.arrival == new_arrival
        @test r2.birth == birth + Day(shift_days)
        @test r2.dialysis == dialysis + Day(shift_days)
        @test r2.expiration_date == expiry + Day(shift_days)

        # Preserve non-date fields
        @test r2.blood == r.blood
        @test r2.a1 == r.a1
        @test r2.b2 == r.b2
        @test r2.dr1 == r.dr1
        @test r2.cpra == r.cpra

        # --- Backward shift ---
        earlier_arrival = Date(2018, 1, 1)
        r3 = shift_recipient_timeline(r, earlier_arrival)

        shift_days_back = days_between(arrival, earlier_arrival)

        @test r3.arrival == earlier_arrival
        @test r3.birth == birth + Day(shift_days_back)
        @test r3.dialysis == dialysis + Day(shift_days_back)
        @test r3.expiration_date == expiry + Day(shift_days_back)

        # --- No expiration date ---
        r_noexp = Recipient(
            birth,
            dialysis,
            arrival,
            O,
            24, 26,
            44, 51,
            1, 4,
            10
        )

        r4 = shift_recipient_timeline(r_noexp, new_arrival)

        @test r4.expiration_date === nothing
        @test r4.birth == birth + Day(shift_days)
        @test r4.dialysis == dialysis + Day(shift_days)
    end


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

        import KidneyAllocation: get_arrival, get_HLA, get_bloodtype

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
        @test get_bloodtype(d) == A

        birth = Date(1980, 1, 1)
        dialysis = Date(2020,1,1)
        cpra = 0

        r = Recipient(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra)
        
        @test get_arrival(r) == arrival
        @test get_HLA(r) == (a1, a2, b1, b2, dr1, dr2)
        @test get_bloodtype(r) == A

    end

end
