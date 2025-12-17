abstract type TransplantEntity end

const HLA = UInt16


# Valid HLA-A serotypes
df_HLA_A = CSV.read("src/types/valid_HLA-A.csv", DataFrame)
const VALID_HLA_A = Set{HLA}(df_HLA_A.VALID_HLA_A)

# Valid HLA-B serotypes
df_HLA_B = CSV.read("src/types/valid_HLA-B.csv", DataFrame)
const VALID_HLA_B = Set{HLA}(df_HLA_B.VALID_HLA_B)

# Valid HLA-DRB1 serotypes
df_HLA_DR = CSV.read("src/types/valid_HLA-DR.csv", DataFrame)
const VALID_HLA_DR = Set{HLA}(df_HLA_DR.VALID_HLA_DR)


"""
    Donor(...) <: TransplantEntity

Represents a donor with demographic, clinical, and immunologic attributes.

Fields
- `arrival::DateTime`: Arrival date and time of the donor.
- `age::Int64`: Age of the donor at the time of donation.
- `blood::ABOGroup`: ABO blood group of the donor.
- `a1::HLA`, `a2::HLA`: HLA-A antigens.
- `b1::HLA`, `b2::HLA`: HLA-B antigens.
- `dr1::HLA`, `dr2::HLA`: HLA-DR antigens.
- `kdri::Float64`: Kidney Donor Risk Index.
"""
struct Donor <: TransplantEntity
    arrival::DateTime
    age::Int64
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    kdri::Float64

    function Donor(arrival::DateTime,
        age::Int64,
        blood::ABOGroup,
        a1::HLA, a2::HLA,
        b1::HLA, b2::HLA,
        dr1::HLA, dr2::HLA,
        kdri::Float64)

        # Validate HLA alleles by locus
        a1 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a1 = $a1"))
        a2 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a2 = $a2"))

        b1 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b1 = $b1"))
        b2 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b2 = $b2"))

        dr1 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr1 = $dr1"))
        dr2 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr2 = $dr2"))

        # Validate age and kdri
        age > 0 || throw(ArgumentError("Donor age must be > 0, got $age"))
        kdri > 0 || throw(ArgumentError("KDRI must be > 0, got $kdri"))

        return new(arrival, age, blood,
            a1, a2, b1, b2,
            dr1, dr2, kdri)
    end
end

# Outer constructors

function Donor(arrival::Union{Date,DateTime},
               age::Int64,
               blood::ABOGroup,
               a1::Integer, a2::Integer,
               b1::Integer, b2::Integer,
               dr1::Integer, dr2::Integer,
               kdri::Float64)

    return Donor(_dt(arrival),
                 age,
                 blood,
                 HLA(a1), HLA(a2),
                 HLA(b1), HLA(b2),
                 HLA(dr1), HLA(dr2),
                 kdri)
end



function Base.show(io::IO, ::MIME"text/plain", d::Donor)
    print(io,
        "Donor\n",
        "  Arrival     : $(d.arrival)\n",
        "  Age         : $(d.age)\n",
        "  Blood Type  : $(d.blood)\n",
        "  HLA-A       : $(d.a1), $(d.a2)\n",
        "  HLA-B       : $(d.b1), $(d.b2)\n",
        "  HLA-DR      : $(d.dr1), $(d.dr2)\n",
        "  KDRI        : $(round(d.kdri, digits=4))"
    )
end

function Base.summary(io::IO, d::Donor)
    print(io,
        "Donor(age=$(d.age), blood=$(d.blood), ",
        "HLA-A=($(d.a1),$(d.a2)), ",
        "HLA-B=($(d.b1),$(d.b2)), ",
        "HLA-DR=($(d.dr1),$(d.dr2)), ",
        "KDRI=$(round(d.kdri, digits=3)))"
    )
end

"""
    Recipient(...) <: TransplantEntity

Represents a recipient in the kidney transplantation system.

This structure is intended for internal use only.

# Fields
- `birth::DateTime`: Birth date of the recipient.
- `dialysis::DateTime`: Start date of dialysis.
- `arrival::DateTime`: Date the recipient was added to the waitlist.
- `blood::ABOGroup`: ABO blood type of the recipient.
- `a1::HLA`, `a2::HLA`: HLA-A antigens.
- `b1::HLA`, `b2::HLA`: HLA-B antigens.
- `dr1::HLA`, `dr2::HLA`: HLA-DR antigens.
- `CPRA::Float64`: Calculated Panel Reactive Antibody (0–100).
- `expiration_date::Union{DateTime, Nothing}`: Date at which the
    recipient’s eligibility expires, or `nothing` if no expiration applies.
"""
struct Recipient <: TransplantEntity
    birth::DateTime
    dialysis::DateTime
    arrival::DateTime
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    CPRA::Int64
    expiration_date::Union{DateTime,Nothing}

    function Recipient(birth::DateTime,
        dialysis::DateTime,
        arrival::DateTime,
        blood::ABOGroup,
        a1::HLA, a2::HLA,
        b1::HLA, b2::HLA,
        dr1::HLA, dr2::HLA,
        CPRA::Int64;
        expiration_date::Union{DateTime,Nothing}=nothing)

        # Validate HLA alleles by locus
        a1 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a1 = $a1"))
        a2 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a2 = $a2"))

        b1 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b1 = $b1"))
        b2 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b2 = $b2"))

        dr1 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr1 = $dr1"))
        dr2 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr2 = $dr2"))

        # Validate CPRA range
        if CPRA < 0 || CPRA > 100
            throw(ArgumentError("CPRA must be in [0, 100], got $CPRA"))
        end

        return new(birth, dialysis, arrival,
            blood,
            a1, a2, b1, b2,
            dr1, dr2,
            CPRA, expiration_date)
    end
end

# Outer constructors

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::HLA, a2::HLA,
    b1::HLA, b2::HLA,
    dr1::HLA, dr2::HLA,
    CPRA::Int64;
    expiration_date::Union{Date,DateTime,Nothing}=nothing)

    return Recipient(_dt(birth),
        _dt(dialysis),
        _dt(arrival),
        blood,
        a1, a2, b1, b2,
        dr1, dr2,
        CPRA;
        expiration_date=expiration_date === nothing ? nothing : _dt(expiration_date))
end

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::Integer, a2::Integer,
    b1::Integer, b2::Integer,
    dr1::Integer, dr2::Integer,
    CPRA::Int64;
    expiration_date::Union{Date,DateTime,Nothing}=nothing)

    return Recipient(_dt(birth),
        _dt(dialysis),
        _dt(arrival),
        blood,
        HLA(a1), HLA(a2),
        HLA(b1), HLA(b2),
        HLA(dr1), HLA(dr2),
        CPRA;
        expiration_date=expiration_date === nothing ? nothing : _dt(expiration_date))
end

function Base.show(io::IO, ::MIME"text/plain", r::Recipient)
    print(io,
        "Recipient\n",
        "  Birth Date     : $(r.birth)\n",
        "  Dialysis Start : $(r.dialysis)\n",
        "  Arrival Date   : $(r.arrival)\n",
        "  Blood Type     : $(r.blood)\n",
        "  HLA-A          : $(r.a1), $(r.a2)\n",
        "  HLA-B          : $(r.b1), $(r.b2)\n",
        "  HLA-DR         : $(r.dr1), $(r.dr2)\n",
        "  CPRA           : $(r.CPRA)\n",
        "  Expiration     : ",
            r.expiration_date === nothing ? "none" : string(r.expiration_date)
    )
end

function Base.summary(io::IO, r::Recipient)
    print(io,
        "Recipient(blood=$(r.blood), CPRA=$(r.CPRA), ",
        "A=($(r.a1),$(r.a2)), ",
        "B=($(r.b1),$(r.b2)), ",
        "DR=($(r.dr1),$(r.dr2)))"
    )
end


"""
    has_expiration(r::Recipient)::Bool

Returns `true` if the recipient has a defined expiration date,
and `false` if `expiration_date` is `nothing`.
"""
has_expiration(r::Recipient) = r.expiration_date !== nothing

"""
    is_expired(r::Recipient, t::DateTime) -> Bool

Returns `true` if the recipient has an expiration date and it is
strictly earlier than the given time `t`. Returns `false` if
there is no expiration date.
"""
is_expired(r::Recipient, t::DateTime) =
    has_expiration(r) && r.expiration_date < t

"""
    is_active(r::Recipient, t::DateTime) -> Bool

Returns `true` if the recipient is active on the waitlist at time `t`.

A recipient is considered active if:
- the current time `t` is on or after their arrival date, and
- they have no expiration date, or the expiration date is on or after `t`.
"""
is_active(r::Recipient, t::DateTime) =
    t >= r.arrival && !is_expired(r, t)

"""
    is_hetero(t::TransplantEntity)

Check if the entity is heterozygous at the HLA-DR locus.

## Details 
Returns true if `dr1` and `dr2` are different; otherwise return `false`.
"""
is_hetero(t::TransplantEntity) = t.dr1 != t.dr2


