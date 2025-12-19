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

# Fields
- `arrival::Date`: Arrival date of the donor.
- `age::Int64`: Age of the donor at the time of donation.
- `blood::ABOGroup`: ABO blood group of the donor.
- `a1::HLA`, `a2::HLA`: HLA-A antigens.
- `b1::HLA`, `b2::HLA`: HLA-B antigens.
- `dr1::HLA`, `dr2::HLA`: HLA-DR antigens.
- `kdri::Float64`: Kidney Donor Risk Index.
"""
struct Donor <: TransplantEntity
    arrival::Date
    age::Int64
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    kdri::Float64

    function Donor(arrival::Date,
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
    age::Int,
    blood::ABOGroup,
    a1::Integer, a2::Integer,
    b1::Integer, b2::Integer,
    dr1::Integer, dr2::Integer,
    kdri::Float64)

    return Donor(Date(arrival),
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

Represents a recipient in the kidney transplantation system (internal use).

# Fields
- `birth::Date`: Birth date of the recipient.
- `dialysis::Date`: Start date of dialysis.
- `arrival::Date`: Date the recipient was added to the waitlist.
- `blood::ABOGroup`: ABO blood type of the recipient.
- `a1::HLA`, `a2::HLA`: HLA-A antigens.
- `b1::HLA`, `b2::HLA`: HLA-B antigens.
- `dr1::HLA`, `dr2::HLA`: HLA-DR antigens.
- `cpra::Int64`: Calculated Panel Reactive Antibody (0–100).
- `expiration_date::Union{Date,Nothing}`: Eligibility expiration date, or `nothing`.
"""
struct Recipient <: TransplantEntity
    birth::Date
    dialysis::Date
    arrival::Date
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    cpra::Int64
    expiration_date::Union{Date,Nothing}

    function Recipient(birth::Date,
        dialysis::Date,
        arrival::Date,
        blood::ABOGroup,
        a1::HLA, a2::HLA,
        b1::HLA, b2::HLA,
        dr1::HLA, dr2::HLA,
        cpra::Int64;
        expiration_date::Union{Date,Nothing}=nothing)

        a1 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a1 = $a1"))
        a2 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a2 = $a2"))
        b1 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b1 = $b1"))
        b2 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b2 = $b2"))
        dr1 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr1 = $dr1"))
        dr2 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr2 = $dr2"))

        (0 <= cpra <= 100) || throw(ArgumentError("cpra must be in [0, 100], got $cpra"))

        return new(birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra, expiration_date)
    end
end

# Outer constructors (mixed Date/DateTime)

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::HLA, a2::HLA,
    b1::HLA, b2::HLA,
    dr1::HLA, dr2::HLA,
    cpra::Integer;
    expiration_date::Union{Date,DateTime,Nothing}=nothing)

    return Recipient(Date(birth), Date(dialysis), Date(arrival),
        blood,
        a1, a2, b1, b2, dr1, dr2,
        Int64(cpra);
        expiration_date=expiration_date === nothing ? nothing : Date(expiration_date))
end

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::Integer, a2::Integer,
    b1::Integer, b2::Integer,
    dr1::Integer, dr2::Integer,
    cpra::Integer;
    expiration_date::Union{Date,DateTime,Nothing}=nothing)

    return Recipient(Date(birth), Date(dialysis), Date(arrival),
        blood,
        HLA(a1), HLA(a2),
        HLA(b1), HLA(b2),
        HLA(dr1), HLA(dr2),
        Int64(cpra);
        expiration_date=expiration_date === nothing ? nothing : Date(expiration_date))
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
        "  CPRA           : $(r.cpra)\n",
        "  Expiration     : ",
        r.expiration_date === nothing ? "none" : string(r.expiration_date)
    )
end

function Base.summary(io::IO, r::Recipient)
    print(io,
        "Recipient(blood=$(r.blood), CPRA=$(r.cpra), ",
        "A=($(r.a1),$(r.a2)), ",
        "B=($(r.b1),$(r.b2)), ",
        "DR=($(r.dr1),$(r.dr2)))"
    )
end

"""
    shift_recipient_timeline(recipient::Recipient, new_arrival::Date) -> Recipient

Shift the recipient's `birth`, `dialysis`, and `expiration_date` so that the
recipient's timeline is consistent with a new `arrival` date.

## Details

All dates are shifted by the same number of days: the difference between the
current `recipient.arrival` and the new `arrival`. This preserves age and
waiting-time durations relative to the (new) arrival date.
"""
function shift_recipient_timeline(recipient::Recipient, new_arrival::Date)::Recipient
    # Signed shift (in days) from old arrival to new arrival
    shift_days = days_between(recipient.arrival, new_arrival)

    birth    = recipient.birth + Day(shift_days)
    dialysis = recipient.dialysis + Day(shift_days)

    expiration_date = recipient.expiration_date === nothing ? nothing :
                      recipient.expiration_date + Day(shift_days)

    return Recipient(birth,
                     dialysis,
                     new_arrival,
                     recipient.blood,
                     recipient.a1, recipient.a2,
                     recipient.b1, recipient.b2,
                     recipient.dr1, recipient.dr2,
                     recipient.cpra;
                     expiration_date = expiration_date)
end

"""
    has_expiration(r::Recipient)::Bool

Returns `true` if the recipient has a defined expiration date,
and `false` if `expiration_date` is `nothing`.
"""
has_expiration(r::Recipient) = r.expiration_date !== nothing

"""
    is_expired(r::Recipient, t::Date) -> Bool

Returns `true` if the recipient has an expiration date and it is
strictly earlier than the given date `t`. Returns `false` if
there is no expiration date.
"""
is_expired(r::Recipient, t::Date) =
    has_expiration(r) && r.expiration_date < t

"""
    is_active(r::Recipient, t::Date) -> Bool

Returns `true` if the recipient is active on the waitlist at time `t`.

A recipient is considered active if:
- the current time `t` is on or after their arrival date, and
- they have no expiration date, or the expiration date is on or after `t`.
"""
is_active(r::Recipient, t::Date) =
    t >= r.arrival && !is_expired(r, t)

"""
    is_hetero(t::TransplantEntity)

Check if the entity is heterozygous at the HLA-DR locus.

## Details 
Returns true if `dr1` and `dr2` are different; otherwise return `false`.
"""
is_hetero(t::TransplantEntity) = t.dr1 != t.dr2


