abstract type TransplantEntity end

# For broadcasting
Base.broadcastable(x::TransplantEntity) = Ref(x)

const HLA = UInt16


# Valid HLA-A serotypes
const VALID_HLA_A = Set{HLA}(HLA[
    1, 2, 3, 11, 23, 24, 25, 26,
    29, 30, 31, 32, 33, 34, 36,
    66, 68, 69, 74, 80,
    203, 2403, 3401, 6601
])

# Valid HLA-B serotypes
const VALID_HLA_B = Set{HLA}(HLA[
    7, 8, 13, 14, 18, 27, 35,
    37, 38, 39, 41, 42,
    44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58,
    60, 61, 62, 63, 64, 65, 67,
    70, 71, 72, 73, 75, 76, 77, 78,
    81, 82,
    3901, 4402, 4403, 5102, 8201
])

# Valid HLA-DRB1 serotypes
const VALID_HLA_DR = Set{HLA}(HLA[
    1, 3, 4, 7, 8, 9,
    10, 11, 12, 13, 14, 15, 16, 17, 18,
    103, 1404
])

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
        age ≥ 0 || throw(ArgumentError("Donor age must be ≥ 0, got $age"))
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
    set_donor_arrival(donor::Donor, new_arrival::Date) -> Donor

Return a copy of `donor` with the arrival date replaced by `new_arrival`.

## Details

All other donor attributes (age, blood group, HLA antigens, KDRI) are preserved.
"""
function set_donor_arrival(donor::Donor, new_arrival::Date)::Donor
    return Donor(new_arrival,
                 donor.age,
                 donor.blood,
                 donor.a1, donor.a2,
                 donor.b1, donor.b2,
                 donor.dr1, donor.dr2,
                 donor.kdri)
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

    birth = recipient.birth + Day(shift_days)
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
        expiration_date=expiration_date)
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

"""
    get_arrival(t::TransplantEntity)

Get arrival time of transplant entity `t`.
"""
function get_arrival(t::TransplantEntity)
    return t.arrival
end

"""
    get_HLA(t::TransplantEntity)

Get HLA of transplant entity `t` (a1, a2, b1, b2, dr1, dr2).
"""
function get_HLA(t::TransplantEntity)
    return (t.a1, t.a2, t.b1, t.b2, t.dr1, t.dr2)
end

"""
    get_bloodtype(t::TransplantEntity)

Get blood type of transplant entity `t`.
"""
function get_bloodtype(t::TransplantEntity)
    return t.blood
end

"""
    is_abo_compatible(donor::Donor, recipient::Recipient)

Returns `true` if the donor is ABO-compatible with the recipient under these kidney allocation rules.

## Details

See @ref(is_abo_compatible(d::ABOGroup, r::ABOGroup)) for ABO-compatibility rules.
"""
function is_abo_compatible(donor::Donor, recipient::Recipient)
    return is_abo_compatible(donor.blood, recipient.blood)
end


"""
    get_HLA_A(t::TransplantEntity) -> Tuple

Return the pair of HLA-A alleles of `t`.
"""
function get_HLA_A(t::TransplantEntity)
    return (t.a1, t.a2)
end

"""
    get_HLA_B(t::TransplantEntity) -> Tuple

Return the pair of HLA-B alleles of `t`.
"""
function get_HLA_B(t::TransplantEntity)
    return (t.b1, t.b2)
end

"""
    get_HLA_DR(t::TransplantEntity) -> Tuple

Return the pair of HLA-DR alleles of `t`.
"""
function get_HLA_DR(t::TransplantEntity)
    return (t.dr1, t.dr2)
end

"""
    mismatch_locus(l₁, l₂) -> Int

Return the number of allele mismatches between two HLA loci.

Counts how many alleles in `l₂` are not present in `l₁`.
"""
function mismatch_locus(l₁::Tuple{HLA, HLA}, l₂::Tuple{HLA, HLA})::Int
    return count(x -> !(x in l₁), l₂)
end

"""
    mismatch_count(t₁, t₂) -> Int

Return the total number of HLA mismatches between two transplant entities across loci A, B, and DR.
"""
function mismatch_count(t₁::TransplantEntity, t₂::TransplantEntity)::Int
    mm_A  = mismatch_locus(get_HLA_A(t₁),  get_HLA_A(t₂))
    mm_B  = mismatch_locus(get_HLA_B(t₁),  get_HLA_B(t₂))
    mm_DR = mismatch_locus(get_HLA_DR(t₁), get_HLA_DR(t₂))

    return mm_A + mm_B + mm_DR
end



function Base.show(io::IO, ::MIME"text/plain", donors::AbstractVector{<:Donor})
    n = length(donors)
    println(io, "$n-element Vector{Donor}:")

    n == 0 && return

    limit = get(io, :limit, false)
    rows, _ = displaysize(io)
    max_lines = max(rows - 1, 1)  # lines available after the header

    # Decide how many elements to show
    if !limit || n <= max_lines
        # Show everything
        for d in donors
            print(io, " ")
            summary(io, d)
            println(io)
        end
        return
    end

    # Long vector display: show head, ellipsis, tail
    # Reserve 1 line for " ⋮"
    avail = max_lines - 1
    head = max(1, avail ÷ 2)
    tail = max(1, avail - head)

    # Ensure we don't accidentally show all elements when n is small-ish
    if head + tail >= n
        head = min(head, n)
        tail = 0
    end

    for i in 1:head
        print(io, " ")
        summary(io, donors[i])
        println(io)
    end

    println(io, " ⋮")

    for i in (n-tail+1):n
        print(io, " ")
        summary(io, donors[i])
        println(io)
    end
end

function Base.show(io::IO, ::MIME"text/plain", recipients::AbstractVector{<:Recipient})
    n = length(recipients)
    println(io, "$n-element Vector{Recipient}:")

    n == 0 && return

    limit = get(io, :limit, false)
    rows, _ = displaysize(io)
    max_lines = max(rows - 1, 1)  # available lines after the header

    if !limit || n <= max_lines
        # Show everything
        for r in recipients
            print(io, " ")
            summary(io, r)
            println(io)
        end
        return
    end

    # Long vector display: show head, ellipsis, tail
    avail = max_lines - 1          # reserve one line for " ⋮"
    head = max(1, avail ÷ 2)
    tail = max(1, avail - head)

    if head + tail >= n
        head = min(head, n)
        tail = 0
    end

    for i in 1:head
        print(io, " ")
        summary(io, recipients[i])
        println(io)
    end

    println(io, " ⋮")

    for i in (n-tail+1):n
        print(io, " ")
        summary(io, recipients[i])
        println(io)
    end
end