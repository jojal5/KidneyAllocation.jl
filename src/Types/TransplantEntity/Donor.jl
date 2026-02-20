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