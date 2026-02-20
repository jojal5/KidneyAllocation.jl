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

include("TransplantEntity/Donor.jl")
include("TransplantEntity/Recipient.jl")


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
    get_bloodtype(t::TransplantEntity)

Get blood type of transplant entity `t`.
"""
function get_bloodtype(t::TransplantEntity)
    return t.blood
end

"""
    get_HLA(t::TransplantEntity)

Get HLA of transplant entity `t` (a1, a2, b1, b2, dr1, dr2).
"""
function get_HLA(t::TransplantEntity)
    return (t.a1, t.a2, t.b1, t.b2, t.dr1, t.dr2)
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
    is_abo_compatible(donor::Donor, recipient::Recipient)

Returns `true` if the donor is ABO-compatible with the recipient under these kidney allocation rules.

## Details

See @ref(is_abo_compatible(d::ABOGroup, r::ABOGroup)) for ABO-compatibility rules.
"""
function is_abo_compatible(donor::Donor, recipient::Recipient)
    return is_abo_compatible(donor.blood, recipient.blood)
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

