"""
Enumeration of ABO blood groups (without Rh factor).

# Values
- `O`: Blood group O.
- `A`: Blood group A.
- `B`: Blood group B.
- `AB`: Blood group AB.
"""
@enum ABOGroup O A B AB

"""
    is_abo_compatible(d::ABOGroup, r::ABOGroup)::Bool

Returns `true` if a donor with blood group `d` is ABO-compatible with
a recipient `r` under these kidney allocation rules:

- O → O only
- A → A or AB
- B → B or AB
- AB → AB only
"""
function is_abo_compatible(d::ABOGroup, r::ABOGroup)::Bool
    if d == O
        return r == O
    elseif d == A
        return r == A || r == AB
    elseif d == B
        return r == B || r == AB
    elseif d == AB
        return r == AB
    else
        return false
    end
end



"""
    parse_abo(s::AbstractString) -> ABOGroup

Parse an ABO blood group from a string like "O", "A", "B", or "AB".
"""
function parse_abo(s::AbstractString)::ABOGroup
    t = uppercase(strip(s))
    if t == "O"
        return O
    elseif t == "A"
        return A
    elseif t == "B"
        return B
    elseif t == "AB"
        return AB
    else
        throw(ArgumentError("Invalid ABO group string: '$s' (expected O, A, B, AB)"))
    end
end