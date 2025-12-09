
# Normalize Date / DateTime to DateTime
_dt(x::DateTime) = x
_dt(x::Date) = DateTime(x)


"""
    days_between(d1::Union{Date,DateTime}, d2::Union{Date,DateTime})

Compute the number of days between two dates.

## Details

- Returns a positive integer if `d2` occurs after `d1`
- Returns a negative integer if `d2` occurs before `d1`
- Returns 0 if the same calendar day

Both `Date` and `DateTime` inputs are accepted and normalized.
The result is symmetric:

    days_between(a, b) == -days_between(b, a)
"""
function days_between(d1::Union{Date,DateTime}, d2::Union{Date,DateTime})
    t1 = _dt(d1)
    t2 = _dt(d2)
    return (t2 - t1)/Day(1)
end

"""
    fractionalyears_between(d1, d2)

Compute the fractional number of years between two dates.

## Details

The calculation is based on the total number of days divided by 365.25,
providing an approximation that accounts for leap years on average.

Returns a positive value if `d2 > d1`, negative if `d2 < d1`,
and zero if both dates fall on the same day.
"""
function fractionalyears_between(d1::Union{Date,DateTime}, d2::Union{Date,DateTime})
    ndays = days_between(d1, d2)
    return ndays/365.25
end

"""
    years_between(d1::Union{Date,DateTime}, d2::Union{Date,DateTime})

Compute the number of full calendar years between two dates.

## Details 
- Returns a positive integer when `d2` is later than `d1`
- Returns a negative integer when `d2` is earlier than `d1`

Full years are counted only when the anniversary of the earlier date
has occurred relative to the later date.

This function satisfies the symmetry:
    years_between(a, b) == -years_between(b, a)
"""
function years_between(d1::Union{Date,DateTime}, d2::Union{Date,DateTime})::Int

    # Normalize inputs
    t1 = _dt(d1)
    t2 = _dt(d2)

    # Equal dates → zero
    t1 == t2 && return 0

    # Determine direction
    forward = t2 > t1

    # Let a = earlier, b = later
    a, b = forward ? (t1, t2) : (t2, t1)

    # Raw year difference
    y = year(b) - year(a)

    # Anniversary check:
    # Subtract 1 full year if the anniversary has NOT yet occurred in year(b)
    if dayofyear(b)<dayofyear(a)
        y -= 1
    end

    return forward ? y : -y
end

"""
    creatinine_mgdl(creat_umol_L::Real)

Convert creatinine from µmol/L to mg/dL.
"""
function creatinine_mgdl(creat_umol_L::Real)
    @assert creat_umol_L > 0 "Creatinine level should be positive, got $creat_umol_L."

    if !(40 < creat_umol_L < 900)
        @warn "Creatinine value ($creat_umol_L µmol/L) is atypical." maxlog=1
    end

    return creat_umol_L / 88.4
end


