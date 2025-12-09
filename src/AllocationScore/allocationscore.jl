
"""
    dr_mismatch_count(d::Donor, r::Recipient)

Compute the number of HLA-DR mismatches between two the donor`d` and the recipient `r`.

## Details

At the DR locus, each entity has two alleles (`dr1`, `dr2`). A mismatch is
counted when a donor DR allele is not present in either of the recipient's
DR alleles. The maximum number of DR mismatches is 2.
"""
function dr_mismatch_count(d::TransplantEntity, r::TransplantEntity)::Int
    mm = 0

    # First DR allele
    if d.dr1 != r.dr1 && d.dr1 != r.dr2
        mm += 1
    end

    # Second DR allele (avoid double-counting if equal to dr1)
    if d.dr2 != d.dr1 && d.dr2 != r.dr1 && d.dr2 != r.dr2
        mm += 1
    end

    return mm
end

"""
    has_identical_HLA(d::Donor, r::Recipient)

Returns `true` if donor and recipient have identical HLA-A, HLA-B,
and HLA-DR alleles.

# Returns
- `Bool`: `true` if all A, B, and DR alleles match exactly.
"""
function has_identical_HLA(d::Donor, r::Recipient)::Bool
    return (d.a1 == r.a1 && d.a2 == r.a2 &&
            d.b1 == r.b1 && d.b2 == r.b2 &&
            d.dr1 == r.dr1 && d.dr2 == r.dr2)
end

"""
    score(d::Donor, r::Recipient) -> Float64

Compute the overall compatibility score between a donor and a recipient.

## Details
The score is computed by aggregating the individual component scores:

- wait time on dialysis
- HLA-DR mismatches
- cPRA (sensitization)
- age gap
- recipient age at donor arrival

The final score is rounded to two decimal places.
"""
function score(d::Donor, r::Recipient)::Float64
    total = 0.0
    total += score_wait_time(d, r)
    total += score_mis_dr(d, r)
    total += score_CPRA(r)
    total += score_age_gap(d, r)
    total += score_age(d, r)

    return round(total, digits = 2)
end

"""
    score_age(d::Donor, r::Recipient)

Compute the score based on the recipient's age at the donor's arrival date.
"""
function score_age(d::Donor, r::Recipient)
    recipient_age = years_between(r.birth, d.arrival)

    @assert recipient_age > 0 "The recipient age $recipient_age should be greater than 0."
    recipient_age = max(1, recipient_age)

    return 50.0 / recipient_age
end


"""
    score_wait_time(d::Donor, r::Recipient)

Compute the allocation score based on the recipient's wait time on dialysis.

## Details 
Wait time is measured in full calendar years between the dialysis start date (`r.dialysis`) and the donor's arrival date (`d.arrival`).

The scoring scale is defined by yearly categories:

    years:   0   1   2   3   4   5   6   7   8   9   ≥10
    score:   0 0.5   1   2   4   6   8  10  12  14   18
"""
function score_wait_time(d::Donor, r::Recipient)
    scores = [0, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 18]

    # Full years on dialysis before donor arrival
    wait_time = years_between(r.dialysis, d.arrival)

    # Cap at highest category
    if wait_time >= length(scores)
        return scores[end]
    else
        return scores[wait_time + 1]
    end
end

"""
    score_mis_dr(d::Donor, r::Recipient)

Evaluate the score based on HLA-DR mismatches between a donor and a recipient.

Scoring rules:
- If donor and recipient are identical at all HLA-A, HLA-B, and HLA-DR locus
  (`has_identical_HLA(d, r) == true`), the score is `8`.
- Otherwise, the score depends on:
    * the number of DR mismatches (`0`, `1`, or `2`), given by
      `dr_mismatch_count(d, r)`, and
    * whether the recipient is heterozygous at DR (`is_hetero(r)`).

Score table:

- `is_hetero(r) == true` (heterozygous recipient):
    DR mismatches:   0    1    2
    score:           4    1    0

- `is_hetero(r) == false` (homozygous recipient):
    DR mismatches:   0    1    2
    score:           4    4    0
"""
function score_mis_dr(d::Donor, r::Recipient)::Int
    # Special case: donor and recipient have identical A, B, and DR alleles
    if has_identical_HLA(d, r)
        return 8
    end

    # rows: hetero (true), homo (false)
    # columns: DR mismatches = 0, 1, 2
    scores = ((4, 1, 0),   # heterozygous recipient
              (4, 4, 0))   # homozygous recipient

    mis_dr = dr_mismatch_count(d, r)

    type_idx = is_hetero(r) ? 1 : 2    # 1 = hetero row, 2 = homo row
    score_idx = mis_dr + 1            

    return scores[type_idx][score_idx]
end

"""
    score_age_gap(d::Donor, r::Recipient)

Compute the score based on the age gap between donor and recipient at the time of donor arrival.

## Details

Scoring rules:
- gap ≤ 10 years  → score = 4
- 10 < gap ≤ 20   → score = 2
- gap > 20        → score = 0
"""
function score_age_gap(d::Donor, r::Recipient)::Int

    donor_age = d.age

    # Recipient age at donor arrival
    recipient_age = years_between(r.birth, d.arrival)

    # Absolute age difference
    gap = abs(donor_age - recipient_age)

    if gap <= 10
        return 4
    elseif  gap <= 20
        return 2
    else
        return 0
    end
end


"""
    score_CPRA(r::Recipient)

Evaluate the score based on the recipient's cPRA value.

Scoring rules:
- 0 ≤ cPRA < 20   → score = 0
- 20 ≤ cPRA < 80  → score = 3
- 80 ≤ cPRA ≤ 100 → score = 8
"""
function score_CPRA(r::Recipient)::Int
    cpra = r.CPRA

    if cpra < 20
        return 0
    elseif cpra < 80
        return 3
    else
        return 8
    end
end



