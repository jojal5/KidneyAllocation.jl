"""
    calculate_epts(recipient::Recipient, has_diabetes::Bool=false, prior_transplant::Bool=false, current_date::DateTime=Dates.now())

Calculate the Estimated Post Transplant Survival (EPTS) score estimation for a kidney transplant candidate.

## Details
- `recipient::Recipient`: The recipient for whom to calculate the EPTS score
- `has_diabetes::Bool`: Whether the recepient has diabetes (defaults to false)
- `prior_transplant::Bool`: Whether the recipient had a previous solid organ transplant (defaults to false)
- `current_date::DateTime`: The reference date for calculating age and dialysis time (defaults to current date)

Convert raw score to final EPTS score (0-100 scale)

Note: This is a simplified conversion. In practice, this would be calculated 
based on the distribution of raw scores in the recipient pool.
Typically done via lookup tables maintained by OPTN.
"""
function calculate_epts(recipient::Recipient, has_diabetes::Bool=false, prior_transplant::Bool=false, current_date::Date=Dates.today())

    current_age = years_between(recipient.birth, current_date)

    # Calculate years on dialysis
    years_on_dialysis = fractionalyears_between(recipient.dialysis, current_date)

    # Calculate raw score based on the EPTS formula
    # Formula based on OPTN/UNOS guidelines
    raw_score = 0.047 * max(current_age - 25, 0) + 0.237 * min(years_on_dialysis, 10)
    if has_diabetes
        raw_score += 0.315
    end
    if prior_transplant
        raw_score += 0.398
    end

    # For now, we'll use a simple exponential mapping to the 0-100 range
    epts_score = (1. - exp(-raw_score * 1.5)) * 100.

    return epts_score
end

"""
    evaluate_kdri(age, height, weight, ethnicity, hypertension, diabetes,
                  death, creatinine, hcv, dcd) -> Float64

Calculate the Kidney Donor Risk Index (KDRI) based on donor characteristics following the SRTR/OPTN formulation.

# Details

Based on the updated KDRI formulation.

## Arguments
- `age::Real`: Donor's age in years.
- `height::Real`: Donor's height in cm.
- `weight::Real`: Donor's weight in kg.
- `hypertension::Bool`: `true` if donor had a history of hypertension, otherwise `false`.
- `diabetes::Bool`: `true` if donor had a history of diabetes, otherwise `false`.
- `cva::Bool`: `true` if cause of death is cerebrovascular accident, otherwise `false`.
- `creatinine::Real`: Serum creatinine level (mg/dL).
- `dcd::Bool`: `true` for donation after circulatory death, otherwise `false`.

# Reference
https://optn.transplant.hrsa.gov/media/0zamk0dr/mac_kdpi_board-briefing-paper.pdf

"""
function evaluate_kdri(age::Real, height::Real, weight::Real,
                       hypertension::Bool, diabetes::Bool, cva::Bool,
                       creatinine::Real, dcd::Bool)::Float64

    @assert age > 0       "Age should be positive, got $age"
    @assert 1 < height < 241    "Height should be in (1, 241) cm, got $height"
    @assert 0.9072 < weight < 294.835  "Weight should be in (0.9072, 294.835) kg, got $weight"
    @assert 0.01 < creatinine < 9.99   "Creatinine should be in (0.01, 9.99) mg/dL, got $creatinine"

    if creatinine > 8
        @warn "Creatinine value entered ($creatinine) exceeds 8 mg/dL and is capped at 707 in KDRI calculation."
    end

    if age > 100
        @warn "Age value entered ($age) is atypical."
    end

    # log-KDRI
    lkdri = 0.0

    # age effect
    lkdri += 0.0092 * (age - 40)
    lkdri += ifelse(age < 18, 0.0113 * (age - 18), 0.0)
    lkdri += ifelse(age > 50, 0.0067 * (age - 50), 0.0)

    # height effect
    lkdri += -0.0557 * (height - 170) / 10

    # weight effect
    lkdri += ifelse(weight < 80, -0.0333 * (weight - 80) / 5.0, 0.0)

    # hypertension effect
    lkdri += ifelse(hypertension, 0.1106, 0.0)

    # diabetes effect
    lkdri += ifelse(diabetes, 0.2577, 0.0)

    # cause of death effect (e.g., cerebrovascular accident)
    lkdri += ifelse(cva, 0.0743, 0.0)

    # creatinine effect (creatinine in µmol/L, capped at 8 mg/dL equivalent)
    capped_creatinine = min(creatinine, 8.)
    lkdri += 0.2128 * (capped_creatinine - 1.0)
    lkdri += ifelse(capped_creatinine > 1.5, -0.2199 * (capped_creatinine - 1.5), 0.0)

    # DCD effect
    lkdri += ifelse(dcd, 0.1966, 0.0)

    return exp(lkdri)
end
