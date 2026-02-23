using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: load_donor
import KidneyAllocation: evaluate_kdri, parse_abo, creatinine_mgdl

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
df = load_donor(donor_filepath)






"""
    build_donor_registry(filepath::String) -> Vector{Donor}

Construct a registry of `Donor` objects from a CSV file containing donor information.

### Details
This function loads donor data using [`load_donor`](@ref) and applies the following steps:
- Removes donors with missing information required for donor characterization and KDRI computation.
- Computes the Kidney Donor Risk Index (KDRI) for each donor.

### Notes
- The donor arrival date is derived from the donor death date (`DON_DEATH_TM`).
- Rows with incomplete data required for KDRI computation are discarded.

### Returns
A vector of `Donor` objects for each valid donor record.
"""
function build_donor_registry(filepath::String)
    df = load_donor(filepath)

    kdri_required = Symbol[
        :DON_ID, :DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD
    ]
    donor_required = Symbol[
        :DON_DEATH_TM, :DON_AGE, :DON_BLOOD, :DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2
    ]

    dropmissing!(df, kdri_required)
    dropmissing!(df, donor_required)

    G = groupby(df, :DON_ID)
    donors = Vector{Donor}(undef, length(G))

    for (i, g) in enumerate(G)
        r = first(g)

        age = r.DON_AGE
        height = r.HEIGHT
        weight = r.WEIGHT
        hypertension = r.HYPERTENSION == 1
        diabetes = r.DIABETES == 1
        cva = r.DEATH ∈ (4, 16)
        creatinine = creatinine_mgdl(r.CREATININE)
        dcd = r.DCD == 1

        kdri = evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

        arrival = Date(r.DON_DEATH_TM)
        blood = parse_abo(r.DON_BLOOD)

        donors[i] = Donor(arrival, age, blood, r.DON_A1, r.DON_A2, r.DON_B1, r.DON_B2, r.DON_DR1, r.DON_DR2, kdri)
    end

    return donors
end




@time donors = build_donor_registry(donor_filepath)