using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation


import KidneyAllocation: build_last_cpra_registry, load_recipient, infer_recipient_expiration_date, parse_abo

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

"""
    build_recipient_registry(recipient_filepath::String, cpra_filepath::String) -> Vector{Recipient}

Build a registry of `Recipient` objects from recipient and CPRA CSV extracts.

### Arguments
- `recipient_filepath::String`: Path to the recipient CSV file (loaded by `load_recipient`).
- `cpra_filepath::String`: Path to the CPRA CSV file used to compute the latest CPRA per `CAN_ID`.

### Returns
- `Vector{Recipient}`: One `Recipient` per unique `CAN_ID`.

### Notes
- The “most recent” recipient record is defined by the maximum `UPDATE_TM`.
- HLA alleles are taken from the selected most recent record (`CAN_A1/A2`, `CAN_B1/B2`, `CAN_DR1/DR2`).
- The expiration date is inferred from the status history and passed as
  `expiration_date=...` when constructing the `Recipient`.
- Rows with missing required fields are discarded before grouping.
"""
function build_recipient_registry(recipient_filepath::String, cpra_filepath::String)::Vector{Recipient}
    df = load_recipient(recipient_filepath)

    # Keep only rows with required fields
    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]
    dropmissing!(df, required)

    # Compute most recent CPRA per CAN_ID
    cpra_by_CAN_ID = build_last_cpra_registry(cpra_filepath)

    G = groupby(df, :CAN_ID)

    recipients = Vector{Recipient}(undef, length(G))

    for (i , g) in enumerate(G)
        # Infer expiration date from the status history (helper sorts internally)
        exp_date = infer_recipient_expiration_date(g)

        # Sort once here for consistent "most recent" row selection
        sort!(g, :UPDATE_TM, rev=true)

        id = g.CAN_ID[1]

        birth = g.CAN_BTH_DT[1]
        dialysis = g.CAN_DIAL_DT[1]
        arrival = g.CAN_LISTING_DT[1]

        blood = parse_abo(String(g.CAN_BLOOD[1]))

        a1, a2 = g.CAN_A1[1], g.CAN_A2[1]
        b1, b2 = g.CAN_B1[1], g.CAN_B2[1]
        dr1, dr2 = g.CAN_DR1[1], g.CAN_DR2[1]

        cpra = get(cpra_by_CAN_ID, id, 0)

        recipients[i] = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra; expiration_date=exp_date)

    end

    return recipients
end


df = load_recipient(recipient_filepath)

@time recipients = build_recipient_registry(recipient_filepath, cpra_filepath)