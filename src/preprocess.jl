"""
    parse_hla_int(s::Union{AbstractString,Missing}) -> Union{Int64,Missing}

Parse an HLA-like value by extracting the leading numeric component.

## Examples:
- `"24"`     → 24
- `"24L"`    → 24
- `"24Low"`  → 24
- `missing`  → `missing`

Throws an `ArgumentError` if a non-missing value does not start with digits.
"""
function parse_hla_int(s::Union{AbstractString,Missing})::Union{Int64,Missing}
    s === missing && return missing

    m = match(r"^\s*(\d+)", s)
    m === nothing && throw(ArgumentError("Invalid HLA value: \"$s\""))

    return parse(Int64, m.captures[1])
end

"""
    infer_recipient_expiration_date(df::AbstractDataFrame) -> Union{Date,Nothing}

Infer the expiration date for a single recipient from their longitudinal status history.

## Details
The input `df` must correspond to a single recipient (`CAN_ID`) and contain
at least the columns `:OUTCOME` and `:UPDATE_TM`.

The history is sorted internally by `:UPDATE_TM` in descending order
(most recent first).

Rules (case-insensitive):
- If the most recent outcome is `"TX"` or `"1"`, return `nothing`.
- Otherwise:
  - If `"1"` never occurs in the history, return the oldest `UPDATE_TM`.
  - If `"1"` occurs but is not the most recent status, return the `UPDATE_TM`
    immediately preceding the first `"1"` in the descending timeline.
"""
function infer_recipient_expiration_date(df::AbstractDataFrame)::Union{Date,Nothing}
    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"

    # Sort the dataframe lines so that the most recent is on top
    idx = sortperm(df.UPDATE_TM; rev=true)
    outcomes = uppercase.(String.(df.OUTCOME[idx]))
    updates = df.UPDATE_TM[idx]

    if outcomes[1] == "TX" || outcomes[1] == "1" # If the recipient has been transplanted or still active, then there is no expiration date.
        return nothing
    else
        ind_active = findfirst(==("1"), outcomes)
        if ind_active === nothing # Recipient was never active
            return updates[end]   # oldest date
        else
            return updates[ind_active-1] # # Si non transplanté avec un donneur décédé, on prend la dernière date d'attente active pour calculer la date d'expiration
        end
    end
end


"""
    load_donor(filepath::String) -> DataFrame

Load a CSV file containing donor information and return a cleaned `DataFrame`.

## Details
Cleaning steps:
- Remove donors for which the decision status (`DECISION`) is missing.
- Keep only donors attributed to List 5 (Administrative to Transplant Québec).
- Drop administrative columns not required for downstream processing.
- Replace `WEIGHT = 0`` and `HEIGHT = 0` by `missing`.
"""
function load_donor(filepath::String)
    df = CSV.read(filepath, DataFrame; missingstring=["-", "", "NULL"])

    dropmissing!(df, [:DECISION])

    # Only attribution type 5
    filter!(row -> row.ATT_TYPE == 5, df)

    cols_to_drop = [:DON_LAB_DT_TM, :ATT_TYPE, :STATUS]
    select!(df, Not(cols_to_drop))

    allowmissing!(df, [:WEIGHT, :HEIGHT])
    df.WEIGHT[df.WEIGHT.≈0.] .= missing
    df.HEIGHT[df.HEIGHT.≈0.] .= missing

    return df
end

"""
    build_donor_registry(filepath::String) -> Vector{Donor}

Construct a registry of `Donor` objects from a CSV file containing donor information.

## Details
This function loads donor data using [`load_donor`](@ref) and applies the following steps:
- Removes donors with missing information required for donor characterization and KDRI computation.
- Computes the Kidney Donor Risk Index (KDRI) for each donor.
- Constructs a `Donor` object for each valid donor record.

## Notes
- The donor arrival date is derived from the donor death date (`DON_DEATH_TM`).
- Binary clinical indicators (hypertension, diabetes, DCD) are inferred from coded values.
- Rows with incomplete data required for KDRI computation are discarded.
- Keep only accepted offers (i.e., `DECISION == "Acceptation"`).

## Returns
A vector of `Donor` objects representing the donor registry.
"""
function build_donor_registry(filepath::String)

    df = load_donor(filepath)

    # Only consideing accepted offers to avois duplication
    filter!(row -> row.DECISION == "Acceptation", df)

    # Si on a une taille et un poids manquant, on remplace par les valeurs moyennes
    # df.WEIGHT[df.WEIGHT.≈0.] .= 80.
    # df.HEIGHT[df.HEIGHT.≈0.] .= 170.

    dropmissing!(df)

    donors = Donor[]

    for r in eachrow(df)

        age = r.DON_AGE
        height = r.HEIGHT
        weight = r.WEIGHT
        hypertension = r.HYPERTENSION == 1
        diabetes = r.DIABETES == 1
        cva = r.DEATH ∈ [4, 16]
        creatinine = creatinine_mgdl(r.CREATININE)
        dcd = r.DCD == 1 # TODO À VÉRIFIER si c'est bien 1, sinon c'est 2 (Anastasiya a confirmé le code)

        kdri = evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

        arrival = Date(r.DON_DEATH_TM)
        blood = KidneyAllocation.parse_abo(r.DON_BLOOD)

        a1 = r.DON_A1
        a2 = r.DON_A2
        b1 = r.DON_B1
        b2 = r.DON_B2
        dr1 = r.DON_DR1
        dr2 = r.DON_DR2

        d = Donor(arrival, age, blood, a1, a2, b1, b2, dr1, dr2, kdri)

        push!(donors, d)

    end

    return donors

end

"""
    load_recipient(filepath::AbstractString)

Load the CSV file containing recipient information and return a cleaned `DataFrame`.

## Details

Cleaning steps:
- Replace the values `24L` and `24Low` with `24` in `CAN_A2`.
- Keep only kidney transplant recipients (exclude multi-organ outcomes).
- Remove recipients with missing dialysis date (`CAN_DIAL_DT`).
- If listing date is missing, replace it by the dialysis date.
- If listing date is after dialysis date, replace the listing date by the dialysis date.
- Convert `DateTime` columns to `Date`.
"""
function load_recipient(filepath::AbstractString)
    df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

    # Replacing the value 24L and 24Low with 24 for instance
    df.CAN_A2 = parse_hla_int.(df.CAN_A2)

    # Keeping only the recipients for kidney transplant
    filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

    # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_DIAL_DT])

    # If listing time is missing, replacing it by the dialysis time
    df.CAN_LISTING_DT = coalesce.(df.CAN_LISTING_DT, df.CAN_DIAL_DT)

    # If listing is before dialysis, set listing = dialysis
    df.CAN_LISTING_DT = ifelse.(df.CAN_LISTING_DT < df.CAN_DIAL_DT,
                                df.CAN_DIAL_DT, df.CAN_LISTING_DT)

    df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    df.CAN_DIAL_DT    = Date.(df.CAN_DIAL_DT)
    df.UPDATE_TM      = passmissing(Date).(df.UPDATE_TM)

    return df
end

"""
    build_last_cpra_registry(cpra_filepath::String) -> Dict{CAN_ID, Int64}

Build a dictionary mapping each recipient (`CAN_ID`) to their most recent calculated Panel Reactive Antibody (cPRA) value.

## Details
- The CPRA file is read from `cpra_filepath`.
- Rows with missing `CAN_ID` or `CAN_CPRA` are discarded.
- If a patient has many CPRA values, the most recent is retained (by update time according to UPDATE_TM)
- CPRA values are rounded to the nearest integer and returned as `Int64`.

## Returns
A dictionary mapping `CAN_ID` to the most recent CPRA value.
"""
function build_last_cpra_registry(cpra_filepath::String)

    df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])

    ("CAN_ID" in names(df_cpra)) || throw(ArgumentError("Missing column :CAN_ID in CPRA file"))
    ("UPDATE_TM" in names(df_cpra)) || throw(ArgumentError("Missing column :UPDATE_TM in CPRA file"))
    ("CAN_CPRA" in names(df_cpra)) || throw(ArgumentError("Missing column :CAN_CPRA in CPRA file"))

    # Drop rows missing the essentials
    dropmissing!(df_cpra, [:CAN_ID, :CAN_CPRA, :UPDATE_TM])

    cpra_by_CAN_ID = Dict{eltype(df_cpra.CAN_ID),Int64}()

    for g in groupby(df_cpra, :CAN_ID)
        sort!(g, :UPDATE_TM, rev=true)  # most recent first
        cpra_by_CAN_ID[g.CAN_ID[1]] = Int64(round(g.CAN_CPRA[1]))
    end

    return cpra_by_CAN_ID
end

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

    recipients = Recipient[]

    for g in groupby(df, :CAN_ID)
        # Infer expiration date from the status history (helper sorts internally)
        exp_date = infer_recipient_expiration_date(g)

        # Sort once here for consistent "most recent" row selection
        sort!(g, :UPDATE_TM, rev=true)

        id = g.CAN_ID[1]

        birth = g.CAN_BTH_DT[1]
        dialysis = g.CAN_DIAL_DT[1]
        arrival = g.CAN_LISTING_DT[1]

        blood = KidneyAllocation.parse_abo(String(g.CAN_BLOOD[1]))

        a1, a2 = g.CAN_A1[1], g.CAN_A2[1]
        b1, b2 = g.CAN_B1[1], g.CAN_B2[1]
        dr1, dr2 = g.CAN_DR1[1], g.CAN_DR2[1]

        cpra = get(cpra_by_CAN_ID, id, 0)

        push!(recipients,
            Recipient(birth, dialysis, arrival, blood,
                a1, a2, b1, b2, dr1, dr2, cpra;
                expiration_date=exp_date))
    end

    return recipients
end