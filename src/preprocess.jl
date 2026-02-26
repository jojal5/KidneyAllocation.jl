"""
    parse_hla_int(s::Union{AbstractString,Missing}) -> Union{Int64,Missing}

Parse an HLA-like value by extracting the leading numeric component.

### Examples:
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

### Details
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
    recipient_arrival_departure(df; future_date=Date(2100,1,1)) -> (arrival, departure)

Infer the arrival and departure dates of a single recipient from its status history.
If the recipient is still active at the most recent update, `future_date` is used
as the departure date.
"""
function recipient_arrival_departure(df::AbstractDataFrame, future_date::Date=Date(2100, 1, 1))

    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"
    @assert "CAN_LISTING_DT" in names(df) "Missing column :CAN_LISTING_DT"
    @assert all(==(df.CAN_ID[1]), df.CAN_ID) "All rows must correspond to the same :CAN_ID"

    # Sort the dataframe rows so that the most recent is on top
    idx = sortperm(df.UPDATE_TM; rev=true)
    outcomes = df.OUTCOME[idx]
    updates = df.UPDATE_TM[idx]

    arrival = df.CAN_LISTING_DT[1]

    if outcomes[1] == "1"
        departure = future_date # Arbitrary date after the end of the historic period
    else
        departure = updates[1] # Si transplanté ou retiré
    end

    return arrival, departure
end


"""
    fill_hla_pair!(df, col1, col2)

Replace `missing` values in `col1` (resp. `col2`) by the value in `col2`
(resp. `col1`). Operates in place.
"""
function fill_hla_pair!(df::AbstractDataFrame, col1::Symbol, col2::Symbol)
    df[!, col1] = coalesce.(df[!, col1], df[!, col2])
    df[!, col2] = coalesce.(df[!, col2], df[!, col1])
    return df
end

"""
    fill_hla_pairs!(df, prefix)

Fill missing HLA allele values (A, B, DR) from their paired column in place.
"""
function fill_hla_pairs!(df::AbstractDataFrame, prefix::AbstractString)
    loci = ["A", "B", "DR"]
    for locus in loci
        col1 = Symbol("$(prefix)_$(locus)1")
        col2 = Symbol("$(prefix)_$(locus)2")
        fill_hla_pair!(df, col1, col2)
    end
    return df
end

"""
    load_donor(filepath::String) -> DataFrame

Load a CSV file containing donor information and return a cleaned `DataFrame`.

### Details
Cleaning steps:
- Remove donors for which the decision status (`DECISION`) is missing.
- Keep only donors attributed to List 5 (Administrative to Transplant Québec).
- Drop administrative columns not required for downstream processing.
- Replace `WEIGHT = 0`` and `HEIGHT = 0` by `missing`.
"""
function load_donor(filepath::String)
    df = CSV.read(filepath, DataFrame; missingstring=["-", "", "NULL"])

    fill_hla_pairs!(df, "DON")

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
    donor_from_row(r::DataFrameRow) -> Donor

Construct a `Donor` from a donor `DataFrameRow`. Assume non-missing value.
"""
function donor_from_row(r::DataFrameRow)

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

    donor = Donor(arrival, age, blood, r.DON_A1, r.DON_A2, r.DON_B1, r.DON_B2, r.DON_DR1, r.DON_DR2, kdri)

    return donor
end

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

        donors[i] = donor_from_row(r)

    end

    return donors
end

"""
    load_recipient(filepath::AbstractString)

Load the CSV file containing recipient information and return a cleaned `DataFrame`.

### Details

Cleaning steps:
- Replace the values `24L` and `24Low` with `24` in `CAN_A2`.
- Keep only kidney transplant recipients (exclude multi-organ outcomes).
- Remove recipients with missing dialysis date (`CAN_DIAL_DT`).
- If listing date is missing, replace it by the dialysis date.
- If listing date is after dialysis date, replace the listing date by the dialysis date.
- Keeping only adult recipients (18 years old and older)
- Convert `DateTime` columns to `Date`.
"""
function load_recipient(filepath::AbstractString)
    df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

    # Replacing the value 24L and 24Low with 24 for instance
    df.CAN_A2 = parse_hla_int.(df.CAN_A2)

    fill_hla_pairs!(df, "CAN")

    # Keeping only the recipients for kidney transplant
    filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

    # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_DIAL_DT])

    # If listing time is missing, replacing it by the dialysis time
    df.CAN_LISTING_DT = coalesce.(df.CAN_LISTING_DT, df.CAN_DIAL_DT)

    # If listing is before dialysis, set listing = dialysis
    df.CAN_LISTING_DT = ifelse.(df.CAN_LISTING_DT < df.CAN_DIAL_DT,
        df.CAN_DIAL_DT, df.CAN_LISTING_DT)

    # Keeping only adult recipients
    filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

    df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
    df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

    return df
end

"""
    recipient_from_row(r, cpra=0, expiration_date=nothing) -> Recipient

Construct a `Recipient` from a recipient `DataFrameRow`. Assume non-missing value.
"""
function recipient_from_row(r::DataFrameRow, cpra::Int=0, expiration_date::Union{Nothing,Date}=nothing)

    birth = r.CAN_BTH_DT
    dialysis = r.CAN_DIAL_DT
    arrival = r.CAN_LISTING_DT

    blood = parse_abo(String(r.CAN_BLOOD))

    a1, a2 = r.CAN_A1, r.CAN_A2
    b1, b2 = r.CAN_B1, r.CAN_B2
    dr1, dr2 = r.CAN_DR1, r.CAN_DR2

    recipient = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra; expiration_date=expiration_date)

    return recipient

end

"""
    build_last_cpra_registry(cpra_filepath::String) -> Dict{CAN_ID, Int64}

Build a dictionary mapping each recipient (`CAN_ID`) to their most recent calculated Panel Reactive Antibody (cPRA) value.

### Details
- The CPRA file is read from `cpra_filepath`.
- Rows with missing `CAN_ID` or `CAN_CPRA` are discarded.
- If a patient has many CPRA values, the most recent is retained (by update time according to UPDATE_TM)
- CPRA values are rounded to the nearest integer and returned as `Int64`.

### Returns
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

        # Get last CPRA value
        sort!(g, :UPDATE_TM, rev=true)
        id = g.CAN_ID[1]
        cpra = get(cpra_by_CAN_ID, id, 0)

        recipients[i] = recipient_from_row(first(g), cpra, exp_date)

    end

    return recipients
end