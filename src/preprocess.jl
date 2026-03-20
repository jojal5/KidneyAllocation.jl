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

    col_to_harmonize = [:DON_AGE, :DON_BLOOD, :DON_DEATH_TM, :WEIGHT, :HEIGHT, :ETHNICITY, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :HCV, :DCD]

    for col in col_to_harmonize
        harmonize_col!(df, col=col, idcol=:DON_ID)
    end

    dropmissing!(df, [:DECISION])

    # Only attribution type 5
    filter!(row -> row.ATT_TYPE == 5, df)

    cols_to_drop = [:DON_LAB_DT_TM, :ATT_TYPE]
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
    build_donor_registry(filepath::String) -> Dict{Int,Donor}

Build a dictionary of `Donor` objects keyed by `DON_ID` from a donor CSV file.

Rows with missing values required to construct a donor and compute KDRI are
discarded. For each `DON_ID`, one `Donor` is constructed from the corresponding
row using `donor_from_row`.

### Notes
- The donor arrival date is derived from the donor death date (`DON_DEATH_TM`).
- Rows with incomplete data required for KDRI computation are discarded.
"""
function build_donor_registry(filepath::String)::Dict{Int,Donor}
    df = load_donor(filepath)

    kdri_required = Symbol[
        :DON_ID, :DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD
    ]
    donor_required = Symbol[
        :DON_DEATH_TM, :DON_AGE, :DON_BLOOD, :DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2
    ]

    dropmissing!(df, kdri_required)
    dropmissing!(df, donor_required)

    donor_by_don_id = Dict{Int,Donor}()

    for g in groupby(df, :DON_ID)
        r = first(g)
        don_id = r.DON_ID
        donor_by_don_id[don_id] = donor_from_row(r)
    end

    return donor_by_don_id
end

"""
    kidneys_given_by_donor(df_donors) -> Dict{Int,Int}

Return the number of given kidneys for each `DON_ID`.
"""
function kidneys_given_by_donor(df_donors::AbstractDataFrame)
    @assert "DON_ID" in names(df_donors) "Missing column :DON_ID"
    @assert "STATUS" in names(df_donors) "Missing column :STATUS"

    kidney_by_don_id = Dict{Int,Int}()

    for g in groupby(df_donors, :DON_ID)
       kidney_by_don_id[g.DON_ID[1]] = count(==("TX"), skipmissing(g.STATUS))
    end

    return kidney_by_don_id

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

    # Harmonize columns when some rows have missig value for a same recipient
    col_to_harmonize = [:CAN_BTH_DT, :CAN_BLOOD, :CAN_DIAL_DT]
    for col in col_to_harmonize
        harmonize_col!(df, col=col, idcol=:CAN_ID )
    end

    # Keeping only the recipients for kidney transplant
    filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

    # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_DIAL_DT])

    # If listing time is missing, replacing it by the dialysis time. If listing is before dialysis, set listing = dialysis
    enforce_listing_after_dialysis!(df)

    # Transform DateTime in Date
    df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
    df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

    # Keeping only adult recipients
    filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

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
    build_recipient_registry(recipient_filepath, cpra_filepath) -> Dict{Int,Recipient}

Build a dictionary of `Recipient` objects keyed by `CAN_ID`.

### Arguments
- `recipient_filepath::String`: Path to the recipient CSV file (loaded by `load_recipient`).
- `cpra_filepath::String`: Path to the CPRA CSV file used to compute the latest CPRA per `CAN_ID`.

### Notes
- The expiration date is inferred from the status history and passed as
  `expiration_date=...` when constructing the `Recipient`.
- Rows with missing required fields are discarded.
"""
function build_recipient_registry(recipient_filepath::String, cpra_filepath::String)::Dict{Int,Recipient}
    df = load_recipient(recipient_filepath)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]
    dropmissing!(df, required)

    cpra_by_can_id = build_last_cpra_registry(cpra_filepath)
    recipient_by_can_id = Dict{Int,Recipient}()

    for g in groupby(df, :CAN_ID)
        exp_date = infer_recipient_expiration_date(g)

        r = first(g)
        can_id = r.CAN_ID

        cpra = get(cpra_by_can_id, can_id, 0)
        recipient_by_can_id[can_id] = recipient_from_row(r, cpra, exp_date)
    end

    return recipient_by_can_id
end

"""
    enforce_listing_after_dialysis!(df) -> AbstractDataFrame

Modify `df` in place to ensure `CAN_LISTING_DT ≥ CAN_DIAL_DT`. If
`CAN_LISTING_DT` is missing or earlier than `CAN_DIAL_DT`, it is replaced by
`CAN_DIAL_DT`.

# Notes
The DataFrame is expected to have the structure returned by
[`load_recipient`](@ref).
"""
function enforce_listing_after_dialysis!(df::AbstractDataFrame)
    @assert "CAN_LISTING_DT" in names(df) "Missing column :CAN_LISTING_DT"
    @assert "CAN_DIAL_DT" in names(df) "Missing column :CAN_DIAL_DT"

    for r in eachrow(df)
        dial = r.CAN_DIAL_DT
        list = r.CAN_LISTING_DT

        # If dialysis date is missing, we cannot enforce the constraint
        if ismissing(dial)
            continue
        end

        if ismissing(list) || list < dial
            r.CAN_LISTING_DT = dial
        end
    end

    return df
end

"""
    harmonize_col!(df; col, idcol=:CAN_ID) -> DataFrame

Ensure that `col` has at most one non-missing value per `idcol` group and
fill missing entries with that value. Throws an error if multiple distinct
non-missing values are found within a group.
"""
function harmonize_col!(df::DataFrame; col::Symbol, idcol::Symbol=:CAN_ID)

    val_by_id = Dict{Int,Union{Missing,eltype(skipmissing(df[!, col]))}}()

    for g in groupby(df, idcol)
        id = g[1, idcol]
        vals = unique(skipmissing(g[!, col]))
        if length(vals) == 0
            val_by_id[id] = missing
        elseif length(vals) == 1
            val_by_id[id] = first(vals)
        else
            throw(ArgumentError("Multiple $(col) values for $(idcol)=$id: $(collect(vals))"))
        end
    end

    df[!, col] = coalesce.(df[!, col], getindex.(Ref(val_by_id), df[!, idcol]))
    return df
end