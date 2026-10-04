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
    check_df_columns(df::AbstractDataFrame, cols::Symbol...) -> Nothing

Verify that `df` contains every column in `cols` and that none contains `missing`. Throw an `ArgumentError` otherwise.
"""
function check_df_columns(df::AbstractDataFrame, cols::Symbol...)
    available_cols = propertynames(df)

    for col in cols
        col ∈ available_cols ||
            throw(ArgumentError("Missing column :$col"))

        any(ismissing, df[!, col]) &&
            throw(ArgumentError("Column :$col contains missing values"))
    end

    return nothing
end

"""
    check_df_column_constant(df::AbstractDataFrame, col::Symbol) -> Nothing

Verify that all values in column `col` are identical. Throw an `ArgumentError` otherwise.
"""
function check_df_column_constant(df::AbstractDataFrame, cols::Symbol...)

    isempty(df) && return nothing

    for col in cols
        check_df_columns(df, col)
    
        values = df[!, col]
        all(==(first(values)), values) ||
            throw(ArgumentError(
                "All rows must have the same value in column :$col",
            ))
    end

    return nothing
end

"""
    check_df_at_most_one_exit_outcome(df, outcomes=...) -> Nothing

Verify that a single recipient history contains at most one of the specified
exit outcomes. Outcome codes are matched case-insensitively.
"""
function check_df_at_most_one_exit_outcome(
    df::AbstractDataFrame,
    outcomes::AbstractVector{<:AbstractString}=[
        "X", "TX VIVANT", "DCD", "TX",
    ],
)
    check_df_columns(df, :OUTCOME)
    check_df_column_constant(df, :CAN_ID)

    exit_outcomes = Set(
        uppercase(strip(string(outcome))) for outcome in outcomes
    )

    n_exit_outcomes = count(
        row -> uppercase(strip(string(row.OUTCOME))) ∈ exit_outcomes,
        eachrow(df),
    )

    n_exit_outcomes ≤ 1 ||
        throw(ArgumentError(
            "More than one exit outcome found for this recipient",
        ))

    return nothing
end

"""
    get_last_update(df::AbstractDataFrame) -> Tuple{DateTime, String}

Return the time and outcome of the most recent update retained by
`filter_outcomes` for a single recipient. 
"""
function get_last_update(df::AbstractDataFrame)::Tuple{DateTime, String}

    isempty(df) && throw(ArgumentError("DataFrame is empty"))

    check_df_columns(df, :OUTCOME, :UPDATE_TM)
    check_df_column_constant(df, :CAN_ID)

    filtered_df = filter_outcomes(df)
    last_index = argmax(filtered_df.UPDATE_TM)

    return (filtered_df.UPDATE_TM[last_index], String(filtered_df.OUTCOME[last_index]))
end

"""
    get_expiration_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

Return the time of the latest retained update if its outcome indicates a
permanent exit from the waiting list; otherwise return `nothing`.
"""
function get_expiration_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
    expiration_outcomes = ("X", "TX VIVANT", "DCD")

    date, outcome = get_last_update(df)

    return uppercase(strip(outcome)) ∈ expiration_outcomes ? date : nothing
end

"""
    get_exit_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

Return the time of the first recorded exit from the waiting list for a single
recipient, or `nothing` if no exit outcome is recorded.
"""
function get_exit_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
    exit_outcomes = ("TX", "TX VIVANT", "X", "DCD")

    date, outcome = get_last_update(df)

    return uppercase(strip(outcome)) ∈ exit_outcomes ? date : nothing
end

"""
    first_exit_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

Return the earliest update time associated with an exit outcome for a single
recipient, or `nothing` if no exit outcome is recorded.
"""
function first_exit_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
    isempty(df) && return nothing

    check_df_column_constant(df, :CAN_ID)
    check_df_columns(df, :OUTCOME, :UPDATE_TM)

    exit_outcomes = ("X", "TX VIVANT", "DCD", "TX")
    earliest_exit_time = nothing

    for row in eachrow(df)
        outcome = uppercase(strip(string(row.OUTCOME)))

        if outcome ∈ exit_outcomes
            update_time = DateTime(row.UPDATE_TM)

            earliest_exit_time = isnothing(earliest_exit_time) ?
                update_time :
                min(earliest_exit_time, update_time)
        end
    end

    return earliest_exit_time
end

"""
    filter_outcomes(df::AbstractDataFrame) -> DataFrame

Return records from the first waiting-list episode of a single recipient.

Rows after the first exit outcome are removed, while the exit record is
retained. Exact duplicate retained rows are removed.
(Happen for recipients 501, 827 and 2051)

If every retained `UPDATE_TM` precedes the recorded `CAN_LISTING_DT`,
`CAN_LISTING_DT` is replaced by `CAN_DIAL_DT` as a proxy for the unavailable
original listing date.

`:UPDATE_TM` must contain `DateTime` values so that outcomes occurring on the
same calendar date can be ordered correctly.
"""
function filter_outcomes(df::AbstractDataFrame)::DataFrame
    isempty(df) && return DataFrame(df)

    check_df_column_constant(df, :CAN_ID, :CAN_LISTING_DT, :CAN_DIAL_DT)
    check_df_columns(df, :OUTCOME, :UPDATE_TM)

    update_type = eltype(skipmissing(df[!, :UPDATE_TM]))

    update_type <: DateTime || throw(ArgumentError("Column :UPDATE_TM must contain DateTime values; got $update_type"))

    exit_time = first_exit_date(df)

    if exit_time === nothing
        filtered_df = DataFrame(df)
    else
        filtered_df = filter(:UPDATE_TM => date -> date ≤ exit_time, df)

        listing_date = first(filtered_df.CAN_LISTING_DT)
        dialysis_date = convert(typeof(listing_date), first(filtered_df.CAN_DIAL_DT))

        if maximum(filtered_df.UPDATE_TM) < listing_date
            filtered_df.CAN_LISTING_DT .= dialysis_date
        end

    end

    # Removing duplicated rows (happens for recipient 501, 827 and 2051)
    unique!(filtered_df)

    return filtered_df
end

"""
    recipient_active_waiting_proportion(df::AbstractDataFrame) -> Float64

Return the proportion in `[0, 1]` of the observed waiting period during which
a single recipient is active. The observed period runs from `CAN_LISTING_DT` to the
most recent `UPDATE_TM`.
"""
function recipient_active_waiting_proportion(df::AbstractDataFrame)::Float64
    isempty(df) && return 0.0

    check_df_columns(df, :OUTCOME, :UPDATE_TM)
    check_df_column_constant(df, :CAN_ID, :CAN_LISTING_DT, :CAN_DIAL_DT)
    check_df_at_most_one_exit_outcome(df)

    recipient_df = sort(df, :UPDATE_TM)

    listing_time = DateTime(first(recipient_df.CAN_LISTING_DT))
    last_update = DateTime(last(recipient_df.UPDATE_TM))
    total_waiting = last_update - listing_time

    total_waiting ≤ Millisecond(0) && return 0.0

    active_waiting = Millisecond(0)
    previous_time = listing_time
    is_active = true

    for row in eachrow(recipient_df)
        update_time = DateTime(row.UPDATE_TM)

        if update_time > previous_time
            is_active && (active_waiting += update_time - previous_time)
            previous_time = update_time
        end

        is_active = strip(string(row.OUTCOME)) == "1"
    end

    return Dates.value(active_waiting) / Dates.value(total_waiting)
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

    # # Manually removing duplicated rows for TX (TODO check in new dataset version if these errors are still present)
    # if (df.OUTCOME[3104] == df.OUTCOME[3105]) & (df.UPDATE_TM[3104] == df.UPDATE_TM[3105])
    #     deleteat!(df, 3104)
    # end

    # Replacing the value 24L and 24Low with 24 for instance
    df.CAN_A2 = parse_hla_int.(df.CAN_A2)

    fill_hla_pairs!(df, "CAN")

    # Harmonize columns when some rows have missig value for a same recipient
    col_to_harmonize = [:CAN_BTH_DT, :CAN_BLOOD, :CAN_DIAL_DT]
    for col in col_to_harmonize
        harmonize_col!(df, col=col, idcol=:CAN_ID )
    end

    # Keeping only the recipients for kidney transplant
    filter!(row -> uppercase(row.OUTCOME) ∈ ("TX", "1", "0", "X", "DCD", "TX VIVANT"), df)

    # # Keeping only the first outcome for each recipient. Some recipients are registered more than once and key information are lost for the additional registration
    # filter_outcomes

    # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_DIAL_DT])

    # # If listing time is missing, replacing it by the dialysis time. If listing is before dialysis, set listing = dialysis
    # enforce_listing_after_dialysis!(df)

    # If listing time is missing, replacing it by the dialysis time.
    coalesce_listing!(df)

    # Transform DateTime in Date
    # df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    # df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
    # df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

    # Keeping only adult recipients
    filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

    return df
end

"""
    recipient_from_row(
        r::DataFrameRow,
        cpra::Integer=0;
        expiration_date=nothing,
        active_waiting_proportion=1.0,
    ) -> Recipient

Construct a `Recipient` from one row of a recipient's history. Required values
in `r` are assumed non-missing.
"""
function recipient_from_row(
    r::DataFrameRow,
    cpra::Integer=0;
    expiration_date::Union{Date,DateTime,Nothing}=nothing,
    active_waiting_proportion::Real=1.0,
)
    birth = r.CAN_BTH_DT
    dialysis = r.CAN_DIAL_DT
    arrival = r.CAN_LISTING_DT

    blood = parse_abo(String(r.CAN_BLOOD))

    return Recipient(
        birth, dialysis, arrival, blood,
        r.CAN_A1, r.CAN_A2,
        r.CAN_B1, r.CAN_B2,
        r.CAN_DR1, r.CAN_DR2,
        cpra;
        expiration_date=expiration_date,
        active_waiting_proportion=active_waiting_proportion,
    )
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

        filtered_df = filter_outcomes(g)

        # exp_date = infer_recipient_expiration_date(filtered_df)
        exp_date = get_expiration_date(filtered_df)
        active_waiting_proportion = recipient_active_waiting_proportion(filtered_df)

        r = first(filtered_df)
        can_id = r.CAN_ID

        cpra = get(cpra_by_can_id, can_id, 0)
        recipient_by_can_id[can_id] = recipient_from_row(r, cpra; expiration_date = exp_date, active_waiting_proportion = active_waiting_proportion)
    end

    return recipient_by_can_id
end

"""
    coalesce_listing!(df::AbstractDataFrame) -> AbstractDataFrame

Replace each missing value in `:CAN_LISTING_DT` with the corresponding value
in `:CAN_DIAL_DT`, modifying `df` in place.

`df` must contain `:CAN_LISTING_DT` and `:CAN_DIAL_DT`. Values in
`:CAN_DIAL_DT` must be non-missing.
"""
function coalesce_listing!(df::AbstractDataFrame)
    for col in (:CAN_LISTING_DT, :CAN_DIAL_DT)
        col ∈ propertynames(df) ||
            throw(ArgumentError("Missing column :$col"))
    end

    any(ismissing, df[!, :CAN_DIAL_DT]) &&
        throw(ArgumentError("Column :CAN_DIAL_DT contains missing values"))

    for row in eachrow(df)
        if ismissing(row.CAN_LISTING_DT)
            row.CAN_LISTING_DT = row.CAN_DIAL_DT
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

"""
    retrieve_observed_waiting_list(recipient_filepath, date) -> Vector{Int64}

Return the IDs of recipients registered on the waiting list at `date`, based
on the recipient history stored in `recipient_filepath`.
"""
function retrieve_observed_waiting_list(recipient_filepath::String, date::DateTime)::Vector{Int64}

    df = load_recipient(recipient_filepath)
    registered_ids = Int64[]

    for group in groupby(df, :CAN_ID)
        filtered_group = filter_outcomes(group)

        arrival_time = DateTime(first(filtered_group.CAN_LISTING_DT))
        exit_time = get_exit_date(filtered_group)

        if arrival_time ≤ date &&
           (isnothing(exit_time) || date ≤ exit_time)
            push!(registered_ids, first(filtered_group.CAN_ID))
        end
    end

    return registered_ids
end

retrieve_observed_waiting_list(recipient_filepath::String, date::Date) = retrieve_observed_waiting_list(recipient_filepath, DateTime(date))
