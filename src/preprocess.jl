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
    load_recipient(filepath::String)

Load the CSV file containing recipient information and return a cleaned DataFrame

##Details

Here are the cleaning steps:
- Replacing the value 24L and 24Low with 24 in CAN_A2
- Keeping only the recipient for kidney transplant, not multi-organs
- Removing the recipient where the dialysis time or the listing time is missing. It happens for recipient that received a kidney from a living donor.
- Converting DateTime format to Date
"""
function load_recipient(filepath::String)

    df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

    # Replacing the value 24L and 24Low with 24 for instance
    df.CAN_A2 = parse_hla_int.(df.CAN_A2)

    # Keeping only the recipients for kidney transplant
    filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

    # Removing the recipient where the dialysis time or the listing time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_LISTING_DT, :CAN_DIAL_DT])

    # Converting DateTime format to Date
    df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
    df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

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

    ("CAN_ID" in names(df_cpra))   || throw(ArgumentError("Missing column :CAN_ID in CPRA file"))
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






# function construct_recipient_database(recipient_filepath::String, cpra_filepath::String)

#     df = load_recipient(recipient_filepath)
#     df_cpra = CSV.read(cpra_filepath, DataFrame)

#     select!(df, Not(:NB_PAST_TRANS))
#     dropmissing!(df)

#     recipients = Recipient[]

#     for id in unique(df.CAN_ID)

#         df_id = filter(row -> row.CAN_ID == id, df)
#         sort!(df_id, :UPDATE_TM, rev=true)

#         df_cpra_id = filter(row -> row.CAN_ID == id, df_cpra)
#         sort!(df_cpra_id, :UPDATE_TM, rev=true)

#         if (df_id.OUTCOME[1] == "TX") || (df_id.OUTCOME[1] == "1")
#             exp_date = nothing # Si transplanté avec un donneur décédé, aucune date d'expiration et on néglige le temps inactif. Même chose si le patient est toujours actif en attente de transplantation.
#         else
#             # Si non transplanté avec un donneur décédé, on prend la dernière date d'attente active pour calculer la date d'expiration
#             ind = findfirst(df_id.OUTCOME .== "1")
#             if ind !== nothing
#                 d = KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1])
#                 exp_date = df_id.UPDATE_TM[end] + Day(d)
#             else # Le patient a été inscrit mais n'a jamais été actif.
#                 exp_date = df_id.UPDATE_TM[end]
#             end

#         end

#         birth = df_id.CAN_BTH_DT[1]
#         dialysis = df_id.CAN_DIAL_DT[1]
#         arrival = df_id.CAN_LISTING_DT[1]
#         blood = KidneyAllocation.parse_abo(df_id.CAN_BLOOD[1])
#         a1 = df_id.CAN_A1[1]
#         a2 = df_id.CAN_A2[1]
#         b1 = df_id.CAN_B1[1]
#         b2 = df_id.CAN_B2[1]
#         dr1 = df_id.CAN_DR1[1]
#         dr2 = df_id.CAN_DR2[1]

#         if isempty(df_cpra_id)
#             cpra = 0 # Si le CPRA est manquant, on prend 0. 
#         else
#             cpra = round(Int64, df_cpra_id.CAN_CPRA[end]) # S'il y a plusieurs CPRA, on prend le plus récent.
#         end

#         r = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra, expiration_date=exp_date)

#         push!(recipients, r)
#     end

#     return recipients

# end