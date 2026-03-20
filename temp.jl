using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation


"""
    indices_in_registry(recipients, registry) -> Vector{Int}

Return the index of each recipient in `registry`, or `0` if not found.
"""
function indices_in_registry(
    recipients::AbstractVector{Recipient},
    registry::AbstractVector{Recipient},
)
    idx = Vector{Int}(undef, length(recipients))
    for (i, r) in enumerate(recipients)
        j = findfirst(==(r), registry)
        idx[i] = j === nothing ? 0 : j
    end
    return idx
end

@testset "indices_in_registry()" begin
    import KidneyAllocation.indices_in_registry
    # Registry (tiny and fakes)
    registry = [
        Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), A, 69, 2403, 7, 35, 4, 103, 10),
        Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),
    ]

    # Two first recipients in registry, but not the third
    recipients = vcat(registry[2], registry[1], Recipient(Date(1965, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),)

    @test indices_in_registry(recipients, registry) == [2, 1, 0]
end


"""
    active_recipient_ids(recipient_filepath, date) -> Vector{Int}

Return the `CAN_ID`s of recipients active on `date` after removing rows with
missing values in required fields.
"""
function active_recipient_ids(recipient_filepath::String, date::Date)
    df = load_recipient(recipient_filepath)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]
    dropmissing!(df, required)

    out = Int[]
    for g in groupby(df, :CAN_ID)
        a, d = recipient_arrival_departure(g)
        if a ≤ date < d
            push!(out, g.CAN_ID[1])
        end
    end

    return out
end

# Registry (tiny and fakes)
registry = [
    Recipient(Date(1979, 1, 1), Date(1995, 1, 1), Date(1998, 1, 1), O, 68, 203, 39, 77, 15, 17, 0),
    Recipient(Date(1981, 1, 1), Date(1997, 1, 1), Date(2000, 6, 1), A, 69, 2403, 7, 35, 4, 103, 10),
    Recipient(Date(1963, 1, 1), Date(1998, 1, 1), Date(2001, 5, 1), B, 25, 68, 67, 5102, 11, 16, 20),
]


"""
    recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids) -> Vector{Recipient}

Load recipient history and return one `Recipient` per `CAN_ID` in `can_ids`.

Throws an error if a requested `CAN_ID` is not found.
"""
function recipients_from_can_ids(
    recipient_filepath::String,
    cpra_filepath::String,
    can_ids::AbstractVector{<:Int},
)::Vector{Recipient}

    df_recipient = load_recipient(recipient_filepath)
    cpra_by_can_id = build_last_cpra_registry(cpra_filepath)

    # Group once
    G = groupby(df_recipient, :CAN_ID)

    # Build mapping CAN_ID -> SubDataFrame
    group_by_id = Dict(first(g.CAN_ID) => g for g in G)

    recipients = Vector{Recipient}(undef, length(can_ids))

    for (i, id) in enumerate(can_ids)

        haskey(group_by_id, id) ||
            throw(ArgumentError("CAN_ID $id not found in recipient file"))

        g = group_by_id[id]

        expiration_date = infer_recipient_expiration_date(g)

        sort!(g, :UPDATE_TM, rev=true)

        cpra = get(cpra_by_can_id, id, 0)

        recipients[i] = recipient_from_row(first(g), cpra, expiration_date)
    end

    return recipients
end


recipients_from_can_ids(recipient_filepath, cpra_filepath, [1, 3])


function get_active_recipients(recipient_filepath::String, cpra_filepath::String, date::Date)

    can_ids = active_recipient_ids(recipient_filepath::String, date::Date)

    recipients = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)

    return recipients
end

can_ids = active_recipient_ids(recipient_filepath, Date(2014, 1, 1,))
r = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)

r = get_active_recipients(recipient_filepath, cpra_filepath, Date(2014, 1, 1))


"""
    transplant_dates_by_recipient(df_donors) -> Dict{Int,Date}

Return a dictionary mapping each transplanted `CAN_ID` to its transplant date
(`DON_DEATH_TM`), based on accepted donor–recipient pairs.
"""
function transplant_dates_by_recipient(df_donors::AbstractDataFrame)

    df = unique(df_donors, [:CAN_ID, :DON_ID])
    filter!(row -> row.DECISION == "Acceptation", df)

    transplant_date_by_id = Dict{Int,Date}()

    for r in eachrow(df)
        transplant_date_by_id[r.CAN_ID] = r.DON_DEATH_TM
    end

    return transplant_date_by_id
end

df_donors = load_donor(donor_filepath)
transplant_dates_by_recipient(df_donors)






