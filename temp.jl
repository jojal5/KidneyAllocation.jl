using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: reconstruct_recipients, build_recipient_registry, load_recipient, build_donor_registry, recipient_arrival_departure, build_last_cpra_registry
import KidneyAllocation: infer_recipient_expiration_date, recipient_from_row

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

recipients = build_recipient_registry(recipient_filepath, cpra_filepath)
donors = build_donor_registry(donor_filepath)


df_recipient = load_recipient(recipient_filepath)

G = groupby(df_recipient, :CAN_ID)

arrival = Vector{Date}(undef, length(G))
departure = Vector{Date}(undef, length(G))
can_id = Vector{Int64}(undef, length(G))

for (i, g) in enumerate(G)
    can_id[i] = g.CAN_ID[1]
    arrival[i], departure[i] = recipient_arrival_departure(g)
end

active_recipient_id = can_id[(arrival.≤Date(2014, 1, 1)).&&(Date(2014, 1, 1).<departure)]


df = filter(row -> row.CAN_ID ∈ active_recipient_id, df_recipient)
dropmissing!(df, :CAN_A1) # Drop recipients with missing HLA, happens when transplanted with a living donor

G = groupby(df, :CAN_ID)

active_recipients = Vector{Recipient}(undef, length(G))

cpra_by_CAN_ID = build_last_cpra_registry(cpra_filepath)

for (i, g) in enumerate(G)

    # Infer expiration date from the status history (helper sorts internally)
    expiration_date = infer_recipient_expiration_date(g)

    # Most_recent CPRA
    id = g.CAN_ID[1]
    cpra = get(cpra_by_CAN_ID, id, 0)

    # Sort once here for consistent "most recent" row selection
    sort!(g, :UPDATE_TM, rev=true)

    active_recipients[i] = recipient_from_row(first(g), cpra, expiration_date)

end

waiting_registry_indices = Vector{Int64}(undef, length(active_recipients))

findfirst(recipients .== active_recipients[i])


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



## Simulation of the allocation process








function infer_recipient_departure(df::AbstractDataFrame, future_date=Date(2100, 1, 1))

    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"
    @assert all(==(df.CAN_ID[1]), df.CAN_ID) "All rows must correspond to the same :CAN_ID"

    # Sort the dataframe rows so that the most recent is on top
    idx = sortperm(df.UPDATE_TM; rev=true)
    outcomes = df.OUTCOME[idx]
    updates = df.UPDATE_TM[idx]

    if outcomes[1] == "1"
        departure = future_date # Arbitrary date after the end of the historic period
    else
        departure = updates[1] # Si transplanté ou retiré
    end

    return departure

end




df_recipient = load_recipient(recipient_filepath)

G = groupby(df_recipient, :CAN_ID)

infer_recipient_departure(G[1000])

dep_by_id = Dict{Int,Date}()

for g in G
    id = g.CAN_ID[1]
    dep = infer_recipient_departure(g)
    dep_by_id[id] = dep
end

df_recipient.DEPARTURE_DT = getindex.(Ref(dep_by_id), df_recipient.CAN_ID)





## Fix load_recipient

import KidneyAllocation: parse_hla_int, fill_hla_pairs!, infer_recipient_expiration_date, load_donor

filepath = recipient_filepath

df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

# Keeping only the recipients for kidney transplant
filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

# Replacing the value 24L and 24Low with 24 for instance
df.CAN_A2 = parse_hla_int.(df.CAN_A2)

fill_hla_pairs!(df, "CAN")

col_to_harmonize = [:CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2, :CAN_LISTING_DT, :CAN_DIAL_DT]

for col in col_to_harmonize
    harmonize_col!(df, col=col)
end

dropmissing!(df, :CAN_DIAL_DT)

enforce_listing_after_dialysis!(df)








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
function infer_recipient_expiration_date(df::AbstractDataFrame)
    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"
    @assert all(==(df.CAN_ID[1]), df.CAN_ID) "All rows must correspond to the same :CAN_ID"

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

# TODO
function infer_recipient_departure(df_recipient, donor_filepath)

    # date de transplantation si transplanté
    # date de retrait si retiré
    # date future si toujours en attente

end



"""
    offered_recipients(df) -> AbstractDataFrame

For a single donor, return the score-ranked recipients that would be offered an
organ. Rows are kept until all accepted offers are reached, up to a maximum of
two acceptances (two kidneys). If fewer acceptances occur, all rows up to the
last acceptance are returned.
"""
function offered_recipients(df::AbstractDataFrame)
    @assert "DON_ID" in names(df) "Missing column :DON_ID"
    @assert "DECISION" in names(df) "Missing column :DECISION"
    @assert "DON_CAN_SCORE" in names(df) "Missing column :DON_CAN_SCORE"

    nrow(df) == 0 && return df
    @assert all(==(df.DON_ID[1]), df.DON_ID) "All rows must correspond to the same :DON_ID"

    df_sort = sort(df, :DON_CAN_SCORE, rev=true)

    accpos = findall(==("Acceptation"), df_sort.DECISION)
    if isempty(accpos)
        return df_sort
    elseif length(accpos) == 1
        return df_sort[1:accpos[1], :]
    else
        return df_sort[1:accpos[2], :]
    end
end

@testset "offered_recipients()" begin

    # Missing columns
    df = DataFrame(DECISION = ["Refus", "Acceptation"], DON_CAN_SCORE = [30, 29])
    @test_throws AssertionError offered_recipients(df)

    df = DataFrame(DON_ID = [1, 1], DON_CAN_SCORE = [30, 29])
    @test_throws AssertionError offered_recipients(df)

    df = DataFrame(DON_ID = [1, 1], DECISION = ["Refus", "Acceptation"])
    @test_throws AssertionError offered_recipients(df)

    # Not a single donor
    df = DataFrame(DON_ID = [1, 2], DECISION = ["Refus", "Acceptation"], DON_CAN_SCORE = [30, 29])
    @test_throws AssertionError offered_recipients(df)

    # Single row accepted
    df = DataFrame(DON_ID = 1, DECISION = "Acceptation", DON_CAN_SCORE = 30)
    df_offered = offered_recipients(df)
    @test df_offered.DON_CAN_SCORE == [30]

    # The first offer is refused, but the second is accepted (not sorted)
    df = DataFrame(DON_ID=[1,1], DECISION=["Acceptation","Refus"], DON_CAN_SCORE=[29,30])
    df_offered = offered_recipients(df)
    @test df_offered.DON_CAN_SCORE == [30,29] 

    # All the offers are refused
    df = DataFrame(DON_ID = [1, 1], DECISION = ["Refus", "Refus"], DON_CAN_SCORE = [30, 29])
    df_offered = offered_recipients(df)
    @test df_offered.DON_CAN_SCORE == [30, 29]

    # The fourth offers is the last accepted.
    df = DataFrame(DON_ID = [1, 1, 1, 1, 1], DECISION = ["Refus", "Acceptation", "Refus", "Acceptation", "Refus"], DON_CAN_SCORE = [30, 29, 28, 27, 26])
    df_offered = offered_recipients(df)
    @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27]

    # Three offers are marked as accepted.
    df = DataFrame(DON_ID = [1, 1, 1, 1, 1], DECISION = ["Refus", "Acceptation", "Refus", "Acceptation", "Acceptation"], DON_CAN_SCORE = [30, 29, 28, 27, 26])
    df_offered = offered_recipients(df)
    @test df_offered.DON_CAN_SCORE == [30, 29, 28, 27]
end





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


df_donors






function number_of_given_kidney(df_donors::AbstractDataFrame)

    kidney_by_don_id = Dict{Int,Int}()

    for g in groupby(df_donors, :DON_ID)
        kidney_by_don_id[g.DON_ID[1]] = count(g.DECISION .== "Acceptation")
    end

    return kidney_by_don_id

end

@time number_of_given_kidney(df_donors)


filter(row->row.DON_ID ==1175, df_donors)





##

function simulate_initial_state_indexed(
    donors::Vector{Donor},
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel;
    start_date::Date=Date(2014, 1, 1),
    nyears::Int=10,
    donor_rate::Real=148.0,
    recipient_rate::Real=272.83,
    origin_date::Date=Date(2000, 1, 1),
    rng::AbstractRNG=Random.default_rng(),
)

    simulation_end = start_date + Year(nyears)

    # Recipients active at start_date: registry indices + their arrival dates
    active_at_start_mask = is_active.(recipients, start_date)
    waiting_registry_indices = findall(active_at_start_mask)
    waiting_arrival_dates = KidneyAllocation.get_arrival.(recipients[waiting_registry_indices])

    # New recipient arrivals: registry indices + simulated arrival dates
    sampled_recipient_indices, sampled_recipient_arrival_dates =
        generate_arrivals(eachindex(recipients), recipient_rate;
            origin=start_date, nyears=nyears, rng=rng)

    append!(waiting_registry_indices, sampled_recipient_indices)
    append!(waiting_arrival_dates, sampled_recipient_arrival_dates)

    # New donor arrivals (used only internally for allocation)
    sampled_donor_indices, sampled_donor_arrival_dates =
        generate_arrivals(eachindex(donors), donor_rate;
            origin=start_date, nyears=nyears, rng=rng)

    # Reconstruct temporary objects for allocation only
    waiting_recipients =
        reconstruct_recipients(recipients, waiting_registry_indices, waiting_arrival_dates)

    arriving_donors =
        reconstruct_donors(donors, sampled_donor_indices, sampled_donor_arrival_dates)

    # Allocate donors (positions refer to waiting_recipients)
    allocated_positions = allocate(arriving_donors, waiting_recipients, dm)

    # Filter non-attributed organ (i.e. ind ==0)
    filter!(>(0), allocated_positions)
    sort!(allocated_positions)

    # Remove transplanted recipients from the indexed representation
    deleteat!(waiting_registry_indices, allocated_positions)
    deleteat!(waiting_arrival_dates, allocated_positions)
    deleteat!(waiting_recipients, allocated_positions)

    # Keep recipients active at end of simulation
    active_at_end_mask = is_active.(waiting_recipients, simulation_end)

    final_recipient_indices = waiting_registry_indices[active_at_end_mask]
    final_arrival_dates = waiting_arrival_dates[active_at_end_mask]

    # Shift arrivals so that simulation_end maps to origin_date
    time_waited = simulation_end .- final_arrival_dates
    shifted_arrival_dates = origin_date .- time_waited

    return final_recipient_indices, shifted_arrival_dates
end




