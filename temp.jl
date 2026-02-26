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

active_recipient_id = can_id[arrival.<Date(2014, 1, 1).<departure]


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

Return the `CAN_ID`s of recipients active on `date`, based on the
recipient status history in the CSV file.
"""
function active_recipient_ids(recipient_filepath::String, date::Date)
    df = load_recipient(recipient_filepath)

    out = Int64[]
    for g in groupby(df, :CAN_ID)
        a, d = recipient_arrival_departure(g)
        if a ≤ date < d
            push!(out, g.CAN_ID[1])
        end
    end

    return out
end


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


recipients_from_can_ids(recipient_filepath,cpra_filepath,[1,3])



# TODO : gérer les données manquantes (peut-être dans active_recipient_ids)
function get_active_recipients(recipient_filepath::String, cpra_filepath::String, date::Date)

    can_ids = active_recipient_ids(recipient_filepath::String, date::Date)

    recipients = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)

    return recipients
end

can_ids = active_recipient_ids(recipient_filepath, Date(2014,1,1,))
r = recipients_from_can_ids(recipient_filepath, cpra_filepath, can_ids)


get_active_recipients(recipient_filepath, cpra_filepath, Date(2014,1,1))













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




