using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: reconstruct_recipients, build_recipient_registry, load_recipient, build_donor_registry


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

recipients = build_recipient_registry(recipient_filepath, cpra_filepath)
donors = build_donor_registry(donor_filepath)

@load "src/SyntheticData/TreeDecisionModel.jld2"



import KidneyAllocation: is_active, build_last_cpra_registry, infer_recipient_expiration_date, parse_abo

count(is_active.(recipients, Date(2014, 1, 1)))

df_recipient = load_recipient(recipient_filepath)

G = groupby(df_recipient, :CAN_ID)

arrival = Vector{Date}(undef, length(G))
departure = Vector{Date}(undef, length(G))
can_id = Vector{Int64}(undef, length(G))

for (i, g) in enumerate(G)
    can_id[i] = g.CAN_ID[1]
    arrival[i], departure[i] = candidate_arrival_departure(g)
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

for i in eachindex(active_recipients)
    waiting_registry_indices[i] = findfirst(recipients .== active_recipients[i])
end










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



