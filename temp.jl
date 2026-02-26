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
    arrival[i], departure[i] = get_arrival_departure(g)
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

@testset "recipient_from_row()" begin
    
    df = DataFrame(CAN_ID = 1, CAN_BTH_DT = Date(1970,1,1), CAN_DIAL_DT = Date(1999,1,1), CAN_LISTING_DT=Date(2000,1,1), CAN_BLOOD = "O",
        CAN_A1 = 3, CAN_A2 = 3, CAN_B1 = 7, CAN_B2 = 8, CAN_DR1 = 7, CAN_DR2 = 8)

    r = first(df)

    r = recipient_from_row(r)

    @test r.birth == Date(1970,1,1)
    @test r.dialysis == Date(1999,1,1)
    @test r.arrival == Date(2000,1,1)
    @test r.blood == O
    @test get_HLA(r) == (3, 3, 7, 8, 7, 8)
    @test r.cpra == 0
    @test isnothing(r.expiration_date)

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


@testset "donor_from_row()" begin
    
    import KidneyAllocation: evaluate_kdri, parse_abo, creatinine_mgdl, get_HLA

    df = DataFrame(DON_ID = 1, DON_DEATH_TM = Date(2000,1,1), DON_AGE = 60, DON_BLOOD="O", HEIGHT = 1.8, WEIGHT = 60., HYPERTENSION=1, DIABETES =0, DEATH = 4, CREATININE =8., DCD = 0,
        DON_A1 = 3, DON_A2 = 3, DON_B1 = 7, DON_B2 = 8, DON_DR1 = 7, DON_DR2 = 8)

    r = first(df)

    d = donor_from_row(r)

    @test d.arrival == Date(2000,1,1)
    @test d.age == 60
    @test d.blood == O
    @test get_HLA(d) == (3, 3, 7, 8, 7, 8)
    @test d.kdri ≈ 3.7152 atol=1e-4

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

function get_arrival_departure(df::AbstractDataFrame, future_date::Date=Date(2024, 1, 1))

    @assert "OUTCOME" in names(df) "Missing column :OUTCOME"
    @assert "UPDATE_TM" in names(df) "Missing column :UPDATE_TM"

    # Sort the dataframe lines so that the most recent is on top
    idx = sortperm(df.UPDATE_TM; rev=true)
    outcomes = uppercase.(String.(df.OUTCOME[idx]))
    updates = df.UPDATE_TM[idx]

    arrival = df.CAN_LISTING_DT[1]

    if outcomes[1] == "1"
        departure = future_date # Arbitrary date after the end of the historic period
    else
        departure = updates[1] # Si transplanté ou retiré
    end

    return arrival, departure
end


