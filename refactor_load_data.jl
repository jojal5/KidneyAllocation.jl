using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Random, Statistics, Test

using Gadfly

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

import KidneyAllocation: parse_hla_int, fill_hla_pairs!, harmonize_col!, coalesce_listing!, years_between, filter_dialysis_episode, filter_offers_until_last_acceptance
import KidneyAllocation: creatinine_mgdl, evaluate_kdri, parse_abo

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
        harmonize_col!(df, col=col, idcol=:CAN_ID)
    end

    # Keeping only the recipients for kidney transplant
    filter!(row -> uppercase(strip(row.OUTCOME)) ∈ ("TX", "1", "0", "X", "DCD", "TX VIVANT"), df)

    # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_DIAL_DT])

    # If listing time is missing, replacing it by the dialysis time.
    coalesce_listing!(df)

    # Keeping only adult recipients
    filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

    # Keeping only the first outcome after the registered dialysis date for each recipient. 
    # Some recipients are registered more than once and key information are lost for the additional registrations.
    # We keep the history corresponding to the recorded dialysis date.
    filtered_df = similar(df, 0)
    for g in groupby(df, :CAN_ID)
        append!(filtered_df, filter_dialysis_episode(g))
    end

    return filtered_df
end

"""
    load_donor(filepath::String) -> DataFrame

Load a CSV file containing donor information and return a cleaned `DataFrame`.

### Details
Cleaning steps:
- Remove donors for which the decision status (`DECISION`) is missing.
- Keep only donors attributed to List 5 (Administrative to Transplant Québec).
- Keep only recipients who would have been offered a kidney using `filter_offers_until_last_acceptance`,
- Drop administrative columns not required for downstream processing.
- Replace `WEIGHT = 0`` and `HEIGHT = 0` by `missing`.
"""
function load_donor(filepath::String)
    df_donors = CSV.read(filepath, DataFrame; missingstring=["-", "", "NULL"])

    dropmissing!(df_donors, [:DECISION])

    # Only attribution type 5
    filter!(row -> row.ATT_TYPE == 5, df_donors)

    df = filter_offers_until_last_acceptance(df_donors)

    fill_hla_pairs!(df, "DON")

    col_to_harmonize = [:DON_AGE, :DON_BLOOD, :DON_DEATH_TM, :WEIGHT, :HEIGHT, :ETHNICITY, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :HCV, :DCD]

    for col in col_to_harmonize
        harmonize_col!(df, col=col, idcol=:DON_ID)
    end

    cols_to_drop = [:DON_LAB_DT_TM, :ATT_TYPE]
    select!(df, Not(cols_to_drop))

    allowmissing!(df, [:WEIGHT, :HEIGHT])
    df.WEIGHT[df.WEIGHT .≈ 0.] .= missing
    df.HEIGHT[df.HEIGHT .≈ 0.] .= missing

    return df
end



function load_data(recipient_filepath::AbstractString, cpra_filepath::AbstractString, donor_filepath::AbstractString)

    df_recipients = load_recipient(recipient_filepath)
    df_donors = load_donor(donor_filepath)
    df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])

    G = groupby(df_recipients, :CAN_ID)

    filtered_df_recipients = similar(df_recipients, 0)
    filtered_df_donors = similar(df_donors, 0)

    cpra = Vector{Int64}(undef, length(G))
    time_to_first_offer = Vector{Union{Int64,Missing}}(undef, length(G))
    time_to_transplant = Vector{Union{Int64,Missing}}(undef, length(G))
    active_waiting_proportion = Vector{Union{Float64,Missing}}(undef, length(G))

    exit_status = ["TX", "TX VIVANT", "X", "DCD"]

    for (i, df_recipient) in enumerate(G)

        r = last(sort(df_recipient, :UPDATE_TM))
        push!(filtered_df_recipients, r)

        CAN_ID = r.CAN_ID
        listing_date = r.CAN_LISTING_DT
        exit_date = uppercase(strip(r.OUTCOME)) ∈ exit_status ? r.UPDATE_TM : nothing


        df_cpra_recipient = filter(row -> row.CAN_ID == CAN_ID, df_cpra)

        if isempty(df_cpra_recipient)
            cpra[i] = 0
        else
            target = r.UPDATE_TM
            nearest_index = argmin(abs.(df_cpra_recipient.UPDATE_TM .- target))
            cpra[i] = Int(df_cpra_recipient.CAN_CPRA[nearest_index])
        end


        if uppercase(strip(r.OUTCOME)) == "TX"
            time_to_transplant[i] = Dates.days(r.UPDATE_TM - r.CAN_LISTING_DT)
        end

        active_waiting_proportion[i] = KidneyAllocation.recipient_active_waiting_proportion(df_recipient)


        df_offers = filter(row -> row.CAN_ID == CAN_ID && row.DON_DEATH_TM ≥ listing_date, df_donors)

        if !isnothing(exit_date)
            filter!(row -> row.DON_DEATH_TM ≤ exit_date, df_offers)
        end

        # Filtre les df_donors pour ne garder que les lignes qui correspondent à la bonne date de dialyse des receveurs
        append!(filtered_df_donors, df_offers)

        if !isempty(df_offers)
            time_to_first_offer[i] = Dates.days(minimum(df_offers.DON_DEATH_TM) - listing_date)
        else
            time_to_first_offer[i] = missing
        end

    end

    filtered_df_recipients.CPRA = cpra
    filtered_df_recipients.active_waiting_proportion = active_waiting_proportion
    filtered_df_recipients.time_to_first_offer = time_to_first_offer
    filtered_df_recipients.time_to_transplant = time_to_transplant

    # Si time_to_first_offer == missing mais qu'on a time_to_transplant != missing, time_to_first_offer est fixé à time_to_transplant
    filtered_df_recipients[!, :time_to_first_offer] = coalesce.(filtered_df_recipients[!, :time_to_first_offer], filtered_df_recipients[!, :time_to_transplant])

    return filtered_df_recipients, filtered_df_donors

end



df_recipients, df_donors = load_data(recipient_filepath, cpra_filepath, donor_filepath)

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
function recipient_from_row(r::DataFrameRow)

    blood = parse_abo(String(r.CAN_BLOOD))
    expiration_date = get_expiration_date(r)

    return Recipient(
        r.CAN_BTH_DT, r.CAN_DIAL_DT, r.CAN_LISTING_DT, blood,
        r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
        r.CPRA,
        expiration_date=expiration_date,
        active_waiting_proportion=r.active_waiting_proportion,
    )
end

"""
    get_expiration_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

Return the time of the latest retained update if its outcome indicates a
permanent exit from the waiting list; otherwise return `nothing`.
"""
function get_expiration_date(r::DataFrameRow)::Union{DateTime,Nothing}
    expiration_outcomes = ("X", "TX VIVANT", "DCD")

    return uppercase(strip(r.OUTCOME)) ∈ expiration_outcomes ? r.UPDATE_TM : nothing
end

r = eachrow(filtered_df_recipients)[4]

get_expiration_date(r)
recipient_from_row(r)


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

r = eachrow(filtered_df_donors)[2000]

donor_from_row(r)







function build_donor_registry(df_donors::AbstractDataFrame)

    kdri_required = Symbol[
        :DON_ID, :DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD
    ]
    donor_required = Symbol[
        :DON_DEATH_TM, :DON_AGE, :DON_BLOOD, :DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2
    ]

    filtered_df_donors = dropmissing(df_donors, kdri_required)
    dropmissing!(filtered_df_donors, donor_required)

    donor_registry = Dict{Int,Donor}()
    given_kidneys = Dict{Int,Int}()

    for g in groupby(filtered_df_donors, :DON_ID)
        r = first(g)
        donor = donor_from_row(r)
        n = count(g.DECISION .== "Acceptation")

        donor_registry[r.DON_ID] = donor
        given_kidneys[r.DON_ID] = n == 0 ? 1 : n  # Parfois, le nombre de reins est égal à 0 car on a retiré le receveur qui la reçu en raison des multiples inscriptions 
    end

    return donor_registry, given_kidneys
end


d, g = build_donor_registry(filtered_df_donors)

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
function build_recipient_registry(df_recipients::AbstractDataFrame)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]
    
    df = dropmissing(df_recipients, required)

    recipient_registry = Dict{Int,Recipient}()
    time_to_first_offer_registry = Dict{Int,Union{Int, Missing}}()
    time_to_transplant_registry = Dict{Int,Union{Int, Missing}}()

    for r in eachrow(df)
        recipient_registry[r.CAN_ID] = recipient_from_row(r)
        time_to_first_offer_registry[r.CAN_ID] = r.time_to_first_offer
        time_to_transplant_registry[r.CAN_ID] = r.time_to_transplant
    end

    return recipient_registry, time_to_first_offer_registry, time_to_transplant_registry
end

d, t1, t2 = build_recipient_registry(df_recipients)

