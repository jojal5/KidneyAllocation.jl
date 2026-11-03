
# Function for Recipients and Donors

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



# Functions for Recipents

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
    filter_dialysis_episode(df) -> DataFrame

Return records from the recorded dialysis time through the first exit status
on or after dialysis, including the exit record. The listing date must lie
within the retained episode; otherwise it is replaced by the dialysis date.
"""
function filter_dialysis_episode(df::AbstractDataFrame)::DataFrame
    check_df_columns(df, :CAN_ID, :CAN_DIAL_DT, :CAN_LISTING_DT, :OUTCOME, :UPDATE_TM)
    check_df_column_constant(df, :CAN_ID, :CAN_DIAL_DT, :CAN_LISTING_DT)

    isempty(df) && return DataFrame(df)

    dialysis_date = first(df.CAN_DIAL_DT)
    listing_date = first(df.CAN_LISTING_DT)
    dialysis_time = DateTime(dialysis_date)
    exit_time = first_exit_after_dialysis(df)

    episode_df = DataFrame(filter(df) do row
        update_time = DateTime(row.UPDATE_TM)

        update_time ≥ dialysis_time &&
            (isnothing(exit_time) || update_time ≤ exit_time)
    end)

    isempty(episode_df) && return episode_df

    episode_end = isnothing(exit_time) ?
                  maximum(DateTime.(episode_df.UPDATE_TM)) :
                  exit_time

    listing_time = DateTime(listing_date)

    if !(dialysis_time ≤ listing_time ≤ episode_end)
        corrected_listing_date = convert(typeof(listing_date), dialysis_date)
        episode_df.CAN_LISTING_DT .= corrected_listing_date
    end

    # Removing duplicated rows (happens for recipient 501, 827 and 2051)
    return unique!(episode_df)
end

"""
    first_exit_after_dialysis(df) -> Union{DateTime,Nothing}

Return the first exit time on or after the recorded dialysis time for a
single recipient, or `nothing` if no such exit exists.
"""
function first_exit_after_dialysis(df::AbstractDataFrame)::Union{DateTime,Nothing}

    check_df_columns(df, :CAN_ID, :CAN_DIAL_DT, :OUTCOME, :UPDATE_TM)
    check_df_column_constant(df, :CAN_ID, :CAN_DIAL_DT)

    isempty(df) && return nothing

    dialysis_time = DateTime(first(df.CAN_DIAL_DT))
    exit_outcomes = ("X", "TX VIVANT", "DCD", "TX")
    exit_time = nothing

    for row in eachrow(df)
        update_time = DateTime(row.UPDATE_TM)
        outcome = uppercase(strip(string(row.OUTCOME)))

        update_time ≥ dialysis_time && outcome ∈ exit_outcomes || continue

        exit_time = isnothing(exit_time) ? update_time : min(exit_time, update_time)
    end

    return exit_time
end

"""
    filter_first_waitlist_episode(df_recipient) -> DataFrame

Return records from the first waiting-list episode of one recipient, ordered by
`:UPDATE_TM`. Rows after the first exit outcome are removed, while the exit
record itself is retained. If no exit outcome is recorded, all rows are kept.
"""
function filter_first_waitlist_episode(df_recipient::AbstractDataFrame)::DataFrame
    check_df_column_constant(df_recipient, :CAN_ID)
    check_df_columns(df_recipient, :OUTCOME, :UPDATE_TM)

    isempty(df_recipient) && return DataFrame(df_recipient)

    exit_outcomes = ("TX", "TX VIVANT", "X", "DCD")
    is_exit(outcome) = uppercase(strip(string(outcome))) ∈ exit_outcomes

    sorted_df = sort(df_recipient, :UPDATE_TM)
    exit_index = findfirst(is_exit, sorted_df.OUTCOME)
    last_index = something(exit_index, nrow(sorted_df))

    if last_index < nrow(sorted_df)
        if sorted_df.CAN_LISTING_DT[last_index] > sorted_df.UPDATE_TM[last_index]
            sorted_df.CAN_LISTING_DT .= missing
        end
        if sorted_df.CAN_DIAL_DT[last_index] > sorted_df.UPDATE_TM[last_index]
            sorted_df.CAN_DIAL_DT .= missing
        end
    end

    return reverse(sorted_df[1:last_index, :])
end

"""
    get_exit_date(r::DataFrameRow) -> Union{DateTime,Nothing}

Return the time of the exit from the waiting list for a single
recipient, or `nothing` if no exit outcome is recorded.
"""
function get_exit_date(r::DataFrameRow)::Union{DateTime,Nothing}
    exit_outcomes = ("TX", "TX VIVANT", "X", "DCD")

    date = r.UPDATE_TM
    outcome = uppercase(strip(r.OUTCOME))

    return outcome ∈ exit_outcomes ? date : nothing
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

    # # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
    # dropmissing!(df, [:CAN_DIAL_DT])

    # If listing time is missing, replacing it by the dialysis time.
    # coalesce_listing!(df)

    # Keeping only adult recipients
    # filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

    # Keeping only the first outcome after the registered dialysis date for each recipient. 
    # Some recipients are registered more than once and key information are lost for the additional registrations.
    # We keep the history corresponding to the recorded dialysis date.
    # filtered_df = similar(df, 0)
    # for g in groupby(df, :CAN_ID)
    #     append!(filtered_df, filter_dialysis_episode(g))
    # end

    # return filtered_df

    return df
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
    check_df_column_constant(df, :CAN_ID, :CAN_LISTING_DT)
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


# Function for Donors

"""
    filter_offers_until_last_acceptance(df_donor::AbstractDataFrame) -> DataFrame

For one donor, sort offers by decreasing `DON_CAN_SCORE` and retain rows
through the last `"ACCEPTATION"` decision, inclusive. If no offer is accepted,
retain all rows.
"""
function filter_offers_until_last_acceptance(df_donor::AbstractDataFrame)::DataFrame
    check_df_columns(df_donor, :DON_ID, :DON_CAN_SCORE, :DECISION)

    isempty(df_donor) && return DataFrame(df_donor)

    check_df_column_constant(df_donor, :DON_ID)

    ranked_df = sort(df_donor, :DON_CAN_SCORE, rev=true)

    last_acceptance = findlast(isequal("ACCEPTATION"), uppercase.(strip.(ranked_df.DECISION)))

    if isnothing(last_acceptance)
        return ranked_df
    else
        return first(ranked_df, last_acceptance)
    end
end

"""
    infer_number_of_kidneys_offered(df_donor::AbstractDataFrame) -> Int

Infer the number of kidneys offered from a single donor history. Return the
number of recorded `"TX"` status, or `2` if no transplantation is recorded.
In the latter, it is assumed that it was two kidneys that were refused by the recipients.
"""
function infer_number_of_kidneys_offered(df_donor::AbstractDataFrame)::Int

    check_df_column_constant(df_donor, :DON_ID)

    :STATUS ∈ propertynames(df_donor) ||
        throw(ArgumentError("Missing column :STATUS"))

    status = uppercase.(strip.(skipmissing(df_donor.STATUS)))

    n_transplanted = count(status .== "TX")

    n_offered = n_transplanted > 0 ? n_transplanted : 2

    return n_offered
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

    # dropmissing!(df_donors, [:DECISION])

    # Only attribution type 5
    # filter!(row -> row.ATT_TYPE == 5, df_donors)

    fill_hla_pairs!(df_donors, "DON")

    col_to_harmonize = [:DON_AGE, :DON_BLOOD, :DON_DEATH_TM, :WEIGHT, :HEIGHT, :ETHNICITY, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :HCV, :DCD]

    for col in col_to_harmonize
        harmonize_col!(df_donors, col=col, idcol=:DON_ID)
    end

    # cols_to_drop = [:DON_LAB_DT_TM, :ATT_TYPE]
    cols_to_drop = [:DON_LAB_DT_TM]
    select!(df_donors, Not(cols_to_drop))

    allowmissing!(df_donors, [:WEIGHT, :HEIGHT])
    df_donors.WEIGHT[df_donors.WEIGHT .≈ 0.] .= missing
    df_donors.HEIGHT[df_donors.HEIGHT .≈ 0.] .= missing

    # filtered_df = similar(df_donors, 0)
    # for g in groupby(df_donors, :DON_ID)
    #     append!(filtered_df, filter_offers_until_last_acceptance(g))
    # end

    # return filtered_df

    return df_donors
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

# Load preprocessed data


"""
    load_data(recipient_filepath::AbstractString, cpra_filepath::AbstractString, donor_filepath::AbstractString)

Load and clean the recipient and donor database

# Details
- To avoid errors in the recipient raw file, recipients are reconstruted using the donor file. A recipient is retrieve if they has received at least an offer from List 5.
- Si le patient est transplanté, sa sortie s'effectue au moment du décès du donneur.
- Si le patient n'est pas transplanté, sa date de dernière mise à jour corespond à la date du dernier score sur la liste 5 calculé dans le fichier (ça se peut que le patient n'ait jamais eu l'opportunité d'accepter ou non).
- Le dataframe des recveurs contient le flag `SCORE_MATCHING_FLAG`. 
  La valeur 1 indique que le score calculé avec les données du receveur correspond au score enregistré dans le fichier, suggérant une absence d'erreur dans les données du receveur.
  La valeur 0 indique que le score calculé avec les données du receveur ne correspond pas au score enregistré dans le fichier, suggérant une erreur dans les données du receveur.
  La valeur missing indique que des données sont manquantes pour le calcul du score.
"""
function load_data(recipient_filepath::AbstractString, cpra_filepath::AbstractString, donor_filepath::AbstractString)

    df_recipients = KidneyAllocation.load_recipient(recipient_filepath)
    df_donors = KidneyAllocation.load_donor(donor_filepath)
    df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])


    # Quelques observations
    # Le patient 7181, on a que la transplantation est enregistrée le 2017-07-14 dans le fichiers des receveurs, mais le donneur est décédé le 2016-01-14
    # Le patient CAN_ID = 142 a reçu plusieurs offres de la liste 5, mais a accepté une offre de la liste 10. Il était compétitif jusqu'à ce moment.
    # Le patient CAN_ID = 14183 a un CPRA de 95, alors c'est pour ça que sa date de dialyse est inconnue.
    # Le patient 7387 revient sur la liste tout de suite après sa TX 
    # Le patient 284 est bon exemple où toutes les données sont mélangées lors de la deuxième inscription


    # Récupération des candidats dont au moins un score a été calculé pour la liste 5
    # Reconstruction des receveurs par le fichier des donneurs

    is_tx(outcome::Union{Missing,AbstractString}) = uppercase(strip(string(outcome))) == "TX"

    filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
    filter!(row -> row.DON_ID ≠ 0, filtered_donors)

    df_offers_recipients = similar(filtered_donors, 0)
    offered_recipients = unique(df_donors.CAN_ID)

    filtered_recipients = DataFrame(CAN_ID=Int[], CAN_BTH_DT=Date[], CAN_BLOOD=String[],
        CAN_A1=Int[], CAN_A2=Int[], CAN_B1=Int[], CAN_B2=Int[], CAN_DR1=Int[], CAN_DR2=Int[],
        CAN_LISTING_DT=DateTime[], CAN_DIAL_DT=Date[],
        OUTCOME=String[], UPDATE_TM=DateTime[],
        CPRA=Int[], active_waiting_proportion=Float64[])

    allowmissing!(filtered_recipients)


    for CAN_ID in offered_recipients

        df_offers_recipient = sort(filter(row -> row.CAN_ID == CAN_ID, filtered_donors), :DON_DEATH_TM)

        df_list5 = filter(row -> row.ATT_TYPE == 5, df_offers_recipient)

        # Skip if patient has no history on List 5
        if isempty(df_list5)
            continue
        end

        # Récupération de la ligne correspondante à la première transplantation (liste 5)
        ind = findfirst(is_tx, skipmissing(df_list5.STATUS))

        # Si receveur non transplanté, on récupére la ligne de la dernière offre (liste 5)
        if isnothing(ind)
            ind = nrow(df_list5)
        end

        target_date = df_list5.DON_DEATH_TM[ind]

        # Récupération des informations sur le receveur
        df_recipient = filter(row -> row.CAN_ID == CAN_ID, df_recipients)
        df_recipient_first_registration = filter_first_waitlist_episode(df_recipient)
        p = KidneyAllocation.recipient_active_waiting_proportion(df_recipient_first_registration)

        # Récupération du CPRA au moment de la target_date
        df_cpra_recipient = sort(filter(row -> row.CAN_ID == CAN_ID, df_cpra), :UPDATE_TM, rev=true)
        latest_index = findfirst(t -> t ≤ target_date, df_cpra_recipient.UPDATE_TM)
        cpra = isnothing(latest_index) ? 0 : Int(df_cpra_recipient.CAN_CPRA[latest_index])

        registered_outcome = df_recipient_first_registration.OUTCOME[1]

        if is_tx(df_list5.STATUS[ind])
            outcome = "TX"
        elseif registered_outcome == "TX"
            outcome = "X"
        else
            outcome = registered_outcome
        end

        r = first(df_recipient_first_registration)

        push!(filtered_recipients,
            [r.CAN_ID, r.CAN_BTH_DT, r.CAN_BLOOD,
                r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
                r.CAN_LISTING_DT, r.CAN_DIAL_DT,
                outcome, target_date,
                cpra, p])

    end

    # Check registered score against computed score
    # A difference in score indicates an error in the data

    df_scores = DataFrame(CAN_ID=Int[], DON_ID=Int[], REGISTERED_SCORE=Float64[], COMPUTED_SCORE=Float64[])
    allowmissing!(df_scores, :COMPUTED_SCORE)

    for r in eachrow(filtered_recipients)

        df_offers = sort(filter(row -> row.CAN_ID == r.CAN_ID, df_donors), :DON_DEATH_TM)

        ind = findfirst(is_tx, df_offers.STATUS)
        if isnothing(ind)
            ind = findlast(df_offers.DON_DEATH_TM .≤ r.UPDATE_TM)
        end

        registered_score = df_offers.DON_CAN_SCORE[ind]

        computed_score = 
            try
                recipient = KidneyAllocation.recipient_from_row(r)
                donor = KidneyAllocation.donor_from_row(df_offers[ind, :])

                KidneyAllocation.score(donor, recipient)
            catch err
                if err isa MethodError
                    missing
                else
                    rethrow()
                end
            end

        push!(df_scores, [r.CAN_ID, df_offers.DON_ID[ind], registered_score, computed_score])

    end

    score_diff = abs.(df_scores.COMPUTED_SCORE - df_scores.REGISTERED_SCORE)
    # count(skipmissing(score_diff .> .1))
    # count(skipmissing(score_diff .≤ .1))
    # count(ismissing.(score_diff))

    score_matching_flag = Vector{Union{Bool,Missing}}(undef, length(score_diff))

    for i in eachindex(score_diff)
        if ismissing(score_diff[i])
            score_matching_flag[i] = missing
        else
            score_matching_flag[i] = score_diff[i] ≤ 0.1
        end
    end

    filtered_recipients.SCORE_MATCHING_FLAG = score_matching_flag

    # Récupération des donneurs de la liste 5

    df_donors_list5 = similar(filtered_donors, 0)

    for df_offer in groupby(filtered_donors, :DON_ID)

        ind = findall(is_tx, df_offer.STATUS)
        n = count(df_offer.ATT_TYPE[ind] .== 5)

        if n > 0
            append!(df_donors_list5, df_offer[ind, :])
        end

    end

    return filtered_recipients, df_donors_list5

end


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
        given_kidneys[r.DON_ID] = infer_number_of_kidneys_offered(g)
        # given_kidneys[r.DON_ID] = n == 0 ? 1 : n  # Parfois, le nombre de reins est égal à 0 car on a retiré le receveur qui la reçu en raison des multiples inscriptions 
    end

    return donor_registry, given_kidneys
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
function build_recipient_registry(df_recipients::AbstractDataFrame)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]

    df = dropmissing(df_recipients, required)

    recipient_registry = Dict{Int,Recipient}()
    time_to_first_offer_registry = Dict{Int,Union{Int,Missing}}()
    time_to_transplant_registry = Dict{Int,Union{Int,Missing}}()

    for r in eachrow(df)
        recipient_registry[r.CAN_ID] = recipient_from_row(r)
        time_to_first_offer_registry[r.CAN_ID] = r.time_to_first_offer
        time_to_transplant_registry[r.CAN_ID] = r.time_to_transplant
    end

    return recipient_registry, time_to_first_offer_registry, time_to_transplant_registry
end

"""
    retrieve_observed_waiting_list(recipient_filepath, date) -> Vector{Int64}

Return the IDs of recipients registered on the waiting list at `date`, based
on the recipient history stored in `recipient_filepath`.
"""
function retrieve_observed_waiting_list(df_recipients::AbstractDataFrame, date::DateTime)::Vector{Int64}

    registered_ids = Int64[]

    for r in eachrow(df_recipients)
        arrival_time = DateTime(r.CAN_LISTING_DT)
        exit_time = get_exit_date(r)

        if arrival_time ≤ date && (isnothing(exit_time) || date ≤ exit_time)
            push!(registered_ids, r.CAN_ID)
        end
    end

    return registered_ids
end

retrieve_observed_waiting_list(df_recipients::AbstractDataFrame, date::Date) = retrieve_observed_waiting_list(df_recipients, DateTime(date))













# """
#     get_last_update(df::AbstractDataFrame) -> Tuple{DateTime, String}

# Return the time and outcome of the most recent update retained by
# `filter_outcomes` for a single recipient. 
# """
# function get_last_update(df::AbstractDataFrame)::Tuple{DateTime, String}

#     isempty(df) && throw(ArgumentError("DataFrame is empty"))

#     check_df_columns(df, :OUTCOME, :UPDATE_TM)
#     check_df_column_constant(df, :CAN_ID)

#     filtered_df = filter_outcomes(df)
#     last_index = argmax(filtered_df.UPDATE_TM)

#     return (filtered_df.UPDATE_TM[last_index], String(filtered_df.OUTCOME[last_index]))
# end

# """
#     get_expiration_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

# Return the time of the latest retained update if its outcome indicates a
# permanent exit from the waiting list; otherwise return `nothing`.
# """
# function get_expiration_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
#     expiration_outcomes = ("X", "TX VIVANT", "DCD")

#     date, outcome = get_last_update(df)

#     return uppercase(strip(outcome)) ∈ expiration_outcomes ? date : nothing
# end

# """
#     get_exit_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

# Return the time of the first recorded exit from the waiting list for a single
# recipient, or `nothing` if no exit outcome is recorded.
# """
# function get_exit_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
#     exit_outcomes = ("TX", "TX VIVANT", "X", "DCD")

#     date, outcome = get_last_update(df)

#     return uppercase(strip(outcome)) ∈ exit_outcomes ? date : nothing
# end

# """
#     first_exit_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

# Return the earliest update time associated with an exit outcome for a single
# recipient, or `nothing` if no exit outcome is recorded.
# """
# function first_exit_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
#     isempty(df) && return nothing

#     check_df_column_constant(df, :CAN_ID)
#     check_df_columns(df, :OUTCOME, :UPDATE_TM)

#     exit_outcomes = ("X", "TX VIVANT", "DCD", "TX")
#     earliest_exit_time = nothing

#     for row in eachrow(df)
#         outcome = uppercase(strip(string(row.OUTCOME)))

#         if outcome ∈ exit_outcomes
#             update_time = DateTime(row.UPDATE_TM)

#             earliest_exit_time = isnothing(earliest_exit_time) ?
#                 update_time :
#                 min(earliest_exit_time, update_time)
#         end
#     end

#     return earliest_exit_time
# end

# """
#     filter_outcomes(df::AbstractDataFrame) -> DataFrame

# Return records from the first waiting-list episode of a single recipient.

# Rows after the first exit outcome are removed, while the exit record is
# retained. Exact duplicate retained rows are removed.
# (Happen for recipients 501, 827 and 2051)

# If every retained `UPDATE_TM` precedes the recorded `CAN_LISTING_DT`,
# `CAN_LISTING_DT` is replaced by `CAN_DIAL_DT` as a proxy for the unavailable
# original listing date.

# `:UPDATE_TM` must contain `DateTime` values so that outcomes occurring on the
# same calendar date can be ordered correctly.
# """
# function filter_outcomes(df::AbstractDataFrame)::DataFrame
#     isempty(df) && return DataFrame(df)

#     check_df_column_constant(df, :CAN_ID, :CAN_LISTING_DT, :CAN_DIAL_DT)
#     check_df_columns(df, :OUTCOME, :UPDATE_TM)

#     update_type = eltype(skipmissing(df[!, :UPDATE_TM]))

#     update_type <: DateTime || throw(ArgumentError("Column :UPDATE_TM must contain DateTime values; got $update_type"))

#     exit_time = first_exit_date(df)

#     if exit_time === nothing
#         filtered_df = DataFrame(df)
#     else
#         filtered_df = filter(:UPDATE_TM => date -> date ≤ exit_time, df)

#         listing_date = first(filtered_df.CAN_LISTING_DT)
#         dialysis_date = convert(typeof(listing_date), first(filtered_df.CAN_DIAL_DT))

#         if maximum(filtered_df.UPDATE_TM) < listing_date
#             filtered_df.CAN_LISTING_DT .= dialysis_date
#         end

#     end

#     # Removing duplicated rows (happens for recipient 501, 827 and 2051)
#     unique!(filtered_df)

#     return filtered_df
# end


# """
#     first_exit_after_dialysis(df) -> Union{DateTime,Nothing}

# Return the first exit time on or after the recorded dialysis time for a
# single recipient, or `nothing` if no such exit exists.
# """
# function first_exit_after_dialysis(
#     df::AbstractDataFrame,
# )::Union{DateTime,Nothing}

#     check_df_columns(df, :CAN_ID, :CAN_DIAL_DT, :OUTCOME, :UPDATE_TM)
#     check_df_column_constant(df, :CAN_ID, :CAN_DIAL_DT)

#     isempty(df) && return nothing

#     dialysis_time = DateTime(first(df.CAN_DIAL_DT))
#     exit_outcomes = ("X", "TX VIVANT", "DCD", "TX")
#     exit_time = nothing

#     for row in eachrow(df)
#         update_time = DateTime(row.UPDATE_TM)
#         outcome = uppercase(strip(string(row.OUTCOME)))

#         update_time ≥ dialysis_time && outcome ∈ exit_outcomes || continue

#         exit_time = isnothing(exit_time) ? update_time : min(exit_time, update_time)
#     end

#     return exit_time
# end







# """
#     load_donor(filepath::String) -> DataFrame

# Load a CSV file containing donor information and return a cleaned `DataFrame`.

# ### Details
# Cleaning steps:
# - Remove donors for which the decision status (`DECISION`) is missing.
# - Keep only donors attributed to List 5 (Administrative to Transplant Québec).
# - Keep only recipients who would have been offered a kidney using `filter_offers_until_last_acceptance`,
# - Drop administrative columns not required for downstream processing.
# - Replace `WEIGHT = 0`` and `HEIGHT = 0` by `missing`.
# """
# function load_donor(filepath::String)
#     df_donors = CSV.read(filepath, DataFrame; missingstring=["-", "", "NULL"])

#     dropmissing!(df_donors, [:DECISION])

#     # Only attribution type 5
#     filter!(row -> row.ATT_TYPE == 5, df_donors)

#     df = filter_offers_until_last_acceptance(df_donors)

#     fill_hla_pairs!(df, "DON")

#     col_to_harmonize = [:DON_AGE, :DON_BLOOD, :DON_DEATH_TM, :WEIGHT, :HEIGHT, :ETHNICITY, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :HCV, :DCD]

#     for col in col_to_harmonize
#         harmonize_col!(df, col=col, idcol=:DON_ID)
#     end

#     cols_to_drop = [:DON_LAB_DT_TM, :ATT_TYPE]
#     select!(df, Not(cols_to_drop))

#     allowmissing!(df, [:WEIGHT, :HEIGHT])
#     df.WEIGHT[df.WEIGHT.≈0.] .= missing
#     df.HEIGHT[df.HEIGHT.≈0.] .= missing

#     return df
# end

# """
#     donor_from_row(r::DataFrameRow) -> Donor

# Construct a `Donor` from a donor `DataFrameRow`. Assume non-missing value.
# """
# function donor_from_row(r::DataFrameRow)

#     age = r.DON_AGE
#     height = r.HEIGHT
#     weight = r.WEIGHT
#     hypertension = r.HYPERTENSION == 1
#     diabetes = r.DIABETES == 1
#     cva = r.DEATH ∈ (4, 16)
#     creatinine = creatinine_mgdl(r.CREATININE)
#     dcd = r.DCD == 1

#     kdri = evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

#     arrival = Date(r.DON_DEATH_TM)
#     blood = parse_abo(r.DON_BLOOD)

#     donor = Donor(arrival, age, blood, r.DON_A1, r.DON_A2, r.DON_B1, r.DON_B2, r.DON_DR1, r.DON_DR2, kdri)

#     return donor
# end


# """
#     build_donor_registry(filepath::String) -> Dict{Int,Donor}

# Build a dictionary of `Donor` objects keyed by `DON_ID` from a donor CSV file.

# Rows with missing values required to construct a donor and compute KDRI are
# discarded. For each `DON_ID`, one `Donor` is constructed from the corresponding
# row using `donor_from_row`.

# ### Notes
# - The donor arrival date is derived from the donor death date (`DON_DEATH_TM`).
# - Rows with incomplete data required for KDRI computation are discarded.
# """
# function build_donor_registry(filepath::String)::Dict{Int,Donor}
#     df = load_donor(filepath)

#     kdri_required = Symbol[
#         :DON_ID, :DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD
#     ]
#     donor_required = Symbol[
#         :DON_DEATH_TM, :DON_AGE, :DON_BLOOD, :DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2
#     ]

#     dropmissing!(df, kdri_required)
#     dropmissing!(df, donor_required)

#     donor_by_don_id = Dict{Int,Donor}()

#     for g in groupby(df, :DON_ID)
#         r = first(g)
#         don_id = r.DON_ID
#         donor_by_don_id[don_id] = donor_from_row(r)
#     end

#     return donor_by_don_id
# end

# """
#     kidneys_given_by_donor(df_donors) -> Dict{Int,Int}

# Return the number of given kidneys for each `DON_ID`.
# """
# function kidneys_given_by_donor(df_donors::AbstractDataFrame)
#     @assert "DON_ID" in names(df_donors) "Missing column :DON_ID"
#     @assert "STATUS" in names(df_donors) "Missing column :STATUS"

#     kidney_by_don_id = Dict{Int,Int}()

#     for g in groupby(df_donors, :DON_ID)
#        kidney_by_don_id[g.DON_ID[1]] = count(==("TX"), skipmissing(g.STATUS))
#     end

#     return kidney_by_don_id

# end


# """
#     load_recipient(filepath::AbstractString)

# Load the CSV file containing recipient information and return a cleaned `DataFrame`.

# ### Details

# Cleaning steps:
# - Replace the values `24L` and `24Low` with `24` in `CAN_A2`.
# - Keep only kidney transplant recipients (exclude multi-organ outcomes).
# - Remove recipients with missing dialysis date (`CAN_DIAL_DT`).
# - If listing date is missing, replace it by the dialysis date.
# - If listing date is after dialysis date, replace the listing date by the dialysis date.
# - Keeping only adult recipients (18 years old and older)
# - Convert `DateTime` columns to `Date`.
# """
# function load_recipient(filepath::AbstractString)
#     df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

#     # # Manually removing duplicated rows for TX (TODO check in new dataset version if these errors are still present)
#     # if (df.OUTCOME[3104] == df.OUTCOME[3105]) & (df.UPDATE_TM[3104] == df.UPDATE_TM[3105])
#     #     deleteat!(df, 3104)
#     # end

#     # Replacing the value 24L and 24Low with 24 for instance
#     df.CAN_A2 = parse_hla_int.(df.CAN_A2)

#     fill_hla_pairs!(df, "CAN")

#     # Harmonize columns when some rows have missig value for a same recipient
#     col_to_harmonize = [:CAN_BTH_DT, :CAN_BLOOD, :CAN_DIAL_DT]
#     for col in col_to_harmonize
#         harmonize_col!(df, col=col, idcol=:CAN_ID )
#     end

#     # Keeping only the recipients for kidney transplant
#     filter!(row -> uppercase(row.OUTCOME) ∈ ("TX", "1", "0", "X", "DCD", "TX VIVANT"), df)

#     # Removing the recipient where the dialysis time is missing. It happens for recipients that received a kidney from a living donor.
#     dropmissing!(df, [:CAN_DIAL_DT])

#     # If listing time is missing, replacing it by the dialysis time.
#     coalesce_listing!(df)

#     # Keeping only adult recipients
#     filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

#     # Keeping only the first outcome for each recipient. Some recipients are registered more than once and key information are lost for the additional registration
#     filtered_df = similar(df, 0)
#     for g in groupby(df, :CAN_ID)
#         append!(filtered_df, filter_outcomes(g))
#     end

#     # Keeping only the first outcome after the registered dialysis date for each recipient. Some recipients are registered more than once and key information are lost for the additional registration
#     filtered_df = similar(df, 0)
#     for g in groupby(df, :CAN_ID)
#         append!(filtered_df, filter_dialysis_episode(g))
#     end

#     # # If listing time is missing, replacing it by the dialysis time. If listing is before dialysis, set listing = dialysis
#     # enforce_listing_after_dialysis!(df)

#     # Transform DateTime in Date
#     # df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
#     # df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
#     # df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

#     return filtered_df
# end

# """
#     recipient_from_row(
#         r::DataFrameRow,
#         cpra::Integer=0;
#         expiration_date=nothing,
#         active_waiting_proportion=1.0,
#     ) -> Recipient

# Construct a `Recipient` from one row of a recipient's history. Required values
# in `r` are assumed non-missing.
# """
# function recipient_from_row(
#     r::DataFrameRow,
#     cpra::Integer=0;
#     expiration_date::Union{Date,DateTime,Nothing}=nothing,
#     active_waiting_proportion::Real=1.0,
# )
#     birth = r.CAN_BTH_DT
#     dialysis = r.CAN_DIAL_DT
#     arrival = r.CAN_LISTING_DT

#     blood = parse_abo(String(r.CAN_BLOOD))

#     return Recipient(
#         birth, dialysis, arrival, blood,
#         r.CAN_A1, r.CAN_A2,
#         r.CAN_B1, r.CAN_B2,
#         r.CAN_DR1, r.CAN_DR2,
#         cpra;
#         expiration_date=expiration_date,
#         active_waiting_proportion=active_waiting_proportion,
#     )
# end

# """
#     build_last_cpra_registry(cpra_filepath::String) -> Dict{CAN_ID, Int64}

# Build a dictionary mapping each recipient (`CAN_ID`) to their most recent calculated Panel Reactive Antibody (cPRA) value.

# ### Details
# - The CPRA file is read from `cpra_filepath`.
# - Rows with missing `CAN_ID` or `CAN_CPRA` are discarded.
# - If a patient has many CPRA values, the most recent is retained (by update time according to UPDATE_TM)
# - CPRA values are rounded to the nearest integer and returned as `Int64`.

# ### Returns
# A dictionary mapping `CAN_ID` to the most recent CPRA value.
# """
# function build_last_cpra_registry(cpra_filepath::String)

#     df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])

#     ("CAN_ID" in names(df_cpra)) || throw(ArgumentError("Missing column :CAN_ID in CPRA file"))
#     ("UPDATE_TM" in names(df_cpra)) || throw(ArgumentError("Missing column :UPDATE_TM in CPRA file"))
#     ("CAN_CPRA" in names(df_cpra)) || throw(ArgumentError("Missing column :CAN_CPRA in CPRA file"))

#     # Drop rows missing the essentials
#     dropmissing!(df_cpra, [:CAN_ID, :CAN_CPRA, :UPDATE_TM])

#     cpra_by_CAN_ID = Dict{eltype(df_cpra.CAN_ID),Int64}()

#     for g in groupby(df_cpra, :CAN_ID)
#         sort!(g, :UPDATE_TM, rev=true)  # most recent first
#         cpra_by_CAN_ID[g.CAN_ID[1]] = Int64(round(g.CAN_CPRA[1]))
#     end

#     return cpra_by_CAN_ID
# end


# """
#     build_recipient_registry(recipient_filepath, cpra_filepath) -> Dict{Int,Recipient}

# Build a dictionary of `Recipient` objects keyed by `CAN_ID`.

# ### Arguments
# - `recipient_filepath::String`: Path to the recipient CSV file (loaded by `load_recipient`).
# - `cpra_filepath::String`: Path to the CPRA CSV file used to compute the latest CPRA per `CAN_ID`.

# ### Notes
# - The expiration date is inferred from the status history and passed as
#   `expiration_date=...` when constructing the `Recipient`.
# - Rows with missing required fields are discarded.
# """
# function build_recipient_registry(recipient_filepath::String, cpra_filepath::String)::Dict{Int,Recipient}
#     df = load_recipient(recipient_filepath)

#     required = Symbol[
#         :CAN_ID, :UPDATE_TM, :OUTCOME,
#         :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
#         :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
#     ]
#     dropmissing!(df, required)

#     cpra_by_can_id = build_last_cpra_registry(cpra_filepath)
#     recipient_by_can_id = Dict{Int,Recipient}()

#     for g in groupby(df, :CAN_ID)

#         filtered_df = filter_outcomes(g)

#         # exp_date = infer_recipient_expiration_date(filtered_df)
#         exp_date = get_expiration_date(filtered_df)
#         active_waiting_proportion = recipient_active_waiting_proportion(filtered_df)

#         r = first(filtered_df)
#         can_id = r.CAN_ID

#         cpra = get(cpra_by_can_id, can_id, 0)
#         recipient_by_can_id[can_id] = recipient_from_row(r, cpra; expiration_date = exp_date, active_waiting_proportion = active_waiting_proportion)
#     end

#     return recipient_by_can_id
# end




















# """
#     check_df_at_most_one_exit_outcome(df, outcomes=...) -> Nothing

# Verify that a single recipient history contains at most one of the specified
# exit outcomes. Outcome codes are matched case-insensitively.
# """
# function check_df_at_most_one_exit_outcome(
#     df::AbstractDataFrame,
#     outcomes::AbstractVector{<:AbstractString}=[
#         "X", "TX VIVANT", "DCD", "TX",
#     ],
# )
#     check_df_columns(df, :OUTCOME)
#     check_df_column_constant(df, :CAN_ID)

#     exit_outcomes = Set(
#         uppercase(strip(string(outcome))) for outcome in outcomes
#     )

#     n_exit_outcomes = count(
#         row -> uppercase(strip(string(row.OUTCOME))) ∈ exit_outcomes,
#         eachrow(df),
#     )

#     n_exit_outcomes ≤ 1 ||
#         throw(ArgumentError(
#             "More than one exit outcome found for this recipient",
#         ))

#     return nothing
# end