
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
    check_df_columns(df, :CAN_DIAL_DT)

    :CAN_LISTING_DT ∈ propertynames(df) ||
        throw(ArgumentError("Missing column :CAN_LISTING_DT"))

    df[!, :CAN_LISTING_DT] = coalesce.(
        df.CAN_LISTING_DT,
        DateTime.(df.CAN_DIAL_DT),
    )

    return df
end

"""
    count_recipient_arrivals(df, start_time, end_time) -> Int

Count recipients whose listing time lies in `[start_time, end_time)`.
"""
function count_recipient_arrivals(df::AbstractDataFrame, start_time::Date, end_time::Date)::Int
    check_df_columns(df, :CAN_ID, :CAN_LISTING_DT)

    start_time ≤ end_time ||
        throw(ArgumentError("`start_time` must not be after `end_time`"))

    n_recipients = 0

    for recipient_df in groupby(df, :CAN_ID)
        listing_time = first(recipient_df.CAN_LISTING_DT)

        if start_time ≤ listing_time < end_time
            n_recipients += 1
        end
    end

    return n_recipients
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
        if !ismissing(sorted_df.CAN_DIAL_DT[last_index]) && sorted_df.CAN_DIAL_DT[last_index] > sorted_df.UPDATE_TM[last_index]
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
# - Keeping only adult recipients (18 years old and older).
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

    # # Keeping only adult recipients
    # filter!(row -> years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df)

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
    count_donor_arrivals(df, start_time, end_time) -> Int

Count donors whose death time lies in `[start_time, end_time)`.
"""
function count_donor_arrivals(df::AbstractDataFrame, start_time::Date, end_time::Date)::Int
    check_df_columns(df, :DON_ID, :DON_DEATH_TM)

    start_time ≤ end_time ||
        throw(ArgumentError("`start_time` must not be after `end_time`"))

    n_donors = 0

    for donor_df in groupby(df, :DON_ID)
        donor_arrival_time = first(donor_df.DON_DEATH_TM)

        if start_time ≤ donor_arrival_time < end_time
            n_donors += 1
        end
    end

    return n_donors
end

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

        # Si la date de dialyse est manquante, on ne retient pas le receveur
        if any(ismissing.(df_recipient_first_registration.CAN_DIAL_DT))
            continue
        else
            # If listing date is missing, or later than dialysis, replace it by dialysis date
            coalesce_listing!(df_recipient_first_registration)

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

        # ind = findall(is_tx, df_offer.STATUS)
        df = filter(row -> is_tx(row.STATUS), df_offer) 

        ind_not_list5 = findall(df.ATT_TYPE .!= 5)
        if !isempty(ind_not_list5)
            id_to_remove = df.CAN_ID[ind_not_list5]
            filter!(row -> row.CAN_ID ∉ id_to_remove, df)
        end

        if nrow(df) > 0
            append!(df_donors_list5, df)
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
        n = count(uppercase.(strip.(g.STATUS)) .== "TX")

        donor_registry[r.DON_ID] = donor
        given_kidneys[r.DON_ID] = n
        # given_kidneys[r.DON_ID] = infer_number_of_kidneys_offered(g)
        # given_kidneys[r.DON_ID] = n == 0 ? 1 : n  # Parfois, le nombre de reins est égal à 0 car on a retiré le receveur qui la reçu en raison des multiples inscriptions 
    end

    return donor_registry, given_kidneys
end



function build_recipient_registry(df_recipients::AbstractDataFrame)

    required = Symbol[
        :CAN_ID, :UPDATE_TM, :OUTCOME,
        :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_LISTING_DT,
        :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2
    ]

    df = dropmissing(df_recipients, required)

    recipient_registry = Dict{Int,Recipient}()

    for r in eachrow(df)
        recipient_registry[r.CAN_ID] = recipient_from_row(r)
    end

    return recipient_registry
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


