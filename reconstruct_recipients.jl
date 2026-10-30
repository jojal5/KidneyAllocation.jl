using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Gadfly

using KidneyAllocation

## Fonctions utiles pas encore introduites dans la librairie

import KidneyAllocation: check_df_columns, check_df_column_constant, parse_abo

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

    return reverse(sorted_df[1:last_index, :])
end


## Chargement des données brutes

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients = KidneyAllocation.load_recipient(recipient_filepath)
df_donors = KidneyAllocation.load_donor(donor_filepath)
df_cpra = CSV.read(cpra_filepath, DataFrame; missingstring=["-", "", "NULL"])


# Quelques observations
# Le patient 7181, on a que la transplantation est enregistrée le 2017-07-14 dans le fichiers des receveurs, mais le donneur est décédé le 2016-01-14
# Le patient CAN_ID = 142 a reçu plusieurs offres de la liste 5, mais a accepté une offre de la liste 10. Il était compétitif jusqu'à ce moment.
# Le patient CAN_ID = 14183 a un CPRA de 95, alors c'est pour ça que sa date de dialyse est inconnue.
# Le patient 7387 revient sur la liste tout de suite après sa TX 
# Le patient 284 est bon exemple où toutes les données sont mélangées lors de la deuxième inscription



## Récupération des candidats dont au moins un score a été calculé pour la liste 5 entre 2016 et 2020
# Reconstruction des recveurs par le fichier des donneurs

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)
G = groupby(filtered_donors, :DON_ID)

df_offers_recipients = similar(filtered_donors, 0)
given_kidneys_registry = Dict{Int, Int}()

for df_offer in groupby(filtered_donors, :DON_ID)
    df_list5 = filter(row -> row.ATT_TYPE == 5, df_offer)
    if !isempty(df_list5)
        append!(df_offers_recipients, df_list5)
        given_kidneys_registry[first(df_list5.DON_ID)] = count(skipmissing(df_list5.STATUS .== "TX"))
    end
end

offered_recipients = unique(df_offers_recipients.CAN_ID)

filtered_recipients = DataFrame(CAN_ID = Int[], CAN_BTH_DT = Date[], CAN_BLOOD = String[],
    CAN_A1 = Int[], CAN_A2 = Int[], CAN_B1 = Int[], CAN_B2 = Int[], CAN_DR1 = Int[], CAN_DR2 = Int[],
    CAN_LISTING_DT = DateTime[], CAN_DIAL_DT = Date[], 
    OUTCOME = String[], UPDATE_TM = DateTime[],
    CPRA = Int[], active_waiting_proportion = Float64[])

diff_score = DataFrame(CAN_ID = Int[], DON_ID = Int[], REGISTERED_SCORE = Float64[], COMPUTED_SCORE = Float64[])

for CAN_ID in offered_recipients

    df_offers_recipient = sort(filter(row-> row.CAN_ID == CAN_ID, df_offers_recipients), :DON_DEATH_TM, rev = true)

    # Récupération des informations sur le receveur
    df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)
    if ismissing(last(df_recipient.CAN_DIAL_DT))
        continue
    end
    df_recipient_first_registration = filter_first_waitlist_episode(df_recipient)
    p = KidneyAllocation.recipient_active_waiting_proportion(df_recipient_first_registration)

    # Récupération de la ligne correspondante à la première transplantation
    is_tx(outcome) = uppercase(strip(string(outcome))) == "TX"
    ind = findfirst(is_tx, skipmissing(df_offers_recipient.STATUS))

    # Si receveur non transplanté, on récupére la ligne de la dernière offre (1ère ligne dans le sorted dataframe)
    if isnothing(ind)
        ind = 1
    end

    # Date à laquelle on calcule le score
    target_date = df_offers_recipient.DON_DEATH_TM[ind]

    # Récupération du CPRA au moment du calcul du score
    df_cpra_recipient = sort(filter(row -> row.CAN_ID == CAN_ID, df_cpra), :UPDATE_TM, rev=true)
    latest_index = findfirst(t -> t ≤ target_date, df_cpra_recipient.UPDATE_TM)
    cpra = isnothing(latest_index) ? 0 : Int(df_cpra_recipient.CAN_CPRA[latest_index])

    # Transformation en Recipient et Donor structures pour le calcul du score
    r = first(df_recipient_first_registration)

    R = Recipient(r.CAN_BTH_DT, r.CAN_DIAL_DT, r.CAN_LISTING_DT, parse_abo(String(r.CAN_BLOOD)),
            r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
            cpra; expiration_date = nothing, active_waiting_proportion = 1.)

    D = KidneyAllocation.donor_from_row(df_offers_recipient[ind, :])

    computed_score = KidneyAllocation.score(D, R)
    registered_score = df_offers_recipient.DON_CAN_SCORE[ind]

    DON_ID = df_offers_recipient.DON_ID[ind]

    if !isapprox(registered_score, computed_score, atol=.1)
        @warn "Registered score $registered_score and computed score $computed_score mismatches between CAN_ID = $CAN_ID and DON_ID = $DON_ID."
        push!(diff_score, [CAN_ID, DON_ID, registered_score, computed_score])
    else
        push!(filtered_recipients,
                [r.CAN_ID, r.CAN_BTH_DT, r.CAN_BLOOD,
                r.CAN_A1, r.CAN_A2, r.CAN_B1, r.CAN_B2, r.CAN_DR1, r.CAN_DR2,
                r.CAN_LISTING_DT, r.CAN_DIAL_DT,
                r.OUTCOME, target_date,
                cpra, p])
    end
end


# Si le patient est transplanté, sa sortie s'effectue au moment du décès du donneur.
# Si le patient n'est pas transplanté, sa date de dernière mise à jour corespond à la date du dernier score calculé dans le fichier (ça se peut que le patient n'ait jamais eu l'opportunité d'accepter ou non).

# Patients pour lesquels il y a une différence entre le score calculé et le score enregistré dans le fichier.
diff_score

# Patients qui ont eu au moins un score de calculer et dont les scores calculé et enregistré correspondent.
filtered_recipients
