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


df = CSV.read("test/data/unfiltered_outcomes.csv", DataFrame)
G = groupby(df, :CAN_ID)

    @testset "filter_first_waitlist_episode" begin
        # import KidneyAllocation.filter_first_waitlist_episode

        @test nrow(filter_first_waitlist_episode(G[1])) == nrow(G[1])
        @test nrow(filter_first_waitlist_episode(G[2])) == 2
        @test nrow(filter_first_waitlist_episode(G[3])) == nrow(G[3])
        @test nrow(filter_first_waitlist_episode(G[4])) == 2
        @test nrow(filter_first_waitlist_episode(G[5])) == nrow(G[5])
        @test nrow(filter_first_waitlist_episode(G[6])) == nrow(G[6])
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


## Récupération des candidats dont au moins un score a été calculé pour la liste 5
# Reconstruction des receveurs par le fichier des donneurs

is_tx(outcome::Union{Missing, AbstractString}) = uppercase(strip(string(outcome))) == "TX"

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
# filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)
G = groupby(filtered_donors, :DON_ID)

given_kidneys_registry = Dict{Int, Int}()
for df_offer in groupby(filtered_donors, :DON_ID)
    ind = findall(is_tx, df_offer.STATUS)
    n = count(df_offer.ATT_TYPE[ind] .== 5)
    given_kidneys_registry[first(df_offer.DON_ID)] = n
end

# filter(row -> row.DON_ID == 6073, filtered_donors)


df_offers_recipients = similar(filtered_donors, 0)
offered_recipients = unique(df_donors.CAN_ID)

filtered_recipients = DataFrame(CAN_ID = Int[], CAN_BTH_DT = Date[], CAN_BLOOD = String[],
    CAN_A1 = Int[], CAN_A2 = Int[], CAN_B1 = Int[], CAN_B2 = Int[], CAN_DR1 = Int[], CAN_DR2 = Int[],
    CAN_LISTING_DT = DateTime[], CAN_DIAL_DT = Date[], 
    OUTCOME = String[], UPDATE_TM = DateTime[],
    CPRA = Int[], active_waiting_proportion = Float64[])

allowmissing!(filtered_recipients)


for CAN_ID in offered_recipients
# CAN_ID = 5

    df_offers_recipient = sort(filter(row-> row.CAN_ID == CAN_ID, filtered_donors), :DON_DEATH_TM)

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
    df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)
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

# Si le patient est transplanté, sa sortie s'effectue au moment du décès du donneur.
# Si le patient n'est pas transplanté, sa date de dernière mise à jour correspond à la date du dernier score calculé dans le fichier (ça se peut que le patient n'ait jamais eu l'opportunité d'accepter ou non).


## Check registered score against computed score
# A difference in score indicates an error in the data


df_scores = DataFrame(CAN_ID = Int[], DON_ID = Int[], REGISTERED_SCORE = Float64[], COMPUTED_SCORE = Float64[])
allowmissing!(df_scores, :COMPUTED_SCORE)

for r in eachrow(filtered_recipients)

    df_offers = sort(filter(row-> row.CAN_ID == r.CAN_ID, df_donors), :DON_DEATH_TM)

    ind = findfirst(is_tx, df_offers.STATUS)
    if isnothing(ind)
        ind = findlast(df_offers.DON_DEATH_TM .≤ r.UPDATE_TM)
    end

    registered_score = df_offers.DON_CAN_SCORE[ind]

    try
        R = KidneyAllocation.recipient_from_row(r)
        D = KidneyAllocation.donor_from_row(df_offers[ind, :])
    catch
        push!(df_scores, [r.CAN_ID, df_offers.DON_ID[ind], registered_score, missing])
        continue # If there's missing values, score cannot be computed
    end

    computed_score = KidneyAllocation.score(D, R)
    

    push!(df_scores, [r.CAN_ID, df_offers.DON_ID[ind], registered_score, computed_score])

end

score_diff = abs.(df_scores.COMPUTED_SCORE - df_scores.REGISTERED_SCORE)
count(skipmissing(score_diff .> .1))
count(skipmissing(score_diff .≤ .1))
count(ismissing.(score_diff))


score_matching_flag = Vector{Union{Bool, Missing}}(undef, length(score_diff))

for i in eachindex(score_diff)
    if ismissing(score_diff[i])
        score_matching_flag[i] = missing
    else
        score_matching_flag[i] = score_diff[i] ≤ .1
    end
end

filtered_recipients.SCORE_MATCHING_FLAG = score_matching_flag


## Récupération des donneurs de la liste 5

df_donors_list5 = similar(filtered_donors, 0)

for df_offer in groupby(filtered_donors, :DON_ID)

    ind = findall(is_tx, df_offer.STATUS)
    n = count(df_offer.ATT_TYPE[ind] .== 5)
    
    if n > 0
        append!(df_donors_list5, df_offer[ind, :])
    end

end

df_donors_list5



## Fonction de chargement des données

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

    is_tx(outcome::Union{Missing, AbstractString}) = uppercase(strip(string(outcome))) == "TX"

    filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
    filter!(row -> row.DON_ID ≠ 0, filtered_donors)

    df_offers_recipients = similar(filtered_donors, 0)
offered_recipients = unique(df_donors.CAN_ID)

filtered_recipients = DataFrame(CAN_ID = Int[], CAN_BTH_DT = Date[], CAN_BLOOD = String[],
    CAN_A1 = Int[], CAN_A2 = Int[], CAN_B1 = Int[], CAN_B2 = Int[], CAN_DR1 = Int[], CAN_DR2 = Int[],
    CAN_LISTING_DT = DateTime[], CAN_DIAL_DT = Date[], 
    OUTCOME = String[], UPDATE_TM = DateTime[],
    CPRA = Int[], active_waiting_proportion = Float64[])

allowmissing!(filtered_recipients)


    for CAN_ID in offered_recipients

        df_offers_recipient = sort(filter(row-> row.CAN_ID == CAN_ID, filtered_donors), :DON_DEATH_TM)

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
        df_recipient = filter(row-> row.CAN_ID == CAN_ID, df_recipients)
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

    df_scores = DataFrame(CAN_ID = Int[], DON_ID = Int[], REGISTERED_SCORE = Float64[], COMPUTED_SCORE = Float64[])
    allowmissing!(df_scores, :COMPUTED_SCORE)

    for r in eachrow(filtered_recipients)

        df_offers = sort(filter(row-> row.CAN_ID == r.CAN_ID, df_donors), :DON_DEATH_TM)

        ind = findfirst(is_tx, df_offers.STATUS)
        if isnothing(ind)
            ind = findlast(df_offers.DON_DEATH_TM .≤ r.UPDATE_TM)
        end

        registered_score = df_offers.DON_CAN_SCORE[ind]

        try
            R = KidneyAllocation.recipient_from_row(r)
            D = KidneyAllocation.donor_from_row(df_offers[ind, :])
        catch
            push!(df_scores, [r.CAN_ID, df_offers.DON_ID[ind], registered_score, missing])
            continue # If there's missing values, score cannot be computed
        end

        computed_score = KidneyAllocation.score(D, R)
        

        push!(df_scores, [r.CAN_ID, df_offers.DON_ID[ind], registered_score, computed_score])

    end

    score_diff = abs.(df_scores.COMPUTED_SCORE - df_scores.REGISTERED_SCORE)
    # count(skipmissing(score_diff .> .1))
    # count(skipmissing(score_diff .≤ .1))
    # count(ismissing.(score_diff))

    score_matching_flag = Vector{Union{Bool, Missing}}(undef, length(score_diff))

    for i in eachindex(score_diff)
        if ismissing(score_diff[i])
            score_matching_flag[i] = missing
        else
            score_matching_flag[i] = score_diff[i] ≤ .1
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

df_recipients, df_donors = load_data(recipient_filepath, cpra_filepath, donor_filepath)





## Récupération des candidats dont au moins un score a été calculé pour la liste 5 entre 2016 et 2020
# Reconstruction des recveurs par le fichier des donneurs

filtered_donors = dropmissing(df_donors, :DON_DEATH_TM)
filter!(row -> row.DON_ID ≠ 0, filtered_donors)
# filter!(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM ≤ Date(2019,12,31), filtered_donors)
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


    d = df_offers_recipient[ind, :]
    
    # D = KidneyAllocation.donor_from_row(d)

    D = Donor(d.DON_DEATH_TM, d.DON_AGE, parse_abo(String(d.DON_BLOOD)), d.DON_A1, d.DON_A2, d.DON_B1, d.DON_B2, d.DON_DR1, d.DON_DR2, 10.)

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






## Récupération des décisions

df_decision = similar(filtered_donors, 0)

for df_offer in groupby(filtered_donors, :DON_ID)

# df_offer = filter(row -> row.DON_ID == 616, filtered_donors)

    println(df_offer.DON_ID[1])

    ind = findall(is_tx, df_offer.STATUS)

    n = clamp(2 - count(df_offer.ATT_TYPE[ind] .!= 5) - count(df_offer.ATT_TYPE .== 0), 0, 2)
    if n == 0 # No kidney offerred on list 5
        continue
    end 

    df_offer_list5_sorted = sort( filter(row -> row.ATT_TYPE==5, df_offer) ,:DON_CAN_SCORE, rev = true)


    ind = findall(is_tx, df_offer_list5_sorted.STATUS)
    if isnothing(ind) || length(ind) < n   # At least a kidney is refused by all recipients
        append!(df_decision, df_offer_list5_sorted)
    else
        append!(df_decision, df_offer_list5_sorted[1:ind[n], :])
    end

end


## Il faut transférer en varibles ... temps d'attente, mismatch, kdri ...




n_train = round(Int, nrow(df_decision) * 2/3)
rng = MersenneTwister(12345)
train_ind = sort(shuffle(rng, 1:nrow(df_decision))[1:n_train])
validation_ind = sort(setdiff(1:nrow(df_decision), train_ind))

df_train = df_decision[train_ind, :]
df_validation = df_decision[validation_ind, :]

