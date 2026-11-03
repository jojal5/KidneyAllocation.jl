using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Gadfly

using KidneyAllocation


## Chargement des données 

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)




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

