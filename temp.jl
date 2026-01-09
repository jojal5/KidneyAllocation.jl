using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames


import KidneyAllocation: load_recipient, infer_recipient_expiration_date, parse_hla_int

filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"

df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

# Load the recipient file
df = KidneyAllocation.load_recipient(filepath)





# Estimate the recipient arrival rate
df = load_recipient(filepath)
df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)
λᵣ = length(unique(df2.CAN_ID)) / 6

# Construct the database of recipients














df = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/cleaned_candidates.csv", DataFrame)
df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
df.UPDATE_TM = Date.(df.UPDATE_TM)

df_cpra = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv", DataFrame)


filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

can_id = unique(df.CAN_ID)





id = 4 # Transplanté
id = 5 # Retiré


df_id = filter(row -> row.CAN_ID == id, df)
sort!(df_id, :UPDATE_TM, rev=true)

df_cpra_id = filter(row -> row.CAN_ID == id, df_cpra)
sort!(df_cpra_id, :UPDATE_TM, rev=true)

if df_id.OUTCOME[1] == "TX"
    exp_date = nothing # Si transplanté avec un donneur décédé, aucune date d'expiration et on néglige le temps inactif
else
    # Si non transplanté avec un donneur décédé, on prend la dernière date d'attente active pour calculer la date d'expiration
    ind = findfirst(df_id.OUTCOME .== "1")
    d = KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1])
    exp_date = df_id.UPDATE_TM[end] + Day(d)
end

birth = df_id.CAN_BTH_DT[1]
dialysis = df_id.CAN_DIAL_DT[1]
arrival = df_id.CAN_LISTING_DT[1]
blood = KidneyAllocation.parse_abo(df_id.CAN_BLOOD[1])
a1 = df_id.CAN_A1[1]
a2 = df_id.CAN_A2[1]
b1 = df_id.CAN_B1[1]
b2 = df_id.CAN_B2[1]
dr1 = df_id.CAN_DR1[1]
dr2 = df_id.CAN_DR2[1]

if isempty(df_cpra_id)
    cpra = 0 # Si le CPRA est manquant, on prend 0. 
else
    cpra = round(Int64, df_cpra_id.CAN_CPRA[end]) # S'il y a plusieurs CPRA, on prend le plus récent.
end

r = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra, expiration_date=exp_date)






recipients = Recipient[]

for id in can_id

    df_id = filter(row -> row.CAN_ID == id, df)
    sort!(df_id, :UPDATE_TM, rev=true)

    df_cpra_id = filter(row -> row.CAN_ID == id, df_cpra)
    sort!(df_cpra_id, :UPDATE_TM, rev=true)

    if (df_id.OUTCOME[1] == "TX") || (df_id.OUTCOME[1] == "1")
        exp_date = nothing # Si transplanté avec un donneur décédé, aucune date d'expiration et on néglige le temps inactif. Même chose si le patient est toujours actif en attente de transplantation.
    else
        # Si non transplanté avec un donneur décédé, on prend la dernière date d'attente active pour calculer la date d'expiration
        ind = findfirst(df_id.OUTCOME .== "1")
        if ind !== nothing
            d = KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1])
            exp_date = df_id.UPDATE_TM[end] + Day(d)
        else # Le patient a été inscrit mais n'a jamais été actif.
            exp_date = df_id.UPDATE_TM[end]
        end

    end

    birth = df_id.CAN_BTH_DT[1]
    dialysis = df_id.CAN_DIAL_DT[1]
    arrival = df_id.CAN_LISTING_DT[1]
    blood = KidneyAllocation.parse_abo(df_id.CAN_BLOOD[1])
    a1 = df_id.CAN_A1[1]
    a2 = df_id.CAN_A2[1]
    b1 = df_id.CAN_B1[1]
    b2 = df_id.CAN_B2[1]
    dr1 = df_id.CAN_DR1[1]
    dr2 = df_id.CAN_DR2[1]

    if isempty(df_cpra_id)
        cpra = 0 # Si le CPRA est manquant, on prend 0. 
    else
        cpra = round(Int64, df_cpra_id.CAN_CPRA[end]) # S'il y a plusieurs CPRA, on prend le plus récent.
    end

    r = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra, expiration_date=exp_date)
    push!(recipients, r)

end




import KidneyAllocation: evaluate_kdri, creatinine_mgdl

df = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/cleaned_donors.csv", DataFrame)

# Si on a une taille et un poids manquant, on remplace par les valeurs moyennes
df.WEIGHT[df.WEIGHT.≈0.] .= 80.
df.HEIGHT[df.HEIGHT.≈0.] .= 170.

# On retire les lignes dont l'age est manquant


don_id = unique(df.DON_ID)

id = 487

ind = findfirst(df.DON_ID .== id)

age = df.DON_AGE[ind]
height = df.HEIGHT[ind]
weight = df.WEIGHT[ind]
hypertension = df.HYPERTENSION[ind] == 1
diabetes = df.DIABETES[ind] == 1
cva = df.DEATH[ind] == 1
creatinine = creatinine_mgdl(df.CREATININE[ind])
dcd = df.DCD[ind] .== 1 # TODO À VÉRIFIER si c'est bien 1, sinon c'est 2 (Anastasiay a confirmé le code)

kdri = evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

arrival = Date(df.DON_DEATH_TM[ind])
age = df.DON_AGE[ind]
blood = KidneyAllocation.parse_abo(df.DON_BLOOD[ind])

a1 = df.DON_A1[ind]
a2 = df.DON_A2[ind]
b1 = df.DON_B1[ind]
b2 = df.DON_B2[ind]
dr1 = df.DON_DR1[ind]
dr2 = df.DON_DR2[ind]

d = Donor(arrival, age, blood, a1, a2, b1, b2, dr1, dr2, kdri)


donors = Donor[]

for id in don_id

    ind = findfirst(df.DON_ID .== id)

    age = df.DON_AGE[ind]
    height = df.HEIGHT[ind]
    weight = df.WEIGHT[ind]
    hypertension = df.HYPERTENSION[ind] == 1
    diabetes = df.DIABETES[ind] == 1
    cva = df.DEATH[ind] == 1
    creatinine = creatinine_mgdl(df.CREATININE[ind])
    dcd = df.DCD[ind] .== 1 # TODO À VÉRIFIER si c'est bien 1, sinon c'est 2 (Anastasiay a confirmé le code)

    kdri = evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

    arrival = Date(df.DON_DEATH_TM[ind])
    age = df.DON_AGE[ind]
    blood = KidneyAllocation.parse_abo(df.DON_BLOOD[ind])

    a1 = df.DON_A1[ind]
    a2 = df.DON_A2[ind]
    b1 = df.DON_B1[ind]
    b2 = df.DON_B2[ind]
    dr1 = df.DON_DR1[ind]
    dr2 = df.DON_DR2[ind]

    d = Donor(arrival, age, blood, a1, a2, b1, b2, dr1, dr2, kdri)

    push!(donors, d)

end



