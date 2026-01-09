using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames


import KidneyAllocation: build_recipient_registry, load_recipient

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

build_recipient_registry(recipient_filepath, cpra_filepath)


# Estimate the recipient arrival rate
df = load_recipient(recipient_filepath)
df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)
λᵣ = length(unique(df2.CAN_ID)) / 6






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



