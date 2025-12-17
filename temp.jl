using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames

df = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/cleaned_candidates.csv", DataFrame)
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
    d = round(Int64, KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1]))
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

    if df_id.OUTCOME[1] == "TX"
        exp_date = nothing # Si transplanté avec un donneur décédé, aucune date d'expiration et on néglige le temps inactif
    else
        # Si non transplanté avec un donneur décédé, on prend la dernière date d'attente active pour calculer la date d'expiration
        ind = findfirst(df_id.OUTCOME .== "1")
        d = round(Int64, KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1]))
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
    push!(recipients, r)

end






HLA_A = DataFrame(VALID_HLA_A = sort(unique(union(df.CAN_A1,df.CAN_A2))))
CSV.write("valid_HLA-A.csv", HLA_A)

HLA_B = DataFrame(VALID_HLA_B = sort(unique(union(df.CAN_B1,df.CAN_B2))))
CSV.write("valid_HLA-B.csv", HLA_B)

HLA_DR = DataFrame(VALID_HLA_DR = sort(unique(union(df.CAN_DR1,df.CAN_DR2))))
CSV.write("valid_HLA-DR.csv", HLA_DR)