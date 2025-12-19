using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames

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


id = can_id[100]


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
        d = round(Int64, KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1]))
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


"""
    shift_recipient_timeline(recipient::Recipient, new_arrival::Date) -> Recipient

Shift the recipient's `birth`, `dialysis`, and `expiration_date` so that the
recipient's timeline is consistent with a new `arrival` date.

## Details

All dates are shifted by the same number of days: the difference between the
current `recipient.arrival` and the new `arrival`. This preserves age and
waiting-time durations relative to the (new) arrival date.
"""
function shift_recipient_timeline(recipient::Recipient, new_arrival::Date)::Recipient
    # Signed shift (in days) from old arrival to new arrival
    shift_days = KidneyAllocation.days_between(recipient.arrival, new_arrival)

    birth    = recipient.birth + Day(shift_days)
    dialysis = recipient.dialysis + Day(shift_days)

    expiration_date = recipient.expiration_date === nothing ? nothing :
                      recipient.expiration_date + Day(shift_days)

    return Recipient(birth,
                     dialysis,
                     new_arrival,
                     recipient.blood,
                     recipient.a1, recipient.a2,
                     recipient.b1, recipient.b2,
                     recipient.dr1, recipient.dr2,
                     recipient.cpra;
                     expiration_date = expiration_date)
end




shift_recipient_timeline(r, Date(2020,1,1))