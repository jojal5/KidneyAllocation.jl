using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames

filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"

df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])


# Replacing the value 24L and 24Low with 24 
parse_hla_int(x) = x === missing ? missing :
                   (v = tryparse(Int, replace(String(x), r"L$" => "")); v === nothing ? missing : v)

df.CAN_A2 = parse_hla_int.(df.CAN_A2)  # broadcast

# Keeping only the recipients for kidney transplant
filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

# Removing the recipient where the dialysis time or the listing time is missing. It happens for recipient that received a kidney from a living donor.
dropmissing!(df, [:CAN_LISTING_DT, :CAN_DIAL_DT])

df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
df.UPDATE_TM = Date.(df.UPDATE_TM)


"""
    load_recipient(filepath::String)

Load the CSV file containing recipient information and return a cleaned DataFrame

##Details

Here are the cleaning steps:
- Replacing the value 24L and 24Low with 24 in CAN_A2
- Keeping only the recipient for kidney transplant, not multi-organs
- Removing the recipient where the dialysis time or the listing time is missing. It happens for recipient that received a kidney from a living donor.
- Converting DateTime format to Date
"""
function load_recipient(filepath::String)

    df = CSV.read(filepath, DataFrame, missingstring=["-", "", "NULL"])

    # Replacing the value 24L and 24Low with 24 for instance
    df.CAN_A2 = parse_hla_int.(df.CAN_A2)

    # Keeping only the recipients for kidney transplant
    filter!(row -> row.OUTCOME ∈ ("TX", "1", "0", "X", "Dcd", "Tx Vivant"), df)

    # Removing the recipient where the dialysis time or the listing time is missing. It happens for recipients that received a kidney from a living donor.
    dropmissing!(df, [:CAN_LISTING_DT, :CAN_DIAL_DT])

    # Converting DateTime format to Date
    df.CAN_LISTING_DT = Date.(df.CAN_LISTING_DT)
    df.CAN_DIAL_DT = Date.(df.CAN_DIAL_DT)
    df.UPDATE_TM = passmissing(Date).(df.UPDATE_TM)

    return df

end

# Arrival rate of eligible recipients
df = load_recipient(filepath)

df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)

λᵣ = length(unique(df2.CAN_ID)) / 6



df = load_recipient(filepath)
df_cpra = CSV.read(cpra_filepath, DataFrame)

select!(df, Not(:NB_PAST_TRANS))
dropmissing!(df)

id = 1 # Transplanté
# id = 5 # Retiré


df_id = filter(row -> row.CAN_ID == id, df)
sort!(df_id, :UPDATE_TM, rev=true)

df_cpra_id = filter(row -> row.CAN_ID == id, df_cpra)
sort!(df_cpra_id, :UPDATE_TM, rev=true)





function construct_recipient_list(recipient_filepath::String, cpra_filepath::String)

    df = load_recipient(recipient_filepath)
    df_cpra = CSV.read(cpra_filepath, DataFrame)

    select!(df, Not(:NB_PAST_TRANS))
    dropmissing!(df)

    recipients = Recipient[]

    for id in unique(df.CAN_ID)

        # id = 1 # Transplanté
        # id = 5 # Retiré


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

    return recipients

end


r = construct_recipient_list(filepath, cpra_filepath)






cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

df_cpra = CSV.read(cpra_filepath, DataFrame)






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









normalize_outcome(s) = uppercase(strip(String(s)))

normalize_outcome("0")


