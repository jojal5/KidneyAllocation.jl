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

@testset "shift_recipient_timeline" begin
    import KidneyAllocation.days_between
    
    birth    = Date(1980, 1, 1)
    dialysis = Date(2010, 1, 1)
    arrival  = Date(2020, 1, 1)
    expiry   = Date(2025, 1, 1)

    r = Recipient(
        birth,
        dialysis,
        arrival,
        A,
        24, 26,
        44, 51,
        1, 4,
        80;
        expiration_date = expiry
    )

    # --- Forward shift ---
    new_arrival = Date(2022, 1, 1)
    r2 = shift_recipient_timeline(r, new_arrival)

    shift_days = days_between(arrival, new_arrival)

    @test r2.arrival == new_arrival
    @test r2.birth == birth + Day(shift_days)
    @test r2.dialysis == dialysis + Day(shift_days)
    @test r2.expiration_date == expiry + Day(shift_days)

    # Preserve non-date fields
    @test r2.blood == r.blood
    @test r2.a1 == r.a1
    @test r2.b2 == r.b2
    @test r2.dr1 == r.dr1
    @test r2.cpra == r.cpra

    # --- Backward shift ---
    earlier_arrival = Date(2018, 1, 1)
    r3 = shift_recipient_timeline(r, earlier_arrival)

    shift_days_back = days_between(arrival, earlier_arrival)

    @test r3.arrival == earlier_arrival
    @test r3.birth == birth + Day(shift_days_back)
    @test r3.dialysis == dialysis + Day(shift_days_back)
    @test r3.expiration_date == expiry + Day(shift_days_back)

    # --- No expiration date ---
    r_noexp = Recipient(
        birth,
        dialysis,
        arrival,
        O,
        24, 26,
        44, 51,
        1, 4,
        10
    )

    r4 = shift_recipient_timeline(r_noexp, new_arrival)

    @test r4.expiration_date === nothing
    @test r4.birth == birth + Day(shift_days)
    @test r4.dialysis == dialysis + Day(shift_days)
end





using Random, Dates, Test

"""
    sample_days(d1::Date, d2::Date, n::Integer) -> Vector{Date}

Sample `n` dates uniformly (with replacement) between `d1` and `d2` (inclusive).
"""
function sample_days(d1::Date, d2::Date, n::Integer)::Vector{Date}
    dmin, dmax = min(d1, d2), max(d1, d2)

    ndays = Dates.value(dmax - dmin) + 1
    offsets = sort(rand(0:ndays-1, n))

    return dmin .+ Day.(offsets)
end

d1 = Date(2000,1,1)
d2 = Date(2020,12,31)

@time sample_days(d1, d2, 100)


@testset "sample_days" begin
    # --- Basic properties ---
    d1 = Date(2024, 1, 1)
    d2 = Date(2024, 1, 10)

    s = sample_days(d1, d2, 10)

    @test length(s) == 10
    @test all(d -> d1 <= d <= d2, s)
end
