using Pkg
Pkg.activate(".")

using KidneyAllocation
using Test

using Dates, CSV, DataFrames

df = CSV.read("/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/cleaned_candidates.csv", DataFrame)

filter!(row -> row.OUTCOME ∈ ("TX","1","0","X","Dcd","Tx Vivant"), df)

can_id = unique(df.CAN_ID)

recipients = Recipient[]

# for id in [1]
#     df_id = filter(row -> row.CAN_ID ==id )
#     sort!(df_id, :UPDATE_TM, rev=true)

#     if df.OUTCOME[end] != "TX"
#         exp_date = KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[1])

# end

id = 5
df_id = filter(row -> row.CAN_ID ==id, df )
sort!(df_id, :UPDATE_TM, rev=true)

ind = findfirst(df_id.OUTCOME .== "1")
d = round(Int64, KidneyAllocation.days_between(df_id.UPDATE_TM[end], df_id.UPDATE_TM[ind-1]))
expiration_date = df_id.UPDATE_TM[end] + Day(d)

birth = df_id.CAN_BTH_DT[1]
dialysis = df_id.CAN_DIAL_DT[1]
arrival = df_id.CAN_LISTING_DT[1]
# blood = ABOGroup(df_id.CAN_BLOOD[1])
blood = O
a1 = df_id.CAN_A1[1]
a2 = df_id.CAN_A2[1]
b1 = df_id.CAN_B1[1]
b2 = df_id.CAN_B2[1]
dr1 = df_id.CAN_DR1[1]
dr2 = df_id.CAN_DR2[1]
cpra = 10



# r = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra, expiration_date = expiration_date)
r = Recipient(birth, dialysis, arrival, blood, a1, a2, b1, b2, dr1, dr2, cpra, expiration_date = nothing)

birth::DateTime
    dialysis::DateTime
    arrival::DateTime
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    CPRA::Int64
    expiration_date::Union{DateTime,Nothing}