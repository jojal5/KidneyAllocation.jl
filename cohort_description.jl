using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Gadfly, Statistics

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

donor_registry, given_kidneys_registry = KidneyAllocation.build_donor_registry(df_donors) 
recipient_registry, time_to_first_offer_registry, time_to_transplant_registry = KidneyAllocation.build_recipient_registry(df_recipients)


## Candidates on the waiting list

t = Date(2012,1,1):Month(1):Date(2024,1,1)
n = Vector{Int64}(undef, length(t))

for (i,tᵢ) in enumerate(t)
    ids = KidneyAllocation.retrieve_observed_waiting_list(df_recipients, tᵢ)
    n[i] = length(ids)
end

df = DataFrame(Date = t, Candidates = n)

plot(df, x=:Date, y=:Candidates)


## Monthly donor arrivals

df_donor_arrival = DataFrame(DON_ID = Int64[], DON_DEATH_TM = DateTime[])

for g in groupby(df_donors, :DON_ID)
    push!(df_donor_arrival, [first(g.DON_ID), first(g.DON_DEATH_TM)])
end
# groupby(df_donors, :DON_DEATH_TM => yearmonth)

df_donor_arrival.year = year.(df_donor_arrival.DON_DEATH_TM)
df_donor_arrival.month = month.(df_donor_arrival.DON_DEATH_TM)

df = combine(groupby(df_donor_arrival, [:year, :month]), :DON_DEATH_TM => length => :Donors)

d = [Date(r.year, r.month, 1) for r in eachrow(df) ]
df.Date = d

plot(df, x=:Date, y=:Donors)


## Monthly candidate arrivals

df_recipient_arrival = DataFrame(CAN_ID = Int64[], CAN_LISTING_DT = DateTime[])

for g in groupby(df_recipients, :CAN_ID)
    push!(df_recipient_arrival, [first(g.CAN_ID), first(g.CAN_LISTING_DT)])
end

df_recipient_arrival.year = year.(df_recipient_arrival.CAN_LISTING_DT)
df_recipient_arrival.month = month.(df_recipient_arrival.CAN_LISTING_DT)

df = combine(groupby(df_recipient_arrival, [:year, :month]), :CAN_LISTING_DT => length => :Candidates)

d = [Date(r.year, r.month, 1) for r in eachrow(df) ]
df.Date = d

plot(df, x=:Date, y=:Candidates)



## Retrieve transplanted recipients for prediction

TX_ids = Int64[]

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

TX_ids = df_TX.CAN_ID

