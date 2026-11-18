using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, Distributions, Statistics

using Cairo, Fontconfig, Gadfly
using ExtremePlots

using KidneyAllocation

fig_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Figures"


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

# On ne retient que les patients majeurs
filter!(row -> KidneyAllocation.years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df_recipients)

# On ne retient que la première inscription des candidats
# df_donors = unique( sort(df_donors, :DON_DEATH_TM), :CAN_ID)

donor_registry, given_kidneys = KidneyAllocation.build_donor_registry(df_donors)
recipient_registry = KidneyAllocation.build_recipient_registry(df_recipients)

## Candidates on the waiting list

t = Date(2012,1,1):Month(1):Date(2024,1,1)
n = Vector{Int64}(undef, length(t))

for (i,tᵢ) in enumerate(t)
    ids = KidneyAllocation.retrieve_observed_waiting_list(df_recipients, tᵢ)
    n[i] = length(ids)
end

df = DataFrame(Date = t, Candidates = n)

fig = plot(df, x=:Date, y=:Candidates,
    Coord.cartesian(xmax=Date(2022,1,1)),
    Guide.ylabel("Waiting candidates"),
    Theme(default_color = "gray30"))

    
filepath = joinpath(fig_folder, "waiting_candidates.pdf")
draw(PDF(filepath), fig)

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

filter!(row -> row.Date ≥ Date(2012,4, 1) && row.Date < Date(2022,8,1), df)

fig = plot(df, x=:Date, y=:Donors, Theme(default_color = "gray30"),
    Guide.ylabel("Donor monthly arrivals"))

filepath = joinpath(fig_folder, "monthly_donors.pdf")
draw(PDF(filepath), fig)



λ̂ = mean(df.Donors)
pd = Poisson(λ̂)
fig = qqplot(pd, df.Donors)
filepath = joinpath(fig_folder, "fitted_Poisson_donors.pdf")
draw(PDF(filepath), fig)


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

filter!(row -> row.Date ≥ Date(2012,4, 1) && row.Date < Date(2022,8,1), df)

fig = plot(df, x=:Date, y=:Candidates, Theme(default_color = "gray30"),
    Guide.ylabel("Candiate monthly arrivals"))

filepath = joinpath(fig_folder, "monthly_candidates.pdf")
draw(PDF(filepath), fig)


λ̂ = mean(df.Candidates)
pd = Poisson(λ̂)
fig = qqplot(pd, df.Candidates)
filepath = joinpath(fig_folder, "fitted_Poisson_recipients.pdf")
draw(PDF(filepath), fig)



## Retrieve transplanted recipients for prediction

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)

TX_ids = Int.(df_TX.CAN_ID)

## Time to tx
# Computed in elapsed days since dialysis

time_to_transplant = Dict{Int, Int}()

for row in eachrow(df_donors)
    if haskey(recipient_registry, row.CAN_ID)
        recipient = recipient_registry[row.CAN_ID]
        elapsed_days = Dates.days(Dates.Date(row.DON_DEATH_TM) - recipient.dialysis)
        time_to_transplant[row.CAN_ID] = elapsed_days
    end
end

## Time to TX as function of blood group

df = DataFrame(CAN_ID = Int[], ABO = String[], time = Int[], years = Float64[])

n_notfound = 0

for id in TX_ids
    if haskey(time_to_transplant, id)
        blood = recipient_registry[id].blood
        waittime = time_to_transplant[id]
        push!(df, [id, string(blood), waittime, waittime/365.25])
    else
        n_notfound +=1
    end
end

filter!(row -> row.years ≤ 10, df)

# To have the order A, B, AB and O
temp = df[6,:]
df[20, :]
df[6,:] = df[20,:]
df[20,:] = temp

fig = plot(df, x=:ABO, y=:years, Geom.boxplot,
    Guide.xlabel("ABO Groups"),
    Guide.ylabel("Waiting Years"),
    Guide.yticks(ticks=0:10),
    Theme(default_color = "gray30"))

filepath = joinpath(fig_folder, "obs_waittime_before_tx.pdf")
draw(PDF(filepath), fig)



# Waiting time should not be longer than 1461 between inscription and transplant.
# But the waiting time is computed from the dialysis date.

filter(row -> row.CAN_ID == 15472, df_recipients)
filter(row -> row.CAN_ID == 9985, df_recipients)

filter(row -> row.CAN_ID == 9985, df_donors)

recipient_registry[9985]
time_to_transplant[9985]



## Time to first offer
# Computed in elapsed days since dialysis

time_to_first_offer = Dict{Int, Int}()

df_offers = KidneyAllocation.load_donor(donor_filepath)
filter!(row -> row.ATT_TYPE == 5, df_offers)

df_offers_filtered = similar(df_offers, 0)

for df_offer in groupby(df_offers, :DON_ID)

    filtered_df = KidneyAllocation.filter_offers_until_last_acceptance(df_offer)
    append!(df_offers_filtered, filtered_df)

end

for g in groupby(df_offers_filtered, :CAN_ID)

    CAN_ID = first(g.CAN_ID)

    if haskey(recipient_registry, CAN_ID)
        recipient = recipient_registry[CAN_ID]
        elapsed_days = Dates.days(Dates.Date(minimum(g.DON_DEATH_TM)) - recipient.dialysis)
        time_to_first_offer[CAN_ID] = elapsed_days
    end

end

