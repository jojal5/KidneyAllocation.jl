using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Statistics, Test

using Cairo, Fontconfig, Gadfly

using KidneyAllocation

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

## Load data

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

# On ne retient que les patients majeurs
filter!(row -> KidneyAllocation.years_between(row.CAN_BTH_DT, row.CAN_LISTING_DT) > 17, df_recipients)

# On ne retient que la première inscription des candidats
# df_donors = unique( sort(df_donors, :DON_DEATH_TM), :CAN_ID)
# Ça fait qu'on ignore les patients qui reviennent et les arrivées des reins qui leurs sont attribués

donor_registry, given_kidneys_registry = KidneyAllocation.build_donor_registry(df_donors)
recipient_registry = KidneyAllocation.build_recipient_registry(df_recipients)


## Transplanted recipient

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)
TX_ids = Int.(df_TX.CAN_ID)

elapsed_days_since_arrival = Dict{Int,Int}()
elapsed_days_since_dialysis = Dict{Int,Int}()

for r in eachrow(df_TX)
    elapsed_days_since_arrival[r.CAN_ID] = Dates.days(r.UPDATE_TM - r.CAN_LISTING_DT)
    elapsed_days_since_dialysis[r.CAN_ID] = Dates.days(Date(r.UPDATE_TM) - r.CAN_DIAL_DT)
end

## Check calibration for all recipients - Time to transplant

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

inbox = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = elapsed_days_since_dialysis[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

    inbox[i] = t1 ≤ obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)

mean(inbox)

## Transplanted recipient - A

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)
filter!(row->row.CAN_BLOOD == "A", df_TX)
TX_ids = Int.(df_TX.CAN_ID)


## Check calibration for all recipients - A

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

inbox = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = elapsed_days_since_dialysis[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

    inbox[i] = t1 ≤ obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)

mean(inbox)

## Transplanted recipient - B

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)
filter!(row->row.CAN_BLOOD == "B", df_TX)
TX_ids = Int.(df_TX.CAN_ID)


## Check calibration for all recipients - B

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

inbox = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = elapsed_days_since_dialysis[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

    inbox[i] = t1 ≤ obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)

mean(inbox)

## Transplanted recipient - AB

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)
filter!(row->row.CAN_BLOOD == "AB", df_TX)
TX_ids = Int.(df_TX.CAN_ID)


## Check calibration for all recipients - AB

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

inbox = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = elapsed_days_since_dialysis[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

    inbox[i] = t1 ≤ obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)

mean(inbox)

## Transplanted recipient - O

df_TX = filter(row-> row.OUTCOME == "TX", df_recipients)
dropmissing!(df_TX, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_TX)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_TX)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_TX)
filter!(row->row.CPRA < 80, df_TX)
filter!(row->row.CAN_BLOOD == "O", df_TX)
TX_ids = Int.(df_TX.CAN_ID)


## Check calibration for all recipients - O

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

q1 = falses(length(TX_ids))
q2 = falses(length(TX_ids))
q3 = falses(length(TX_ids))

inbox = falses(length(TX_ids))

for (i, id) in enumerate(TX_ids)

    # obs_waiting_time_before_first_offer = time_to_first_offer_registry[id]
    obs_waiting_time_before_first_offer = elapsed_days_since_dialysis[id]

    filepath = joinpath(sim_folder, "$id.csv")
    res = CSV.read(filepath, DataFrame)

    t1 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .25)
    t2 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .5)
    t3 = KidneyAllocation.first_qualifying_offer_quantile_given_offer(res, .75)

    q1[i] = !isnothing(t1) && obs_waiting_time_before_first_offer ≤ t1
    q2[i] = !isnothing(t2) && obs_waiting_time_before_first_offer ≤ t2
    q3[i] = !isnothing(t3) && obs_waiting_time_before_first_offer ≤ t3

    inbox[i] = t1 ≤ obs_waiting_time_before_first_offer ≤ t3

end

p₁ = count(q1)/length(q1)
p₂ = count(q2)/length(q2)
p₃ = count(q3)/length(q3)

mean(inbox)

## Affichage des résultats

# id = 6740
#   7847 # 530 jours A CPRA=0 p=1 61 years old
#   8275 # 4728 jours
#   6740 # 1676 jours B CPRA=0 p=1 56 years old
#   8081 # 755 jours O CPRA = 0 p=.99 45 years old
#   7958 # 624 jours
#   8121 # 355 jours AB CPRA=0 p=1 25 years old
#   8056 # 1129 jours O CPRA = 0 p=.86 73 years old
#   7325 # 1543 jours O CPRA = 19 p=.64 49 years old
#  7897
#  8595
#  4512
#  8159
#  7633
#  7800
#  8586
#  8550
#  7959
#  8538

filter(row->row.CAN_ID == 7847, df_donors)
Dates.days(Date(2016,10,15) - r.arrival)
Dates.days(Date(2016,10,15) - r.dialysis)

donor_registry[2692]

id = 7847
r = recipient_registry[id]
filepath = joinpath(sim_folder, "$id.csv")
res = CSV.read(filepath, DataFrame)
df_A = combine(groupby(res, :simulation_id), :elapsed_days_since_dialysis => minimum => :elapsed_days_since_dialysis)
rename!(df_A, :elapsed_days_since_dialysis => :A)

id = 6740
recipient_registry[id]
filepath = joinpath(sim_folder, "$id.csv")
res = CSV.read(filepath, DataFrame)
df_B = combine(groupby(res, :simulation_id), :elapsed_days_since_dialysis => minimum => :elapsed_days_since_dialysis)
rename!(df_B, :elapsed_days_since_dialysis => :B)

id = 8121
recipient_registry[id]
filepath = joinpath(sim_folder, "$id.csv")
res = CSV.read(filepath, DataFrame)
df_AB = combine(groupby(res, :simulation_id), :elapsed_days_since_dialysis => minimum => :elapsed_days_since_dialysis)
rename!(df_AB, :elapsed_days_since_dialysis => :AB)

id = 7325
recipient_registry[id]
filepath = joinpath(sim_folder, "$id.csv")
res = CSV.read(filepath, DataFrame)
df_O = combine(groupby(res, :simulation_id), :elapsed_days_since_dialysis => minimum => :elapsed_days_since_dialysis)
rename!(df_O, :elapsed_days_since_dialysis => :O)

df = innerjoin(df_A, df_B, df_AB, df_O, on=:simulation_id)

longdf = stack(df, Not(:simulation_id))
rename!(longdf, :variable => :Patient, :value => :Days)

longdf.Years = longdf.Days/365.25

simulated_times = layer(longdf, x=:Patient, y=:Years, Geom.boxplot,
    Theme(default_color = "gray30"))

observed_waiting = DataFrame(Patient = ["A", "B", "AB", "O"], Days = [530, 1676, 355, 1543])

observed_waiting.Years = observed_waiting.Days / 365.25

observed_times =layer(observed_waiting, x=:Patient, y=:Years, Geom.point,
    Theme(default_color = "black", point_shapes = [Shape.cross]))

fig = plot(observed_times, simulated_times, 
    Guide.ylabel("Simulated waiting time (years)"))

fig_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Figures"
filepath = joinpath(fig_folder, "waiting_times.pdf")

draw(PDF(filepath), fig)


## Permanently removed patients

X_ids = Int64[]

removed_status = ["X", "TX VIVANT", "DCD"]

df_X = filter(row-> uppercase(strip(row.OUTCOME)) ∈ removed_status, df_recipients)
dropmissing!(df_X, :SCORE_MATCHING_FLAG)
filter!(row->row.SCORE_MATCHING_FLAG, df_X)
filter!(row->row.CAN_LISTING_DT ≥ Date(2016,1,1), df_X)
filter!(row->row.UPDATE_TM < Date(2020,1,1), df_X)
filter!(row->row.CPRA < 80, df_X)
filter!(row -> row.CAN_ID ∈ keys(recipient_registry), df_X)

X_ids = df_X.CAN_ID

removal_since_dialysis = Dict{Int,Int}()

for r in eachrow(df_X)
    removal_since_dialysis[r.CAN_ID] =  Dates.days(Date(r.UPDATE_TM) - r.CAN_DIAL_DT)
end


## Check calibration for permanently removed patients

sim_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Simulations"

p = Vector{Float64}(undef, length(X_ids))
target = Vector{Int}(undef, length(X_ids))

df_offers = KidneyAllocation.load_donor(donor_filepath)

for (i, id) in enumerate(X_ids)

    df = filter(row->row.CAN_ID == id, df_offers)
    sort!(df, :DON_DEATH_TM)
    target[i] = Dates.days(Date(df.DON_DEATH_TM[1]) - recipient_registry[id].dialysis)
    
    filepath = joinpath(sim_folder, "$id.csv")

    res = CSV.read(filepath, DataFrame)

    # target = removal_since_dialysis[id]

    p[i] = KidneyAllocation.probability_qualifying_offer_by(res, target[i]; n_simulations=1000)

end
p

count(p .< .5)

p

fig = plot(y=p, Geom.point,
    Guide.xlabel("Recipient"),
    Guide.ylabel("Estimated probability of an offer"),
    Guide.xticks(ticks=[1;86]),
    Guide.yticks(ticks=0:.25:1),
    Theme(default_color = "gray30"))

fig_folder = "/Users/jalbert/Dropbox/Files/Supervision/encours/AnastasiyaOlekBasanets/Figures"
filepath = joinpath(fig_folder, "probability_offer_before_removal.pdf")
draw(PDF(filepath), fig)




## Check calibration for probability of first offer before target time

df = [df_X; df_TX]

target_date = Dict{Int, Int}()

for r in eachrow(df)
    target_date[r.CAN_ID] = Dates.days(Date(r.UPDATE_TM) - r.CAN_DIAL_DT)
end

gt = df.OUTCOME .== "TX"

p = Vector{Float64}(undef, nrow(df))

for (i, r) in enumerate(eachrow(df))

    id = r.CAN_ID

    filepath = joinpath(sim_folder, "$id.csv")

    res = CSV.read(filepath, DataFrame)

    target = target_date[id]

    p[i] = KidneyAllocation.probability_qualifying_offer_by(res, target; n_simulations=1000)

end


KidneyAllocation.auc(gt, p)
KidneyAllocation.brier_score(gt, p)



