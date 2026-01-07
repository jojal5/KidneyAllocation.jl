using Pkg
Pkg.activate(".")

using KidneyAllocation

using CSV, DataFrames, Dates, GLM, Missings

"""
    auc(gt::Array{<:Real}, scores::Array{<:Real})

Compute the area under the ROC curve based on the ground truth `gt` and the success probability `scores`.

See also `roc()` of MLBase.
"""
function auc(gt::Array{<:Real},scores::Array{<:Real})

    # Compute the ROC curve for 100 equally spaced thresholds - see `roc()`
    r = roc(gt, scores, 0:.01:1)

    # Compute the true positive rate and false positive rate
    tpr = true_positive_rate.(r)
    fpr = false_positive_rate.(r)

    # Numerical computation of the area under the ROC curve
    p = sortperm(fpr)

    permute!(tpr,p)
    permute!(fpr,p)

    area = 0.0

    for i in 2:length(tpr)
        dx = fpr[i] - fpr[i-1]
        dy = tpr[i] - tpr[i-1]
        area += dx*tpr[i-1] + dx*dy/2
    end

    return area

end

filename = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
donors = CSV.read(filename, DataFrame, missingstring="NULL")

filter!(row->row.ATT_TYPE == 5, donors)
select!(donors, :DON_ID, :DON_DEATH_TM, :CAN_ID, :DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD, :DECISION)

filter!(row->year(row.DON_DEATH_TM) ∈ 2014:2020, donors)

dropmissing!(donors)
filter!(row->row.WEIGHT > 0., donors)

kdri = Float64[]

# r = eachrow(donors)[2]
for r in eachrow(donors)

    age = r.DON_AGE
    height = r.HEIGHT
    weight = r.WEIGHT
    hypertension = r.HYPERTENSION == 1
    diabetes = r.DIABETES == 1
    cva = (r.DEATH == 4) || (r.DEATH ==16)
    creatinine = KidneyAllocation.creatinine_mgdl(r.CREATININE)
    dcd = df.DCD[ind] .== 1 # TODO À VÉRIFIER si c'est bien 1, sinon c'est 2 (Anastasiay a confirmé le code)

    kdri_r = KidneyAllocation.evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

    push!(kdri, kdri_r)
end

donors.KDRI = kdri

data = select!(donors, :DON_ID, :DON_DEATH_TM, :DON_AGE, :KDRI, :CAN_ID, :DECISION)

filename = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
candidates = CSV.read(filename, DataFrame, missingstring="NULL")

select!(candidates, :CAN_ID, :CAN_BTH_DT, :CAN_BLOOD, :CAN_DIAL_DT)

age = Int64[]
waittime = Union{Float64, Missing}[]
blood = String[]

r = eachrow(data)[1]
for r in eachrow(data)

    ind = findfirst(candidates.CAN_ID .== r.CAN_ID)
    age_r = KidneyAllocation.years_between(Date(candidates.CAN_BTH_DT[ind]), Date(r.DON_DEATH_TM))

    if ismissing(candidates.CAN_DIAL_DT[ind])
        waittime_r = missing
    else
        waittime_r = KidneyAllocation.fractionalyears_between(Date(candidates.CAN_DIAL_DT[ind]), Date(r.DON_DEATH_TM))
    end
        blood_r = candidates.CAN_BLOOD[ind]

    push!(age, age_r)
    push!(waittime, waittime_r)
    push!(blood, blood_r)

end

data.CAN_AGE = age
data.CAN_WAIT = waittime
data.CAN_BLOOD = blood

dropmissing!(data)
data.DECISION = data.DECISION .== "Acceptation"

using MLBase, Optim


filter!(row->row.CAN_WAIT > 0, data)

gt = Int64.(data.DECISION)

fm = glm(@formula(DECISION ~ KDRI + CAN_AGE + CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ KDRI + CAN_AGE + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)


data.CAN_YOUTH = data.CAN_AGE .< 50
fm = glm(@formula(DECISION ~ KDRI + CAN_YOUTH + CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ KDRI + CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)


data.LOG_CAN_WAIT = log.(data.CAN_WAIT)
fm = glm(@formula(DECISION ~ KDRI + LOG_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)


data.EXP_KDRI = exp.(data.KDRI)
fm = glm(@formula(DECISION ~ EXP_KDRI + CAN_YOUTH + CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

data.LOG_KDRI = log.(data.KDRI)
fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_YOUTH + CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ KDRI + CAN_WAIT  + CAN_AGE + CAN_AGE*KDRI + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)


data.SQUARE_CAN_AGE = (data.CAN_AGE).^2
fm = glm(@formula(DECISION ~ KDRI + CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

# data.SQUARE_CAN_AGE = (data.CAN_AGE).^2
fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

data.SQUARE_CAN_WAIT = (data.CAN_WAIT).^2
fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI +SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*KDRI + SQUARE_CAN_AGE*KDRI +SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*KDRI + SQUARE_CAN_AGE*KDRI +SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*KDRI*CAN_WAIT + SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*LOG_KDRI + SQUARE_CAN_AGE*LOG_KDRI +SQUARE_CAN_AGE*LOG_KDRI*SQUARE_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*LOG_KDRI + SQUARE_CAN_AGE*LOG_KDRI +SQUARE_CAN_AGE*LOG_KDRI*SQUARE_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)



fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + CAN_AGE + CAN_WAIT*SQUARE_CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)


fm = glm(@formula(DECISION ~ KDRI + CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + SQUARE_CAN_AGE*KDRI + CAN_AGE*KDRI + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)








fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE*KDRI + SQUARE_CAN_AGE*KDRI +SQUARE_CAN_AGE*KDRI*SQUARE_CAN_WAIT + CAN_BLOOD), data, Bernoulli(), LogitLink())
θ̂ = predict(fm)
auc(gt, θ̂)

fobj(u::Real) = -f1score(roc(gt, θ̂, u))

res = optimize(fobj, .1, .5)

u = res.minimizer

r = roc(gt, θ̂, u)




β̂ = coef(fm)

lp = β̂[1]
