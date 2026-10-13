using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

# import KidneyAllocation: build_recipient_registry, build_donor_registry

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"


df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)


# Filter the recipients and donors to the period from 2016 to 2019 inclusively
filter!(row -> row.CAN_LISTING_DT ≥ DateTime(2016,1,1) && row.CAN_LISTING_DT ≥ DateTime(2019,12,31), df_recipients)
filter!(row -> row.DON_DEATH_TM ≥ DateTime(2016,1,1) && row.DON_DEATH_TM ≥ DateTime(2019,12,31), df_donors)


donor_registry, given_kidneys = KidneyAllocation.build_donor_registry(df_donors) 
recipient_registry, time_to_first_offer_registry, time_to_transplant_registry = KidneyAllocation.build_recipient_registry(df_recipients)

don_id = collect(keys(donor_registry))
can_id = collect(keys(recipient_registry))

# Filter the dataframe for retrieving donors with all the information
df = filter(:DON_ID => id -> id ∈ don_id, df_donors)
filter!(:CAN_ID => id -> id ∈ can_id, df)

m = nrow(df)

don_age = Vector{Int64}(undef, m)
kdri = Vector{Float64}(undef, m)

can_age = Vector{Int64}(undef, m)
can_wait = Vector{Float64}(undef, m)
can_blood = Vector{ABOGroup}(undef, m)

cpra = Vector{Int64}(undef, m)

mismatch = Vector{Int64}(undef, m)
score = Vector{Float64}(undef, m)

decision = falses(m)
learning_set = Vector{String}(undef, m)

n_train = round(Int, length(don_id) * 0.8)
rng = MersenneTwister(12345)
train_id = Set(shuffle(rng, don_id)[1:n_train])

for i in 1:m
    d = donor_registry[df.DON_ID[i]]
    don_age[i] = d.age
    kdri[i] = d.kdri

    r = recipient_registry[df.CAN_ID[i]]
    can_age[i] = KidneyAllocation.years_between(r.birth, d.arrival)
    can_wait[i] = KidneyAllocation.fractionalyears_between(r.dialysis, d.arrival)
    can_blood[i] = r.blood
    cpra[i] = r.cpra

    mismatch[i] = KidneyAllocation.mismatch_count(d, r)
    score[i] = df.DON_CAN_SCORE[i]

    decision[i] = uppercase(strip(df.DECISION[i])) == "ACCEPTATION"

    learning_set[i] = df.DON_ID[i] ∈ train_id ? "training" : "validation"

end

data = DataFrame(DON_AGE = don_age,
    KDRI = kdri,
    CAN_AGE = can_age,
    CAN_WAIT = can_wait,
    CAN_BLOOD = can_blood,
    CPRA = cpra,
    MISMATCH = mismatch,
    DON_CAN_SCORE = score,
    DECISION = decision,
    LEARNING_SET = learning_set
)

data_train = filter(row -> row.LEARNING_SET == "training", data)
data_validation = filter(row -> row.LEARNING_SET == "validation", data)

using GLM
import KidneyAllocation: auc, brier_score

# model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + CPRA + DON_CAN_SCORE)
# model = @formula(DECISION ~ log(KDRI) + CAN_BLOOD + CAN_WAIT + CAN_WAIT^2 + CPRA + CPRA^2 + CAN_AGE + CAN_AGE^2 + DON_CAN_SCORE)
# model = @formula(DECISION ~ KDRI + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + CPRA + DON_CAN_SCORE)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE + CAN_BLOOD + DON_CAN_SCORE)

fm = glm(model, data_train, Bernoulli(), LogitLink())

## Performance on the train set

gt = response(fm) .≈ 1.
p = GLM.predict(fm)

auc(gt, p)
brier_score(gt, p)

## Performance on the validation set
gt = data_validation.DECISION
p = float.(GLM.predict(fm, data_validation))

auc(gt, p)
brier_score(gt, p)

## Refit the GLM model on all the data and save it for later use

fm = glm(model, data, Bernoulli(), LogitLink())

gt = data.DECISION
p = float.(GLM.predict(fm, data))

auc(gt, p)
brier_score(gt, p)

u = KidneyAllocation.fit_threshold_f1(data.DECISION, GLM.predict(fm))

dm = GLMDecisionModel(fm, u)

jldsave("src/SyntheticData/GLMDecisionModel.jld2"; dm)

## Fit decision model based on classification tree on the training set

using DecisionTree

features = Symbol.([
    "DON_AGE"
    "KDRI"
    "CAN_AGE"
    "CAN_WAIT"
    "MISMATCH"
    "CPRA"
    "DON_CAN_SCORE"
    "is_bloodtype_O"
    "is_bloodtype_A"
    "is_bloodtype_B"
    "is_bloodtype_AB"])

m = DecisionTreeClassifier(
    max_depth=10, min_samples_leaf=175,
    pruning_purity_threshold=1
)

X = KidneyAllocation.construct_feature_matrix_from_df(data_train, features)
y = data_train.DECISION

DecisionTree.fit!(m, X, y)

## Performance on the train set

gt = y
p = DecisionTree.predict_proba(m, X)[:, 2]

auc(gt, p)
brier_score(gt, p)

## Performance on the validation set

X = KidneyAllocation.construct_feature_matrix_from_df(data_validation, features)
y = data_validation.DECISION

gt = y
p = DecisionTree.predict_proba(m, X)[:, 2]

auc(gt, p)
brier_score(gt, p)

## Refit the tree based model on all the data and save it for later use

features = Symbol.([
    "DON_AGE"
    "KDRI"
    "CAN_AGE"
    "CAN_WAIT"
    "CPRA"
    "MISMATCH"
    "DON_CAN_SCORE"
    "is_bloodtype_O"
    "is_bloodtype_A"
    "is_bloodtype_B"
    "is_bloodtype_AB"])

X = KidneyAllocation.construct_feature_matrix_from_df(data, features)
y = data.DECISION

DecisionTree.fit!(m, X, y)

gt = y
p = DecisionTree.predict_proba(m, X)[:, 2]

auc(gt, p)
brier_score(gt, p)

# u = fit_threshold_prevalence(data.DECISION, DecisionTree.predict_proba(m, X)[:,2])
u = KidneyAllocation.fit_threshold_f1(data.DECISION, DecisionTree.predict_proba(m, X)[:,2])

dm = TreeDecisionModel(m, features, u)

jldsave("src/SyntheticData/TreeDecisionModel.jld2"; dm)