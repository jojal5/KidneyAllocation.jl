using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, DecisionTree, GLM, JLD2, Random, Statistics, Test

using KidneyAllocation

import KidneyAllocation: auc, brier_score


## Chargement des données 

recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

df_recipients, df_donors = KidneyAllocation.load_data(recipient_filepath, cpra_filepath, donor_filepath)

## build_donor_registry

donor_registry, given_kidneys = KidneyAllocation.build_donor_registry(df_donors)
recipient_registry = KidneyAllocation.build_recipient_registry(df_recipients)


## Récupération des décisions

df_offers = KidneyAllocation.load_donor(donor_filepath)
filter!(row -> row.ATT_TYPE == 5, df_offers)

df_offers_filtered = similar(df_offers, 0)

for df_offer in groupby(df_offers, :DON_ID)

    filtered_df = KidneyAllocation.filter_offers_until_last_acceptance(df_offer)
    append!(df_offers_filtered, filtered_df)

end



## Build decision dataset

df_decision = DataFrame(DON_ID = Int[], DON_AGE = Int[], KDRI = Float64[], CAN_AGE = Int[], CAN_BLOOD = ABOGroup[],  CAN_WAIT = Float64[], DON_CAN_SCORE = Float64[], CPRA = Int[], MISMATCH = Int[], DECISION = Bool[])

# Only patient between 2016 and 2020
df = filter(row -> row.DON_DEATH_TM ≥ Date(2016,1,1) && row.DON_DEATH_TM < Date(2020,1,1), df_offers_filtered)


for row in eachrow(df)
    CAN_ID = row.CAN_ID
    DON_ID = row.DON_ID

    if haskey(recipient_registry, CAN_ID) && haskey(donor_registry, DON_ID)

        recipient = recipient_registry[row.CAN_ID]
        donor = donor_registry[row.DON_ID]


        donor_age = donor.age
        kdri = donor.kdri
        recipient_age = KidneyAllocation.years_between(recipient.birth, donor.arrival)
        recipient_blood = recipient.blood
        waittime = KidneyAllocation.fractionalyears_between(recipient.dialysis, donor.arrival)
        s = row.DON_CAN_SCORE
        mismatch = KidneyAllocation.mismatch_count(recipient, donor)
        cpra = recipient.cpra
        decision = row.DECISION == "Acceptation"

        push!(df_decision, [DON_ID, donor_age, kdri, recipient_age, recipient_blood, waittime, s, cpra, mismatch, decision])

    end
end




id = unique(df_decision.DON_ID)
n_train = round(Int, length(id) * 2/3)
rng = MersenneTwister(12345)
train_id = Set(shuffle(rng, id)[1:n_train])

learning_set = Vector{String}(undef, nrow(df_decision))
for (i, row) in enumerate(eachrow(df_decision))
    learning_set[i] = row.DON_ID ∈ train_id ? "train" : "validation"
end

df_decision.LEARNING_SET = learning_set

df_decision.AGE_GAP = df_decision.CAN_AGE - df_decision.DON_AGE

select!(df_decision, :DON_AGE, :KDRI, :CAN_AGE, :CAN_BLOOD, :CAN_WAIT, :DON_CAN_SCORE, :CPRA, :MISMATCH, :AGE_GAP, :DECISION, :LEARNING_SET)

df_decision.is_bloodtype_A = df_decision.CAN_BLOOD .== A
df_decision.is_bloodtype_B = df_decision.CAN_BLOOD .== B
df_decision.is_bloodtype_AB = df_decision.CAN_BLOOD .== AB
df_decision.is_bloodtype_O = df_decision.CAN_BLOOD .== O


## GLM based decision model

data_train = filter(row -> row.LEARNING_SET == "train", df_decision)
data_validation = filter(row -> row.LEARNING_SET == "validation", df_decision)

# model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + CPRA + DON_CAN_SCORE)
# model = @formula(DECISION ~ log(KDRI) + CAN_BLOOD + CAN_WAIT + CAN_WAIT^2 + CPRA + CPRA^2 + CAN_AGE + CAN_AGE^2 + DON_CAN_SCORE)
# model = @formula(DECISION ~ KDRI + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + CPRA + DON_CAN_SCORE)

model = @formula(DECISION ~ log(KDRI)*CAN_AGE + CAN_BLOOD + CPRA + DON_CAN_SCORE)

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

fm = glm(model, df_decision, Bernoulli(), LogitLink())

gt = df_decision.DECISION
p = float.(GLM.predict(fm, df_decision))

auc(gt, p)
brier_score(gt, p)

u = KidneyAllocation.fit_threshold_f1(df_decision.DECISION, GLM.predict(fm))

dm = GLMDecisionModel(fm, u)

jldsave("src/SyntheticData/GLMDecisionModel.jld2"; dm)



## Tree based decision model

features = Symbol.([
    "DON_AGE"
    "KDRI"
    "CAN_AGE"
    "CAN_WAIT"
    "DON_CAN_SCORE"
    "CPRA"
    "MISMATCH"
    "AGE_GAP"
    "is_bloodtype_O"
    "is_bloodtype_A"
    "is_bloodtype_B"
    "is_bloodtype_AB"])

tree = DecisionTreeClassifier(max_depth = 3)

X = KidneyAllocation.construct_feature_matrix_from_df(data_train, features)
y = data_train.DECISION

DecisionTree.fit!(tree, X, y)
# print_tree(tree)

## Performance on the train set

p = DecisionTree.predict_proba(tree, X)[:, 2]

auc(y, p)
brier_score(y, p)

## Performance on the validation set

X = KidneyAllocation.construct_feature_matrix_from_df(data_validation, features)
y = data_validation.DECISION

p = DecisionTree.predict_proba(tree, X)[:, 2]

auc(y, p)
brier_score(y, p)


## Refit the tree model on all the data and save it for later use

X = KidneyAllocation.construct_feature_matrix_from_df(df_decision, features)
y = df_decision.DECISION

DecisionTree.fit!(tree, X, y)

p = DecisionTree.predict_proba(tree, X)[:, 2]

auc(df_decision.DECISION, p)
brier_score(df_decision.DECISION, p)

u = KidneyAllocation.fit_threshold_f1(df_decision.DECISION, p)

dm = TreeDecisionModel(tree, features, u)

jldsave("src/SyntheticData/TreeDecisionModel.jld2"; dm)
