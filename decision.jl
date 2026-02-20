using Pkg
Pkg.activate(".")

using KidneyAllocation

using Dates, DataFrames, DecisionTree, GLM

import KidneyAllocation: retrieve_decision_data, fit_threshold_f1, fit_threshold_prevalence, auc, brier_score


## Load data

recipients_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
donors_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

data = retrieve_decision_data(donors_filepath, recipients_filepath)

data_train = filter(row -> row.LEARNING_SET == "train", data)
data_validation = filter(row -> row.LEARNING_SET == "validation", data)


## Fit decison model based on GLM on the train set
model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

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




## Fit decision model based on classification tree on the training set

features = Symbol.([
    "DON_AGE"
    "KDRI"
    "CAN_AGE"
    "CAN_WAIT"
    "MISMATCH"
    "is_bloodtype_O"
    "is_bloodtype_A"
    "is_bloodtype_B"
    "is_bloodtype_AB"])

m = DecisionTreeClassifier(
    max_depth=10, min_samples_leaf=125,
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

