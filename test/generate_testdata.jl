using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, DecisionTree, JLD2

using KidneyAllocation


import KidneyAllocation: retrieve_decision_data, construct_feature_matrix_from_df, fit_threshold_prevalence

# Fit decision model (Classification tree based)
donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
data = retrieve_decision_data(donor_filepath, recipient_filepath)

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

 X = KidneyAllocation.construct_feature_matrix_from_df(data, features)

y = data.DECISION

m = DecisionTreeClassifier(
        max_depth=10, min_samples_leaf=50,
        pruning_purity_threshold=1
        )

DecisionTree.fit!(m, X, y)

gt = data.DECISION
p = DecisionTree.predict_proba(m, X)[:,2]

threshold = fit_threshold_prevalence(gt, p)

dm = TreeDecisionModel(m, features, threshold)

jldsave("test/data/tree_decision_model.jld2"; dm)

