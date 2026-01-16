using Pkg
Pkg.activate(".")

using KidneyAllocation

using Dates, GLM

"""
    auc(gt::Array{<:Real}, scores::Array{<:Real})

Compute the area under the ROC curve based on the ground truth `gt` and the success probability `scores`.

See also `roc()` of MLBase.
"""
function auc(gt::Array{<:Real}, scores::Array{<:Real})

    # Compute the ROC curve for 100 equally spaced thresholds - see `roc()`
    r = roc(gt, scores, 0:0.01:1)

    # Compute the true positive rate and false positive rate
    tpr = true_positive_rate.(r)
    fpr = false_positive_rate.(r)

    # Numerical computation of the area under the ROC curve
    p = sortperm(fpr)

    permute!(tpr, p)
    permute!(fpr, p)

    area = 0.0

    for i in 2:length(tpr)
        dx = fpr[i] - fpr[i-1]
        dy = tpr[i] - tpr[i-1]
        area += dx * tpr[i-1] + dx * dy / 2
    end

    return area

end

# Best logistic regression model
# fm = glm(@formula(DECISION ~ LOG_KDRI + CAN_WAIT + SQUARE_CAN_WAIT + CAN_AGE + SQUARE_CAN_AGE + CAN_AGE * KDRI + SQUARE_CAN_AGE * KDRI + SQUARE_CAN_AGE * KDRI * SQUARE_CAN_WAIT + CAN_BLOOD + DON_AGE), data, Bernoulli(), LogitLink())

import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision

recipients_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
donors_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

data = retrieve_decision_data(donors_filepath, recipients_filepath)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

fm = glm(model, data, Bernoulli(), LogitLink())

u = fit_decision_threshold(fm)

donor = Donor(Date(2010, 1, 1), 22, B, 3, 29, 7, 44, 7, 13, 0.8190073900205177)
recipient = Recipient(Date(1953, 10, 20), Date(2002, 4, 21), Date(2003, 4, 10), B,
    2, 23, 45, 65, 11, 16,
    0)

get_decision(donor, recipient, fm ,u)
