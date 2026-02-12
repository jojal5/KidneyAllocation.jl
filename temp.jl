using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2, Random

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, load_recipient, is_active, is_expired, is_abo_compatible
import KidneyAllocation: load_donor, build_donor_registry
import KidneyAllocation: shift_recipient_timeline, set_donor_arrival
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision
import KidneyAllocation: score, years_between, fractionalyears_between
import KidneyAllocation: allocate_one_donor, allocate
import KidneyAllocation: decide, fit_decision_threshold

recipients_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipients = build_recipient_registry(recipients_filepath, cpra_filepath)


donors_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

donors = build_donor_registry(donors_filepath)


# Fit decision model
data = retrieve_decision_data(donors_filepath, recipients_filepath)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE + MISMATCH)

fm = glm(model, data, Bernoulli(), LogitLink())
u = fit_decision_threshold(fm)


dm = GLMDecisionModel(fm)

u = fit_decision_threshold(dm)

r = Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 39, 77, 15, 17, 0)
d = Donor(Date(2000,1,1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5)



@time decide(dm, u, r, d)



# ------------------------------------------------------------------------------------

















