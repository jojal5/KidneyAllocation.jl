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
# ------------------------------------------------------------------------------------


















struct GLMDecisionModel <: AbstractDecisionModel
    fm::StatsModels.TableRegressionModel
    df::DataFrame
    blood_str::Dict{Any,String}
end

function GLMDecisionModel(fm::StatsModels.TableRegressionModel)
    # 1-row table with the predictor columns used by the model
    df = DataFrame(
        KDRI      = Float64[1.0],
        CAN_AGE   = Float64[50.0],
        CAN_WAIT  = Float64[2.0],
        CAN_BLOOD = String["O"],
        DON_AGE   = Float64[50.0],
    )

    # Cache of bloodtype -> string to avoid allocations
    blood_str = Dict{Any,String}(O=>"O", A=>"A", B=>"B", AB=>"AB")

    return GLMDecisionModel(fm, df, blood_str)
end


# You can keep a shared interface for all decision models
decision_probability(dm::GLMDecisionModel, donor::Donor, recipient::Recipient)::Float64 = begin
    arrival = donor.arrival

    @inbounds begin
        dm.df.KDRI[1]     = float(donor.kdri)   # used both as KDRI and log(KDRI) internally
        dm.df.CAN_AGE[1]  = float(years_between(recipient.birth, arrival))
        dm.df.CAN_WAIT[1] = float(fractionalyears_between(recipient.dialysis, arrival))
        dm.df.CAN_BLOOD[1] = dm.blood_str[get_bloodtype(recipient)]
        dm.df.DON_AGE[1]  = float(donor.age)
    end

    # For GLM Bernoulli(LogitLink), predict returns probabilities
    return predict(dm.fm, dm.df)[1]
end

decide(dm::GLMDecisionModel, donor::Donor, recipient::Recipient, u::Real)::Bool =
    decision_probability(dm, donor, recipient) > u