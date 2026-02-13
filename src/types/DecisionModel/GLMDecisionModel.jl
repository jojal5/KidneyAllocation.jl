struct GLMDecisionModel <: AbstractDecisionModel
    fm::StatsModels.TableRegressionModel
    df::DataFrame
    blood_str::Dict{Any,String}
end

function GLMDecisionModel(fm::StatsModels.TableRegressionModel)
    # 1-row table with the predictor columns used by the model
    df = DataFrame(
        KDRI      = Float64[1.0],
        CAN_AGE   = Int64[50],
        CAN_WAIT  = Float64[2.0],
        CAN_BLOOD = String["O"],
        DON_AGE   = Int64[50],
        MISMATCH  = Int64[6]
    )

    # Cache of bloodtype -> string to avoid allocations
    blood_str = Dict{Any,String}(O=>"O", A=>"A", B=>"B", AB=>"AB")

    return GLMDecisionModel(fm, df, blood_str)
end

function acceptance_probability(dm::GLMDecisionModel, recipient::Recipient, donor::Donor)
    arrival = donor.arrival

    mm = mismatch_count(donor, recipient)

    @inbounds begin
        dm.df.KDRI[1]     = donor.kdri
        dm.df.CAN_AGE[1]  = years_between(recipient.birth, arrival)
        dm.df.CAN_WAIT[1] = fractionalyears_between(recipient.dialysis, arrival)
        dm.df.CAN_BLOOD[1] = dm.blood_str[get_bloodtype(recipient)]
        dm.df.DON_AGE[1]  = donor.age
        dm.df.MISMATCH[1] = mm
    end

    # Acceptance probability
    p = predict(dm.fm, dm.df)[1]

    return p
end

function fit_decision_threshold(dm::GLMDecisionModel)

    fm = dm.fm

    gt = Int64.(response(fm))

    θ̂ = predict(fm)

    fobj(u::Real) = -f1score(roc(gt, θ̂, u))

    res = optimize(fobj, 0.01, 0.75)

    u = res.minimizer

    return u

end