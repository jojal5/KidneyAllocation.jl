struct GLMDecisionModel <: AbstractDecisionModel
    fm::StatsModels.TableRegressionModel
    threshold::Real
end


function acceptance_probability(dm::GLMDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    arrival = donor.arrival
    n = length(recipients)

    # Preallocate columns (types matter)
    KDRI     = Vector{Float64}(undef, n)
    DON_AGE  = Vector{Int64}(undef, n)
    CAN_AGE  = Vector{Int64}(undef, n)
    CAN_WAIT = Vector{Float64}(undef, n)
    CAN_BLOOD = Vector{String}(undef, n)
    MISMATCH = Vector{Int64}(undef, n) 

    kdri = donor.kdri
    don_age = donor.age

    @inbounds for i in 1:n
        r = recipients[i]

        KDRI[i]     = kdri
        DON_AGE[i]  = don_age
        CAN_AGE[i]  = years_between(get_birth(r), arrival)
        CAN_WAIT[i] = fractionalyears_between(get_dialysis(r), arrival)
        CAN_BLOOD[i] = string(r.blood)
        MISMATCH[i] = mismatch_count(donor, r)
    end

    df = DataFrame(
        KDRI = KDRI,
        CAN_AGE = CAN_AGE,
        CAN_WAIT = CAN_WAIT,
        CAN_BLOOD = CAN_BLOOD,
        DON_AGE = DON_AGE,
        MISMATCH = MISMATCH,
    )

    return predict(dm.fm, df)
end

function acceptance_probability(dm::GLMDecisionModel, recipient::Recipient, donor::Donor)
    p = acceptance_probability(dm, [recipient], donor::Donor)
    return p[1]
end


function decide(dm::GLMDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    
    # Acceptance probability
    p = acceptance_probability(dm, recipients, donor)

    # Decision
    return p .> dm.threshold
end