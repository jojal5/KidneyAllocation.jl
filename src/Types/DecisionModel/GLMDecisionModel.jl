struct GLMDecisionModel <: AbstractDecisionModel
    fm::StatsModels.TableRegressionModel
    threshold::Real
end


function acceptance_probability(dm::GLMDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    arrival = donor.arrival
    n = length(recipients)

    # Preallocate columns (types matter)
    DON_AGE  = Vector{Int64}(undef, n)
    KDRI     = Vector{Float64}(undef, n)
    CAN_AGE  = Vector{Int64}(undef, n)
    CAN_WAIT = Vector{Float64}(undef, n)
    CAN_BLOOD = Vector{ABOGroup}(undef, n)
    CPRA = Vector{Int64}(undef, n)

    MISMATCH = Vector{Int64}(undef, n)
    DON_CAN_SCORE = Vector{Float64}(undef, n)

    for (i, r) in enumerate(recipients)

        DON_AGE[i]  = donor.age
        KDRI[i]     = donor.kdri

        CAN_AGE[i]  = years_between(r.birth, arrival)
        CAN_WAIT[i] = fractionalyears_between(r.dialysis, arrival)
        CAN_BLOOD[i] = r.blood
        CPRA[i] = r.cpra

        MISMATCH[i] = mismatch_count(donor, r)
        DON_CAN_SCORE[i] = score(donor, r)
        
    end

    df = DataFrame(DON_AGE = DON_AGE,
        KDRI = KDRI,
        CAN_AGE = CAN_AGE,
        CAN_WAIT = CAN_WAIT,
        CAN_BLOOD = CAN_BLOOD,
        CPRA = CPRA,
        MISMATCH = MISMATCH,
        DON_CAN_SCORE = DON_CAN_SCORE,
    )

    return GLM.predict(dm.fm, df)
end

function acceptance_probability(dm::GLMDecisionModel, recipient::Recipient, donor::Donor)
    p = acceptance_probability(dm, [recipient], donor::Donor)
    return p[1]
end