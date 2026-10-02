

struct TreeDecisionModel <: AbstractDecisionModel
    fm::DecisionTree.DecisionTreeClassifier
    features::Vector{Symbol}
    threshold::Real
end

function construct_feature_matrix(dm::TreeDecisionModel, recipients::Vector{Recipient}, donor::Donor)
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

    is_bloodtype_O = CAN_BLOOD .== O
    is_bloodtype_A = CAN_BLOOD .== A
    is_bloodtype_B = CAN_BLOOD .== B
    is_bloodtype_AB = CAN_BLOOD .== AB

    cols = Dict{Symbol,AbstractVector}(
        :DON_AGE => DON_AGE,
        :KDRI => KDRI,
        :CAN_AGE => CAN_AGE,
        :CAN_WAIT => CAN_WAIT,
        :CPRA => CPRA,
        :MISMATCH => MISMATCH,
        :DON_CAN_SCORE => DON_CAN_SCORE,
        :is_bloodtype_O => is_bloodtype_O,
        :is_bloodtype_A => is_bloodtype_A,
        :is_bloodtype_B => is_bloodtype_B,
        :is_bloodtype_AB => is_bloodtype_AB,
    )

    X = hcat((cols[f] for f in dm.features)...)

    return X
end


function acceptance_probability(dm::TreeDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    
    X = construct_feature_matrix(dm, recipients, donor)

    p = DecisionTree.predict_proba(dm.fm, X)

    return p[:,2]
end


