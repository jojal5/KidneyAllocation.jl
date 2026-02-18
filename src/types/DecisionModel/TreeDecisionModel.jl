

struct TreeDecisionModel <: AbstractDecisionModel
    fm::DecisionTree.DecisionTreeClassifier
    features::Vector{Symbol}
end

function construct_feature_matrix(dm::TreeDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    arrival = donor.arrival
    n = length(recipients)

    DON_AGE  = Vector{Float64}(undef, n)
    KDRI     = Vector{Float64}(undef, n)
    CAN_AGE  = Vector{Float64}(undef, n)
    CAN_WAIT = Vector{Float64}(undef, n)
    MISMATCH = Vector{Float64}(undef, n)
    is_bloodtype_O  = zeros(Float64, n)
    is_bloodtype_A  = zeros(Float64, n)
    is_bloodtype_B  = zeros(Float64, n)
    is_bloodtype_AB = zeros(Float64, n)

    kdri = float(donor.kdri)
    don_age = float(donor.age)

    fill!(DON_AGE, don_age)
    fill!(KDRI, kdri)

    @inbounds for i in 1:n
        r = recipients[i]

        CAN_AGE[i]  = float(years_between(get_birth(r), arrival))
        CAN_WAIT[i] = float(fractionalyears_between(get_dialysis(r), arrival))
        MISMATCH[i] = float(mismatch_count(donor, r))

        # One-hot blood type (exactly one is 1.0)
        if r.blood == O
            is_bloodtype_O[i] = 1.0
        elseif r.blood == A
            is_bloodtype_A[i] = 1.0
        elseif r.blood == B
            is_bloodtype_B[i] = 1.0
        else
            is_bloodtype_AB[i] = 1.0
        end
    end

    cols = Dict{Symbol,AbstractVector}(
        :DON_AGE => DON_AGE,
        :KDRI => KDRI,
        :CAN_AGE => CAN_AGE,
        :CAN_WAIT => CAN_WAIT,
        :MISMATCH => MISMATCH,
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

function acceptance_probability(dm::TreeDecisionModel, recipient::Recipient, donor::Donor)::Float64
    return acceptance_probability(dm, Recipient[recipient], donor)[1]
end

function decide(dm::TreeDecisionModel, recipients::Vector{Recipient}, donor::Donor)
    
    X = construct_feature_matrix(dm, recipients, donor)

    acceptation = DecisionTree.predict(dm.fm, X) .== 1

    return acceptation

end
