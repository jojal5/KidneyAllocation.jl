abstract type AbstractDecisionModel end

include("DecisionModel/GLMDecisionModel.jl")
include("DecisionModel/TreeDecisionModel.jl")

function decide(dm::M, threshold::Real, recipient::Recipient, donor::Donor) where {M<:AbstractDecisionModel}
    
    # Acceptance probability
    p = acceptance_probability(dm, recipient, donor)

    # Decision
    return p>threshold
end