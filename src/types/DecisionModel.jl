abstract type AbstractDecisionModel end

include("DecisionModel/GLMDecisionModel.jl")
include("DecisionModel/TreeDecisionModel.jl")

function construct_feature_matrix_from_df(df::DataFrame, features::Vector{Symbol})
    cols = Vector{Vector{Float64}}(undef, length(features))

    for (j, f) in pairs(features)
        if f === :is_bloodtype_O
            cols[j] = Float64.(df.CAN_BLOOD .== "O")
        elseif f === :is_bloodtype_A
            cols[j] = Float64.(df.CAN_BLOOD .== "A")
        elseif f === :is_bloodtype_B
            cols[j] = Float64.(df.CAN_BLOOD .== "B")
        elseif f === :is_bloodtype_AB
            cols[j] = Float64.(df.CAN_BLOOD .== "AB")
        else
            cols[j] = Float64.(df[!, f])   # avoids copying the column first
        end
    end

    return hcat(cols...)
end

