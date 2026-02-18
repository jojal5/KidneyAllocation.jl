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


"""
    fit_threshold_f1(gt, p) -> Float64

Return the classification threshold that maximizes the F1 score.

Candidate thresholds are `u ∈ {0, unique(p), 1}`. Predictions use the rule
`p ≥ u`. The vector `gt` must contain boolean ground-truth labels (`true` is
the positive class).
"""
function fit_threshold_f1(gt::AbstractVector{<:Bool}, p::AbstractVector{<:Real})
    @assert length(gt) == length(p) "Vector lengths must match, got $(length(gt)) ≠ $(length(p))."

    # Validate probabilities
    for pi in p
        (0.0 <= pi <= 1.0) || throw(ArgumentError("All probabilities must be in [0,1]. Got p = $(pi)."))
    end

    u = [1.; Float64.(sort(unique(p), rev=true)); 0.]
    f1 = Vector{Float64}(undef, length(u))

    k = count(gt)

    for i in eachindex(u)

        uᵢ = u[i]

        pred = p .≥ uᵢ

        tp = count(gt .&& pred)
        fp = count(.!(gt) .&& pred)
        fn = k - tp

        denom = 2tp + fp + fn

        f1[i] = denom == 0 ? 0.0 : (2tp / denom)
    end

    ind = argmax(f1)
    return u[ind]
end


"""
    fit_threshold_prevalence(gt, p) -> Float64

Return a threshold `u` that matches the observed prevalence as closely as possible
under the decision rule `p ≥ u`. In the presence of ties, the lower threshold is prefered for equally close prevalences.
"""
function fit_threshold_prevalence(
    gt::AbstractVector{<:Bool},
    p::AbstractVector{<:Real},
)
    @assert length(gt) == length(p) "Vector lengths must match, got $(length(gt)) ≠ $(length(p))."

    # Validate probabilities
    for pi in p
        (0.0 <= pi <= 1.0) || throw(ArgumentError("All probabilities must be in [0,1]. Got p = $(pi)."))
    end

    n = length(gt)
    k = count(gt)

    # Candidate thresholds for rule p ≥ u (tie-aware)
    u = unique!(sort!(Float64.(vcat(0.0, p, 1.0))))

    best_u = u[1]
    best_err = typemax(Int)

    for ui in u
        npos = count(>=(ui), p)      # number predicted positive under p ≥ ui
        err = abs(npos - k)          # closeness to observed prevalence
        if err < best_err
            best_err = err
            best_u = ui
        end
    end

    return best_u
end


