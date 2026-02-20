module KidneyAllocation

    using CSV, DataFrames, Dates, DecisionTree, Distributions, GLM, Random
    using StatsModels # Needed when passing a glm fitted model as a function argument

    include("Types/ABOGroup.jl")
    include("Types/TransplantEntity.jl")
    include("Types/DecisionModel.jl")
    include("AllocationScore/allocationscore.jl")
    include("Quality/quality.jl")
    include("allocation.jl")
    include("decision.jl")
    include("preprocess.jl")
    include("simulation.jl")
    include("utils.jl")

    export 

    ABOGroup, O, A, B, AB,
    HLA,
    TransplantEntity, Donor, Recipient,
    AbstractDecisionModel, GLMDecisionModel, TreeDecisionModel

end
