module KidneyAllocation

    using CSV, DataFrames, Dates, Random
    using StatsModels # Needed when passing a glm fitted model as a function argument
    using MLBase, Optim

    include("types/ABOGroup.jl")
    include("types/TransplantEntity.jl")
    include("AllocationScore/allocationscore.jl")
    include("Quality/quality.jl")
    include("decision.jl")
    include("preprocess.jl")
    include("simulation.jl")
    include("utils.jl")

    export 

    ABOGroup, O, A, B, AB,
    HLA,
    TransplantEntity, Donor, Recipient

end
