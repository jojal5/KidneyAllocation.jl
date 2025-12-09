module KidneyAllocation

    using Dates

    include("types/ABOGroup.jl")
    include("types/TransplantEntity.jl")
    include("AllocationScore/allocationscore.jl")
    include("Quality/quality.jl")
    include("utils.jl")

    export 

    ABOGroup, O, A, B, AB,
    HLA,
    TransplantEntity, Donor, Recipient

end
