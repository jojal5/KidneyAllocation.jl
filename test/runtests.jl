using KidneyAllocation
using Test

using CSV, DataFrames, Dates, JLD2, Random 

@testset "KidneyAllocation.jl" begin
    include("test_types.jl")
    include("test_decision.jl")
    include("test_utils.jl")
    include("test_preprocess.jl")
    include("test_simulation.jl")
    include("Quality/test_quality.jl")
    include("AllocationScore/test_allocationscore.jl")
end
