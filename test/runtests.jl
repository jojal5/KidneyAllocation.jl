using KidneyAllocation
using Test

using Dates

@testset "KidneyAllocation.jl" begin
    include("test_types.jl")
    include("test_utils.jl")
    include("Quality/test_quality.jl")
    include("AllocationScore/test_allocationscore.jl")
end
