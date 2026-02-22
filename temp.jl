using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: load_donor

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
df = load_donor(donor_filepath)


"""
    fill_hla_pair!(df, col1, col2)

Replace `missing` values in `col1` (resp. `col2`) by the value in `col2`
(resp. `col1`). Operates in place.
"""
function fill_hla_pair!(df::AbstractDataFrame, col1::Symbol, col2::Symbol)
    df[!, col1] = coalesce.(df[!, col1], df[!, col2])
    df[!, col2] = coalesce.(df[!, col2], df[!, col1])
    return df
end

"""
    fill_hla_pairs!(df, prefix)

Fill missing HLA allele values (A, B, DR) from their paired column in place.
"""
function fill_hla_pairs!(df::AbstractDataFrame, prefix::AbstractString)
    loci = ["A", "B" , "DR"]
    for locus in loci
        col1 = Symbol("$(prefix)_$(locus)1")
        col2 = Symbol("$(prefix)_$(locus)2")
        fill_hla_pair!(df, col1, col2)
    end
    return df
end

@testset "fill_hla_pairs" begin

    df = DataFrame(DON_A1 = [2, 2 ,2, missing], DON_A2= [3, missing ,3, 3], DON_B1 =[missing, 4, 4, 4], DON_B2=[5, 5, missing, missing], DON_DR1=[missing, 6 ,6 ,6], DON_DR2=[7, 7, 7, missing])

    fill_hla_pairs!(df, "DON")

    @test df.DON_A1 == [2, 2, 2, 3]
    @test df.DON_A2 == [3, 2, 3, 3]
    @test df.DON_B1 == [5, 4, 4, 4]
    @test df.DON_B2 == [5, 5, 4, 4]
    @test df.DON_DR1 == [7, 6, 6, 6]
    @test df.DON_DR2 == [7, 7, 7, 6]
end