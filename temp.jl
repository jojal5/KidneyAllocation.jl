using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM, JLD2, Random

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, load_recipient, is_active, is_expired, is_abo_compatible
import KidneyAllocation: load_donor, build_donor_registry
import KidneyAllocation: shift_recipient_timeline, set_donor_arrival
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision
import KidneyAllocation: score, years_between, fractionalyears_between
import KidneyAllocation: allocate_one_donor, allocate


recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipients = build_recipient_registry(recipient_filepath, cpra_filepath)


donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

donors = build_donor_registry(donor_filepath)


# Fit decision model
data = retrieve_decision_data(donor_filepath, recipient_filepath)



model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

fm = glm(model, data, Bernoulli(), LogitLink())

u = fit_decision_threshold(fm)
# ------------------------------------------------------------------------------------












"""
    mismatch_count(d::TransplantEntity, r::TransplantEntity)::Int

Nombre de mismatches positionnels sur les loci A, B et DR (2 allèles chacun).
Compare (a1,a2,b1,b2,dr1,dr2) du donneur et du receveur.
Retourne un entier entre 0 (tous match) et 6 (aucun match).
"""




recipients = [
        Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 39, 77, 15, 17, 0),
        Recipient(Date(1981,1,1), Date(1997,1,1), Date(2000,6,1), A, 69, 2403, 7, 35, 4, 103, 10),
        Recipient(Date(1963,1,1), Date(1998,1,1), Date(2001,5,1), B, 25, 68, 67, 5102, 11, 16, 20),
    ]
donors = [
        Donor(Date(2000,1,1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5),
        Donor(Date(2001,1,1), 55, A, 2, 33, 37, 53, 4, 11, 1.6),
    ]


function get_HLA_A(t::TransplantEntity)
    return (t.a1, t.a2)
end

function get_HLA_B(t::TransplantEntity)
    return (t.b1, t.b2)
end

function get_HLA_DR(t::TransplantEntity)
    return (t.dr1, t.dr2)
end

"""
    mismatch_locus(l₁, l₂) -> Int

Return the number of allele mismatches between two HLA loci.

Counts how many alleles in `l₂` are not present in `l₁`.
"""
function mismatch_locus(l₁::Tuple{HLA, HLA}, l₂::Tuple{HLA, HLA})::Int
    return count(x -> !(x in l₁), l₂)
end

"""
    mismatch_count(t₁, t₂) -> Int

Return the total number of HLA mismatches between two transplant entities across loci A, B, and DR.
"""
function mismatch_count(t₁::TransplantEntity, t₂::TransplantEntity)::Int
    mm_A  = mismatch_locus(get_HLA_A(t₁),  get_HLA_A(t₂))
    mm_B  = mismatch_locus(get_HLA_B(t₁),  get_HLA_B(t₂))
    mm_DR = mismatch_locus(get_HLA_DR(t₁), get_HLA_DR(t₂))

    return mm_A + mm_B + mm_DR
end


@time mismatch_locus(get_HLA_A(recipients[1]),get_HLA_A(donors[1]))
mismatch_count(recipients[1], donors[1])

using Test

r = Recipient(Date(1979,1,1), Date(1995,1,1), Date(1998,1,1), O, 68, 203, 73, 77, 15, 17, 0)
        d = Donor(Date(2000,1,1), 40, O, 34, 3401, 73, 77, 3, 17, 1.5)

        @test mismatch_locus(get_HLA_A(r),get_HLA_A(d)) == 2
        @test mismatch_locus(get_HLA_B(r),get_HLA_B(d)) == 0
        @test mismatch_locus(get_HLA_DR(r),get_HLA_DR(d)) == 1
        @test mismatch_count(r, d) == 3









struct GLMDecisionModel <: AbstractDecisionModel
    fm::StatsModels.TableRegressionModel
    df::DataFrame
    blood_str::Dict{Any,String}
end

function GLMDecisionModel(fm::StatsModels.TableRegressionModel)
    # 1-row table with the predictor columns used by the model
    df = DataFrame(
        KDRI      = Float64[1.0],
        CAN_AGE   = Float64[50.0],
        CAN_WAIT  = Float64[2.0],
        CAN_BLOOD = String["O"],
        DON_AGE   = Float64[50.0],
    )

    # Cache of bloodtype -> string to avoid allocations
    blood_str = Dict{Any,String}(O=>"O", A=>"A", B=>"B", AB=>"AB")

    return GLMDecisionModel(fm, df, blood_str)
end


# You can keep a shared interface for all decision models
decision_probability(dm::GLMDecisionModel, donor::Donor, recipient::Recipient)::Float64 = begin
    arrival = donor.arrival

    @inbounds begin
        dm.df.KDRI[1]     = float(donor.kdri)   # used both as KDRI and log(KDRI) internally
        dm.df.CAN_AGE[1]  = float(years_between(recipient.birth, arrival))
        dm.df.CAN_WAIT[1] = float(fractionalyears_between(recipient.dialysis, arrival))
        dm.df.CAN_BLOOD[1] = dm.blood_str[get_bloodtype(recipient)]
        dm.df.DON_AGE[1]  = float(donor.age)
    end

    # For GLM Bernoulli(LogitLink), predict returns probabilities
    return predict(dm.fm, dm.df)[1]
end

decide(dm::GLMDecisionModel, donor::Donor, recipient::Recipient, u::Real)::Bool =
    decision_probability(dm, donor, recipient) > u