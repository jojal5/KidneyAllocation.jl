using Pkg
Pkg.activate(".")

using Dates, CSV, DataFrames, Distributions, GLM

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, load_recipient, is_active, is_expired, is_abo_compatible
import KidneyAllocation: load_donor, build_donor_registry
import KidneyAllocation: shift_recipient_timeline, set_donor_arrival
import KidneyAllocation: retrieve_decision_data, fit_decision_threshold, get_decision



recipient_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipients = build_recipient_registry(recipient_filepath, cpra_filepath)

# Estimate the recipient arrival rate
df = load_recipient(recipient_filepath)
df2 = filter(row -> 2014 ≤ year(row.CAN_LISTING_DT) < 2020, df)
λᵣ = length(unique(df2.CAN_ID)) / 6

donor_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"

donors = build_donor_registry(donor_filepath)

# Estimate the donor arrival rate
df = load_donor(donor_filepath)
df2 = filter(row -> 2014 ≤ year(row.DON_DEATH_TM) < 2020, df)
filter!(row -> row.DECISION == "Acceptation", df2)
λₒ = length(unique(df2.CAN_ID)) / 6


# Fit decision model
data = retrieve_decision_data(donor_filepath, recipient_filepath)

model = @formula(DECISION ~ log(KDRI) + CAN_AGE * KDRI * CAN_WAIT + CAN_AGE^2 * KDRI * CAN_WAIT^2 + CAN_BLOOD + DON_AGE)

fm = glm(model, data, Bernoulli(), LogitLink())

u = fit_decision_threshold(fm)
# ------------------------------------------------------------------------------------



ind_active = is_active.(recipients, Date(2013, 12, 31))

waiting_recipients = recipients[ind_active]

nᵣ = rand(Poisson(λᵣ * 10))
tᵣ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nᵣ)

sampled_recipients = rand(recipients, nᵣ)
new_recipients = KidneyAllocation.shift_recipient_timeline.(sampled_recipients, tᵣ)

append!(waiting_recipients, new_recipients)

ind_active = is_active.(waiting_recipients, Date(2013, 12, 31))

waiting_recipients[ind_active]


nₒ = rand(Poisson(λₒ * 10))
tₒ = KidneyAllocation.sample_days(Date(2014, 1, 1), Date(2023, 12, 31), nₒ)

sampled_donors = rand(donors, nₒ)
new_donors = set_donor_arrival.(sampled_donors, tₒ)

KidneyAllocation.get_arrival.(new_donors)



donor = new_donors[1]

@time for donor in new_donors

    filter!(x->!is_expired(x, donor.arrival), waiting_recipients)

    id_recipient = eachindex(waiting_recipients)

    ind_active = is_active.(waiting_recipients, donor.arrival)
    ind_compatible = is_abo_compatible.(KidneyAllocation.get_bloodtype.(waiting_recipients), KidneyAllocation.get_bloodtype(donor))

    eligible_id = id_recipient[ind_active.&&ind_compatible]

    s = KidneyAllocation.score.(donor, waiting_recipients[eligible_id])

    ind = sortperm(s; rev=true)

    permute!(eligible_id, ind)

    accepting_recipient_id = 0

    for iter in eachindex(eligible_id)
        id = eligible_id[iter]
        decision = get_decision(donor, waiting_recipients[id], fm, u)
        if decision
            accepting_recipient_id = id
            break
        end
    end


    # Sanity checks
    # is_active(waiting_recipients[accepting_recipient_id ], donor.arrival)
    # is_abo_compatible.(KidneyAllocation.get_bloodtype.(waiting_recipients[accepting_recipient_id]), KidneyAllocation.get_bloodtype(donor))
    # KidneyAllocation.score.(donor, waiting_recipients[accepting_recipient_id])

    if accepting_recipient_id != 0
        deleteat!(waiting_recipients, accepting_recipient_id)
    end

end

KidneyAllocation.is_expired(waiting_recipients[1], Date(204,1,1))

waiting_recipients

is_active.(waiting_recipients, Date(2024,1,1))

waiting_recipients[is_active.(waiting_recipients, Date(2024,1,1))]

r = filter(x->!is_expired(x, Date(2024,1,1)), waiting_recipients)

arr = KidneyAllocation.get_arrival.(r)

minimum(arr)





eligible_id=Int[]
empty!(eligible_id)
sizehint!(eligible_id, length(waiting_recipients))






import KidneyAllocation.score

function allocate_one_donor!(
    waiting_recipients::Vector{Recipient},
    donor::Donor,
    fm,
    u;
    eligible_id::Vector{Int}=Int[],
    scores::Vector{Float64}=Float64[],
    perm::Vector{Int}=Int[])

    arrival = donor.arrival
    donor_bt = donor.blood

    empty!(eligible_id)
    empty!(scores)

    sizehint!(eligible_id, length(waiting_recipients))
    sizehint!(scores, length(waiting_recipients))

    @inbounds for i in eachindex(waiting_recipients)
        r = waiting_recipients[i]
        is_active(r, arrival) || continue
        is_abo_compatible(r.blood, donor_bt) || continue

        push!(eligible_id, i)
        push!(scores, score(donor, r))
    end

    isempty(eligible_id) && return 0

    # sortperm! with reusable buffer
    resize!(perm, length(scores))
    sortperm!(perm, scores; rev=true)

    accepting_id = 0
    @inbounds for k in perm
        id = eligible_id[k]
        if get_decision(donor, waiting_recipients[id], fm, u)
            accepting_id = id
            break
        end
    end

    # K = min(100, length(scores))
    # topk = partialsortperm(scores, 1:K; rev=true)

    # accepting_id = 0
    # @inbounds for k in topk
    #     id = eligible_id[k]
    #     if get_decision(donor, waiting_recipients[id], fm, u)
    #         accepting_id = id
    #         break
    #     end
    # end

    return accepting_id
end



@time allocate_one_donor!(waiting_recipients, donor, fm, u)



@inline function swapdeleteat!(v::Vector, i::Int)
    @inbounds begin
        v[i] = v[end]
        pop!(v)
    end
    return v
end



eligible_id = Int[]
scores = Float64[]
perm = Int[]

n = 0

@time for donor in new_donors

    # If you must drop expired each donor, still do it, but avoid extra work later
    filter!(r -> !is_expired(r, donor.arrival), waiting_recipients)

    id = allocate_one_donor!(waiting_recipients, donor, fm, u;
                            eligible_id=eligible_id, scores=scores, perm=perm)
    if id != 0
        # deleteat!(waiting_recipients, id)  # see next section for faster option
        swapdeleteat!(waiting_recipients, id)
    else
        n+=1
    end
end

n