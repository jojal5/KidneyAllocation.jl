

function retrieve_decision_data(donors_filepath::String, recipients_filepath::String)

    df_donors = load_donor(donors_filepath)
    df_recipients = load_recipient(recipients_filepath)

    filter!(row -> year(row.DON_DEATH_TM) ∈ 2014:2019, df_donors)

    # columns needed for KDRI computation
    VAR_KDRI = [:DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD]

    dropmissing!(df_donors, VAR_KDRI)

    kdri = Float64[]

    for r in eachrow(df_donors)

        age = r.DON_AGE
        height = r.HEIGHT
        weight = r.WEIGHT
        hypertension = r.HYPERTENSION == 1
        diabetes = r.DIABETES == 1
        cva = (r.DEATH == 4) || (r.DEATH == 16)
        creatinine = KidneyAllocation.creatinine_mgdl(r.CREATININE)
        dcd = r.DCD == 1 # TODO À VÉRIFIER si c'est bien 1, sinon c'est 2 (Anastasiya a confirmé le code)

        kdri_r = KidneyAllocation.evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)

        push!(kdri, kdri_r)
    end

    df_donors.KDRI = kdri

    data = deepcopy(df_donors)

    # Keeping only the lines in df_donors where we have recipient data
    filter!(row -> row.CAN_ID ∈ unique(df_recipients.CAN_ID), data)

    age = Int64[]
    waittime = Union{Float64,Missing}[]
    blood = String[]

    for r in eachrow(data)

        ind = findfirst(df_recipients.CAN_ID .== r.CAN_ID)

        age_r = KidneyAllocation.years_between(Date(df_recipients.CAN_BTH_DT[ind]), Date(r.DON_DEATH_TM))

        if ismissing(df_recipients.CAN_DIAL_DT[ind])
            waittime_r = missing
        else
            waittime_r = KidneyAllocation.fractionalyears_between(Date(df_recipients.CAN_DIAL_DT[ind]), Date(r.DON_DEATH_TM))
        end
        blood_r = df_recipients.CAN_BLOOD[ind]

        push!(age, age_r)
        push!(waittime, waittime_r)
        push!(blood, blood_r)

    end

    data.CAN_AGE = age
    data.CAN_WAIT = waittime
    data.CAN_BLOOD = blood

    data.DECISION = data.DECISION .== "Acceptation"

    select!(data, [:DON_AGE, :KDRI, :CAN_AGE, :CAN_WAIT, :CAN_BLOOD, :DECISION])

    return data

end


function fit_decision_threshold(fm::StatsModels.TableRegressionModel)

    gt = Int64.(response(fm))

    θ̂ = predict(fm)

    fobj(u::Real) = -f1score(roc(gt, θ̂, u))

    res = optimize(fobj, 0.01, 0.75)

    u = res.minimizer

    return u

end

function get_decision(donor::Donor, recipient::Recipient, fm::StatsModels.TableRegressionModel, u::Real)

    if is_abo_compatible(get_bloodtype(donor), get_bloodtype(recipient))

        arrival = get_arrival(donor)
        DON_AGE = donor.age
        KDRI = donor.kdri

        CAN_AGE = years_between(recipient.birth, arrival)
        CAN_WAIT = fractionalyears_between(recipient.dialysis, arrival)
        CAN_BLOOD = string(get_bloodtype(recipient))

        df = DataFrame(DON_AGE=DON_AGE, KDRI=KDRI, CAN_AGE=CAN_AGE, CAN_WAIT=CAN_WAIT, CAN_BLOOD=CAN_BLOOD)

        θ̂ = predict(fm, df)[]

        decision = θ̂ > u 
    else
        decision = false
    end

end


# function allocate_one_donor(
#     donor::Donor,
#     eligible_recipients::Vector{Recipient},
#     fm,
#     u)::Int64

#     ind = rank_recipients(donor, eligible_recipients)

#     accepting_index = 0

#     for i in eachindex(ind)
#         if get_decision(donor, eligible_recipients[ind[i]], fm, u)
#             accepting_index = ind[i]
#             break
#         end
#     end

#     return accepting_index
# end

"""
    allocate_one_donor(donor, recipients, ranked_indices, fm, u) -> Int

Evaluate sequential offers for `donor` to `recipients[ranked_indices]` and
return the first accepting recipient index (into `recipients`). Returns `0`
if no recipient accepts.
"""
function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    ranked_indices::AbstractVector{<:Integer},
    fm,
    u,
)

    for i in eachindex(ranked_indices)
        r_idx = ranked_indices[i]
        if get_decision(donor, recipients[r_idx], fm, u)
            return Int(r_idx)
        end
    end

    return 0
end


"""
    rank_eligible_indices_by_score(donor, recipients, eligible_indices)

Return a new vector of recipient indices ranked by decreasing attribution score
for `donor`.

The returned vector contains the same elements as `eligible_indices`, reordered
so that `score(donor, recipients[i])` is decreasing.

# Arguments
- `donor::Donor`: Donor being allocated.
- `recipients::Vector{Recipient}`: Current waiting list.
- `eligible_indices::AbstractVector{<:Integer}`: Indices into `recipients`
  identifying eligible recipients.
"""
function rank_eligible_indices_by_score(
    donor::Donor,
    recipients::Vector{Recipient},
    eligible_indices::AbstractVector{<:Integer},
)
    scores = score.(Ref(donor), recipients[eligible_indices])
    p = sortperm(scores; rev=true)
    return eligible_indices[p]
end




# function allocate(donors::Vector{Donor}, recipients::Vector{Recipient}, fm, u; until::Int64=-9999)

#     is_unallocated = trues(length(recipients))                 
#     allocated_recipient_index = zeros(Int64,length(donors))

#     eligible_mask = falses(length(recipients))
#     eligible_indices = Int64[]

#     for donor_idx in eachindex(donors)

#         donor = donors[donor_idx]
#         arrival = donor.arrival

#         # Build eligibility mask with minimal temporaries
#         eligible_mask .= is_unallocated                         # start from availability
#         eligible_mask .&= is_active.(recipients, arrival)        # active at arrival
#         eligible_mask .&= is_abo_compatible.(donor, recipients)  # ABO compatible

#         # Reuse eligible_indices instead of allocating a new vector each time
#         empty!(eligible_indices)
#         append!(eligible_indices, findall(eligible_mask))

#         if isempty(eligible_indices)
#             allocated_recipient_index[donor_idx] = 0
#             continue
#         end

#         chosen_recipient = allocate_one_donor(recipients[eligible_mask], donor, fm, u)

#         if chosen_recipient != 0
#             is_unallocated[eligible_indices[chosen_recipient]] = false
#             allocated_recipient_index[donor_idx] = eligible_indices[chosen_recipient]
#             if eligible_indices[chosen_recipient] == until
#                 break
#             end
#         else
#             allocated_recipient_index[donor_idx] = 0
#         end

#     end

#     return allocated_recipient_index
# end


function allocate(donors::Vector{Donor}, recipients::Vector{Recipient}, fm, u; until::Int64=-9999)

    is_unallocated = trues(length(recipients))                 
    allocated_recipient_index = zeros(Int64,length(donors))

    eligible_mask = falses(length(recipients))
    eligible_indices = Int64[]

    for donor_idx in eachindex(donors)

        donor = donors[donor_idx]
        arrival = donor.arrival

        # Build eligibility mask with minimal temporaries
        eligible_mask .= is_unallocated                         # start from availability
        eligible_mask .&= is_active.(recipients, arrival)        # active at arrival
        eligible_mask .&= is_abo_compatible.(donor, recipients)  # ABO compatible

        # Reuse eligible_indices instead of allocating a new vector each time
        empty!(eligible_indices)
        append!(eligible_indices, findall(eligible_mask))

        if isempty(eligible_indices)
            allocated_recipient_index[donor_idx] = 0
            continue
        end

        ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

        accepted_index = allocate_one_donor(donor, recipients, ranked_indices, fm, u)

        if accepted_index != 0
            is_unallocated[ranked_indices[accepted_index]] = false
            allocated_recipient_index[donor_idx] = ranked_indices[accepted_index]
            if ranked_indices[accepted_index] == until
                break
            end
        else
            allocated_recipient_index[donor_idx] = 0
        end

    end

    return allocated_recipient_index
end