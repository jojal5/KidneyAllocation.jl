

function retrieve_decision_data(donors_filepath::String, recipients_filepath::String)

    df_donors = load_donor(donors_filepath)
    df_recipients = load_recipient(recipients_filepath)

    filter!(row -> year(row.DON_DEATH_TM) ∈ 2014:2019, df_donors)

    # columns needed for KDRI computation
    VAR_KDRI = [:DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD]
    dropmissing!(df_donors, VAR_KDRI)

    # columns needed for mismatch counting (donors)
    VAR_MISMATCH_DON = [:DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2]
    dropmissing!(df_donors, VAR_MISMATCH_DON)

    # columns needed for mismatch counting (recipients)
    VAR_MISMATCH_REC = [:CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2]
    dropmissing!(df_recipients, VAR_MISMATCH_REC)

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
    # waittime = Union{Float64,Missing}[]
    waittime = Float64[]
    blood = String[]
    n_mismatch = Int64[]

    for r in eachrow(data)

        ind = findfirst(df_recipients.CAN_ID .== r.CAN_ID)

        age_r = KidneyAllocation.years_between(Date(df_recipients.CAN_BTH_DT[ind]), Date(r.DON_DEATH_TM))

        if ismissing(df_recipients.CAN_DIAL_DT[ind])
            waittime_r = missing
        else
            waittime_r = KidneyAllocation.fractionalyears_between(Date(df_recipients.CAN_DIAL_DT[ind]), Date(r.DON_DEATH_TM))
        end
        blood_r = df_recipients.CAN_BLOOD[ind]

        n = 0;
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_A1, r.DON_A2)), HLA.((df_recipients.CAN_A1[ind], df_recipients.CAN_A2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_B1, r.DON_B2)), HLA.((df_recipients.CAN_B1[ind], df_recipients.CAN_B2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_DR1, r.DON_DR2)), HLA.((df_recipients.CAN_DR1[ind], df_recipients.CAN_DR2[ind])))

        push!(age, age_r)
        push!(waittime, waittime_r)
        push!(blood, blood_r)
        push!(n_mismatch, n)

    end

    data.CAN_AGE = age
    data.CAN_WAIT = waittime
    data.CAN_BLOOD = blood
    data.MISMATCH = n_mismatch

    data.DECISION = data.DECISION .== "Acceptation"

    select!(data, [:DON_AGE, :KDRI, :CAN_AGE, :CAN_WAIT, :CAN_BLOOD, :MISMATCH, :DECISION])

    return data

end



# function get_decision(donor::Donor, recipient::Recipient, fm::StatsModels.TableRegressionModel, u::Real)

#     if is_abo_compatible(get_bloodtype(donor), get_bloodtype(recipient))

#         arrival = get_arrival(donor)
#         DON_AGE = donor.age
#         KDRI = donor.kdri

#         CAN_AGE = years_between(recipient.birth, arrival)
#         CAN_WAIT = fractionalyears_between(recipient.dialysis, arrival)
#         CAN_BLOOD = string(get_bloodtype(recipient))

#         df = DataFrame(DON_AGE=DON_AGE, KDRI=KDRI, CAN_AGE=CAN_AGE, CAN_WAIT=CAN_WAIT, CAN_BLOOD=CAN_BLOOD)

#         θ̂ = predict(fm, df)[]

#         decision = θ̂ > u 
#     else
#         decision = false
#     end

# end



function allocate_one_donor(
    donor::Donor,
    recipients::Vector{Recipient},
    dm::AbstractDecisionModel,
    is_unallocated::BitVector=trues(length(recipients))
)

    eligible_indices = get_eligible_recipient_indices(donor, recipients, is_unallocated)
    ranked_indices = rank_eligible_indices_by_score(donor, recipients, eligible_indices)

    p = acceptance_probability(dm, recipients[ranked_indices], donor)

    acceptation = p .> dm.threshold

    if any(acceptation)
        ind = findfirst(acceptation)
    else
        ind = 0
    end

    return ind

end


"""
    get_eligible_recipient_indices(donor, recipients, is_unallocated) 

Return the indices of recipients eligible to receive an offer from `donor`.

A recipient is considered eligible if it:
- is currently unallocated,
- is active at the donor arrival date, and
- is ABO-compatible with the donor.

# Arguments
- `donor::Donor`: Donor being allocated.
- `recipients::Vector{Recipient}`: Current waiting list.
- `is_unallocated::BitVector`: Optional mask indicating which recipients are
  still available for allocation (default: all `true`).

# Returns
- `Vector{Int}`: Indices into `recipients` identifying eligible recipients.
"""
function get_eligible_recipient_indices(
    donor::Donor,
    recipients::Vector{Recipient},
    is_unallocated::BitVector = trues(length(recipients)),
)

    arrival = donor.arrival

    eligible_mask = copy(is_unallocated)
    eligible_mask .&= is_active.(recipients, arrival)
    eligible_mask .&= is_abo_compatible.(donor, recipients)

    return findall(eligible_mask)
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



# C'est important pour cette fonction qu'on envoie les patients en attente, pas tous les patients de la banque.
function allocate(donors::Vector{Donor}, recipients::Vector{Recipient}, dm::AbstractDecisionModel; until::Int64=-9999)

    is_unallocated = trues(length(recipients))                 
    allocated_recipient_index = zeros(Int64,length(donors))

    for donor_idx in eachindex(donors)

        donor = donors[donor_idx]

        allocated_recipient_index[donor_idx] = allocate_one_donor(donor, recipients, dm, is_unallocated)

        if allocated_recipient_index[donor_idx] != 0
            is_unallocated[allocated_recipient_index[donor_idx]] = false
        end

        if allocated_recipient_index[donor_idx] == until
            break
        end

    end

    return allocated_recipient_index
end