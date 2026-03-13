

"""
    offered_recipients(df) -> AbstractDataFrame

For a single donor, return the score-ranked recipients that would be offered an
organ. Rows are kept up to the second transplanted recipient (`STATUS == "TX"`),
or up to the first if only one recipient is transplanted. If no recipient is
transplanted, all rows are returned.
"""
function offered_recipients(df::AbstractDataFrame)
    @assert "DON_ID" in names(df) "Missing column :DON_ID"
    @assert "STATUS" in names(df) "Missing column :STATUS"
    @assert "DON_CAN_SCORE" in names(df) "Missing column :DON_CAN_SCORE"

    nrow(df) == 0 && return df
    @assert all(==(df.DON_ID[1]), df.DON_ID) "All rows must correspond to the same :DON_ID"

    df_sort = sort(df, :DON_CAN_SCORE, rev=true)

    accpos = findall(isequal("TX"), df_sort.STATUS)
    if isempty(accpos)
        return df_sort
    elseif length(accpos) == 1
        return df_sort[1:accpos[1], :]
    else
        return df_sort[1:accpos[2], :]
    end
end

"""
    build_decision_dataset(donors_filepath::String, recipients_filepath::String) -> DataFrame

Construct a donor–recipient decision dataset from donor and recipient CSV files.

The function:
- loads donor and recipient data with `load_donor` and `load_recipient`,
- restricts donors to death years 2014–2019,
- removes rows with missing values required for KDRI computation and HLA mismatch calculation,
- computes one KDRI value per `DON_ID`,
- keeps only recipients who would have been offered a kidney using `offered_recipients`,
- matches each donor–recipient pair with recipient covariates,
- computes recipient age, waiting time, blood group, and HLA mismatch,
- assigns each donor to a training or validation set.

# Returns
A `DataFrame` with columns

`DON_AGE`, `KDRI`, `CAN_AGE`, `CAN_WAIT`, `CAN_BLOOD`, `MISMATCH`, `DECISION`, `LEARNING_SET`.

# Notes
The train/validation split is performed at the donor level.
"""
function build_decision_dataset(donors_filepath::String, recipients_filepath::String, cpra_filepath::String)

    df_donors = load_donor(donors_filepath)
    df_recipients = load_recipient(recipients_filepath)
    cpra_by_can_id = build_last_cpra_registry(cpra_filepath)

    filter!(row -> year(row.DON_DEATH_TM) in 2014:2019, df_donors)

    # columns needed for KDRI computation
    var_kdri = [:DON_AGE, :HEIGHT, :WEIGHT, :HYPERTENSION, :DIABETES, :DEATH, :CREATININE, :DCD]
    dropmissing!(df_donors, var_kdri)

    # columns needed for mismatch counting (donors)
    var_mismatch_don = [:DON_A1, :DON_A2, :DON_B1, :DON_B2, :DON_DR1, :DON_DR2]
    dropmissing!(df_donors, var_mismatch_don)

    # columns needed for mismatch counting (recipient)
    var_mismatch_rec = [:CAN_ID, :CAN_BTH_DT, :CAN_DIAL_DT, :CAN_BLOOD, :CAN_A1, :CAN_A2, :CAN_B1, :CAN_B2, :CAN_DR1, :CAN_DR2]
    dropmissing!(df_recipients, var_mismatch_rec)

    # Compute kdri for each donor
    kdri_by_don_id = Dict{Int,Float64}()
    G = groupby(df_donors, :DON_ID)

    for g in G
        r = first(g)
        age = r.DON_AGE
        height = r.HEIGHT
        weight = r.WEIGHT
        hypertension = r.HYPERTENSION == 1
        diabetes = r.DIABETES == 1
        cva = r.DEATH in (4, 16)
        creatinine = KidneyAllocation.creatinine_mgdl(r.CREATININE)
        dcd = r.DCD == 1

        kdri_by_don_id[r.DON_ID] =
            KidneyAllocation.evaluate_kdri(age, height, weight, hypertension, diabetes, cva, creatinine, dcd)
    end

    # Retrieve the recipients that have received an offer
    data = similar(df_donors, 0)
    for g in G
        append!(data, offered_recipients(g))
    end

    data.KDRI = getindex.(Ref(kdri_by_don_id), data.DON_ID)

    # Keeping only the lines in df_donors where we have recipient dat
    valid_can_ids = Set(df_recipients.CAN_ID)
    filter!(row -> row.CAN_ID in valid_can_ids, data)

    rec_idx_by_id = Dict(id => i for (i, id) in enumerate(df_recipients.CAN_ID))

    age = Int[]
    waittime = Float64[]
    blood = String[]
    n_mismatch = Int[]
    cpra = Int[]
    learning_set = String[]

    don_id = unique(data.DON_ID)
    n_train = round(Int, length(don_id) * 0.8)
    rng = MersenneTwister(12345)
    train_id = Set(shuffle(rng, don_id)[1:n_train])

    for r in eachrow(data)
        ind = rec_idx_by_id[r.CAN_ID]

        age_r = KidneyAllocation.years_between(Date(df_recipients.CAN_BTH_DT[ind]), Date(r.DON_DEATH_TM))
        waittime_r = KidneyAllocation.fractionalyears_between(Date(df_recipients.CAN_DIAL_DT[ind]), Date(r.DON_DEATH_TM))
        blood_r = df_recipients.CAN_BLOOD[ind]

        n = 0
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_A1, r.DON_A2)), HLA.((df_recipients.CAN_A1[ind], df_recipients.CAN_A2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_B1, r.DON_B2)), HLA.((df_recipients.CAN_B1[ind], df_recipients.CAN_B2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_DR1, r.DON_DR2)), HLA.((df_recipients.CAN_DR1[ind], df_recipients.CAN_DR2[ind])))

        push!(age, age_r)
        push!(waittime, waittime_r)
        push!(blood, blood_r)
        push!(n_mismatch, n)
        push!(cpra, get(cpra_by_can_id, r.CAN_ID, 0))
        push!(learning_set, r.DON_ID in train_id ? "train" : "validation")
    end

    data.CAN_AGE = age
    data.CAN_WAIT = waittime
    data.CAN_BLOOD = blood
    data.MISMATCH = n_mismatch
    data.CPRA = cpra
    data.DECISION = data.DECISION .== "Acceptation"
    data.LEARNING_SET = learning_set

    select!(data, [:DON_AGE, :KDRI, :CAN_AGE, :CAN_WAIT, :CAN_BLOOD, :MISMATCH, :CPRA, :DON_CAN_SCORE, :DECISION, :LEARNING_SET])

    return data
end
