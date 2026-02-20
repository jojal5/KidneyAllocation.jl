

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
    waittime = Float64[]
    blood = String[]
    n_mismatch = Int64[]
    learning_set = String[]

    don_id = unique(data.DON_ID)
    n_train = round(Int64, length(don_id)*.8)
    train_id = rand(MersenneTwister(12345), don_id, n_train)

    for r in eachrow(data)

        ind = findfirst(df_recipients.CAN_ID .== r.CAN_ID)

        age_r = KidneyAllocation.years_between(Date(df_recipients.CAN_BTH_DT[ind]), Date(r.DON_DEATH_TM))

        waittime_r = KidneyAllocation.fractionalyears_between(Date(df_recipients.CAN_DIAL_DT[ind]), Date(r.DON_DEATH_TM))

        blood_r = df_recipients.CAN_BLOOD[ind]

        n = 0;
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_A1, r.DON_A2)), HLA.((df_recipients.CAN_A1[ind], df_recipients.CAN_A2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_B1, r.DON_B2)), HLA.((df_recipients.CAN_B1[ind], df_recipients.CAN_B2[ind])))
        n += KidneyAllocation.mismatch_locus(HLA.((r.DON_DR1, r.DON_DR2)), HLA.((df_recipients.CAN_DR1[ind], df_recipients.CAN_DR2[ind])))

        push!(age, age_r)
        push!(waittime, waittime_r)
        push!(blood, blood_r)
        push!(n_mismatch, n)

        if r.DON_ID ∈ train_id
            push!(learning_set, "train")
        else
            push!(learning_set, "validation")
        end

    end

    data.CAN_AGE = age
    data.CAN_WAIT = waittime
    data.CAN_BLOOD = blood
    data.MISMATCH = n_mismatch

    data.DECISION = data.DECISION .== "Acceptation"

    data.LEARNING_SET = learning_set

    select!(data, [:DON_AGE, :KDRI, :CAN_AGE, :CAN_WAIT, :CAN_BLOOD, :MISMATCH, :DECISION, :LEARNING_SET])

    return data

end

