

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