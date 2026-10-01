using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates, JLD2, Random, Test

using KidneyAllocation

import KidneyAllocation: build_recipient_registry, build_donor_registry

recipients_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Candidates.csv"
donors_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/Donors.csv"
cpra_filepath = "/Users/jalbert/Documents/PackageDevelopment.nosync/kidney-research/kidney_research/KidneyResearch/data/CandidatesCPRA.csv"

recipient_dict = build_recipient_registry(recipients_filepath, cpra_filepath)

recipient_dict[3]

import KidneyAllocation: check_df_column_constant, check_df_columns

"""
    first_exit_date(df::AbstractDataFrame) -> Union{DateTime,Nothing}

Return the earliest update time associated with an exit outcome for a single
recipient, or `nothing` if no exit outcome is recorded.
"""
function first_exit_date(df::AbstractDataFrame)::Union{DateTime,Nothing}
    isempty(df) && return nothing

    check_df_column_constant(df, :CAN_ID)
    check_df_columns(df, :OUTCOME, :UPDATE_TM)

    exit_outcomes = ("X", "TX VIVANT", "DCD", "TX")
    earliest_exit_time = nothing

    for row in eachrow(df)
        outcome = uppercase(strip(string(row.OUTCOME)))

        if outcome ∈ exit_outcomes
            update_time = DateTime(row.UPDATE_TM)

            earliest_exit_time = isnothing(earliest_exit_time) ?
                update_time :
                min(earliest_exit_time, update_time)
        end
    end

    return earliest_exit_time
end

"""
    filter_outcomes(df::AbstractDataFrame) -> DataFrame

Return records from the first waiting-list episode of a single recipient.

Rows after the first exit outcome are removed, while the exit record is
retained. Exact duplicate retained rows are removed.

If every retained `UPDATE_TM` precedes the recorded `CAN_LISTING_DT`,
`CAN_LISTING_DT` is replaced by `CAN_DIAL_DT` as a proxy for the unavailable
original listing date.

`:UPDATE_TM` must contain `DateTime` values so that outcomes occurring on the
same calendar date can be ordered correctly.
"""
function filter_outcomes(df::AbstractDataFrame)::DataFrame
    isempty(df) && return DataFrame(df)

    check_df_column_constant(df, :CAN_ID, :CAN_LISTING_DT, :CAN_DIAL_DT)
    check_df_columns(df, :OUTCOME, :UPDATE_TM)

    update_type = eltype(skipmissing(df[!, :UPDATE_TM]))

    update_type <: DateTime || throw(ArgumentError("Column :UPDATE_TM must contain DateTime values; got $update_type"))

    exit_time = first_exit_date(df)

    if exit_time === nothing
        filtered_df = DataFrame(df)
    else
        filtered_df = filter(:UPDATE_TM => date -> date ≤ exit_time, df)

        listing_date = first(filtered_df.CAN_LISTING_DT)
        dialysis_date = convert(typeof(listing_date), first(filtered_df.CAN_DIAL_DT))

        if maximum(filtered_df.UPDATE_TM) < listing_date
            filtered_df.CAN_LISTING_DT .= dialysis_date
        end

    end

    # Removing duplicated rows (happens for recipient 501, 827 and 2051)
    unique!(filtered_df)

    return filtered_df
end







df = CSV.read("test/data/unfiltered_outcomes.csv", DataFrame)

g = groupby(df, :CAN_ID)

g[2]

first_exit_date(g[2])


df = CSV.read(recipients_filepath, DataFrame, missingstring=["", "NULL", "-"])


#501, 827 and 2051
df2 = filter(:CAN_ID => x->x==3, df)

filter_outcomes(df2)

