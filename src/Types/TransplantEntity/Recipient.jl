
"""
    Recipient(...) <: TransplantEntity

Represents a recipient in the kidney transplantation system (internal use).

# Fields
- `birth::Date`: Birth date of the recipient.
- `dialysis::Date`: Start date of dialysis.
- `arrival::Date`: Date the recipient was added to the waitlist.
- `blood::ABOGroup`: ABO blood type of the recipient.
- `a1::HLA`, `a2::HLA`: HLA-A antigens.
- `b1::HLA`, `b2::HLA`: HLA-B antigens.
- `dr1::HLA`, `dr2::HLA`: HLA-DR antigens.
- `cpra::Int64`: Calculated Panel Reactive Antibody (0–100).
- `expiration_date::Union{Date,Nothing}`: Eligibility expiration date, or `nothing`.
- `active_waiting_proportion::Float64`: Fraction of the observed waiting period spent active, in `[0, 1]`.
"""
struct Recipient <: TransplantEntity
    birth::Date
    dialysis::Date
    arrival::Date
    blood::ABOGroup
    a1::HLA
    a2::HLA
    b1::HLA
    b2::HLA
    dr1::HLA
    dr2::HLA
    cpra::Int64
    expiration_date::Union{Date,Nothing}
    active_waiting_proportion::Float64

    function Recipient(
        birth::Date,
        dialysis::Date,
        arrival::Date,
        blood::ABOGroup,
        a1::HLA, a2::HLA,
        b1::HLA, b2::HLA,
        dr1::HLA, dr2::HLA,
        cpra::Int64;
        expiration_date::Union{Date,Nothing}=nothing,
        active_waiting_proportion::Real=1.0)

        a1 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a1 = $a1"))
        a2 ∈ VALID_HLA_A || throw(ArgumentError("Invalid A allele a2 = $a2"))
        b1 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b1 = $b1"))
        b2 ∈ VALID_HLA_B || throw(ArgumentError("Invalid B allele b2 = $b2"))
        dr1 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr1 = $dr1"))
        dr2 ∈ VALID_HLA_DR || throw(ArgumentError("Invalid DR allele dr2 = $dr2"))

        (0 <= cpra <= 100) || throw(ArgumentError("cpra must be in [0, 100], got $cpra"))

        proportion = Float64(active_waiting_proportion)

        isfinite(proportion) ||
            throw(ArgumentError(
                "`active_waiting_proportion` must be finite, got $active_waiting_proportion",
            ))

        (0.0 <= proportion <= 1.0) ||
            throw(ArgumentError(
                "`active_waiting_proportion` must be in [0, 1], got $active_waiting_proportion",
            ))

        return new(
            birth, dialysis, arrival, blood,
            a1, a2, b1, b2, dr1, dr2,
            cpra, expiration_date, proportion,
        )
    end
end

# Outer constructors (mixed Date/DateTime)

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::HLA, a2::HLA,
    b1::HLA, b2::HLA,
    dr1::HLA, dr2::HLA,
    cpra::Integer;
    expiration_date::Union{Date,DateTime,Nothing}=nothing,
    active_waiting_proportion::Real=1.)

    return Recipient(Date(birth), Date(dialysis), Date(arrival),
        blood,
        a1, a2, b1, b2, dr1, dr2,
        Int64(cpra);
        expiration_date=expiration_date === nothing ? nothing : Date(expiration_date),
        active_waiting_proportion=active_waiting_proportion)
end

function Recipient(birth::Union{Date,DateTime},
    dialysis::Union{Date,DateTime},
    arrival::Union{Date,DateTime},
    blood::ABOGroup,
    a1::Integer, a2::Integer,
    b1::Integer, b2::Integer,
    dr1::Integer, dr2::Integer,
    cpra::Integer;
    expiration_date::Union{Date,DateTime,Nothing}=nothing,
    active_waiting_proportion::Real=1.)

    return Recipient(Date(birth), Date(dialysis), Date(arrival),
        blood,
        HLA(a1), HLA(a2),
        HLA(b1), HLA(b2),
        HLA(dr1), HLA(dr2),
        Int64(cpra);
        expiration_date=expiration_date === nothing ? nothing : Date(expiration_date),
        active_waiting_proportion=active_waiting_proportion)
end


"""
    get_birth(r::Recipient)

Get birth date of recipient `r`.
"""
function get_birth(r::Recipient)
    return r.birth
end

"""
    get_dialysis(r::Recipient)

Get dialysis start date of recipient `r`.
"""
function get_dialysis(r::Recipient)
    return r.dialysis
end

"""
    has_expiration(r::Recipient)::Bool

Returns `true` if the recipient has a defined expiration date,
and `false` if `expiration_date` is `nothing`.
"""
has_expiration(r::Recipient) = r.expiration_date !== nothing

"""
    is_registered(r::Recipient, t::Date) -> Bool

Return whether recipient `r` is registered on the waiting list at date `t`.

A recipient is registered from `r.arrival` through `r.expiration_date`
(inclusively), or indefinitely when `r.expiration_date` is `nothing`.
"""
is_registered(r::Recipient, t::Date) =
    t >= r.arrival && !is_expired(r, t)

"""
    is_active(r::Recipient, t::Date; rng=Random.default_rng()) -> Bool

Return whether recipient `r` is active at date `t`.

An unregistered recipient is never active. For a registered recipient, activity
is sampled independently with probability `r.active_waiting_proportion`.
"""
function is_active(
    r::Recipient,
    t::Date;
    rng::AbstractRNG=Random.default_rng(),
)::Bool
    return is_registered(r, t) &&
           rand(rng) < r.active_waiting_proportion
end

"""
    is_expired(r::Recipient, t::Date) -> Bool

Returns `true` if the recipient has an expiration date and it is
strictly earlier than the given date `t`. Returns `false` if
there is no expiration date.
"""
is_expired(r::Recipient, t::Date) =
    has_expiration(r) && r.expiration_date < t


"""
    shift_recipient_timeline(recipient::Recipient, new_arrival::Date) -> Recipient

Shift the recipient's `birth`, `dialysis`, and `expiration_date` so that the
recipient's timeline is consistent with a new `arrival` date.

## Details

All dates are shifted by the same number of days: the difference between the
current `recipient.arrival` and the new `arrival`. This preserves age and
waiting-time durations relative to the (new) arrival date.
"""
function shift_recipient_timeline(recipient::Recipient, new_arrival::Date)::Recipient
    # Signed shift (in days) from old arrival to new arrival
    shift_days = days_between(recipient.arrival, new_arrival)

    birth = recipient.birth + Day(shift_days)
    dialysis = recipient.dialysis + Day(shift_days)

    expiration_date = recipient.expiration_date === nothing ? nothing :
                      recipient.expiration_date + Day(shift_days)

    return Recipient(birth,
        dialysis,
        new_arrival,
        recipient.blood,
        recipient.a1, recipient.a2,
        recipient.b1, recipient.b2,
        recipient.dr1, recipient.dr2,
        recipient.cpra;
        expiration_date=expiration_date,
        active_waiting_proportion=recipient.active_waiting_proportion)
end

"""
    sim_cpra_compatibility(recipient) -> Bool

Simulate the recipient compatibility based on CPRA for any given donor. Return `true` if the recipient is compatible.

Compatibility is simulated as a Bernoulli trial with probability
`1 - recipient.cpra/100`.
"""
function sim_cpra_compatibility(recipient::Recipient)
    return rand() ≥ recipient.cpra / 100
end


function Base.show(io::IO, ::MIME"text/plain", r::Recipient)
    print(io,
        "Recipient\n",
        "  Birth Date        : $(r.birth)\n",
        "  Dialysis Start    : $(r.dialysis)\n",
        "  Arrival Date      : $(r.arrival)\n",
        "  Blood Type        : $(r.blood)\n",
        "  HLA-A             : $(r.a1), $(r.a2)\n",
        "  HLA-B             : $(r.b1), $(r.b2)\n",
        "  HLA-DR            : $(r.dr1), $(r.dr2)\n",
        "  CPRA              : $(r.cpra)\n",
        "  Expiration        : ",
        r.expiration_date === nothing ? "none\n" : string(r.expiration_date)*"\n",
        "  Active proportion : ",
        string(round(r.active_waiting_proportion, digits=2))
    )
end

function Base.summary(io::IO, r::Recipient)
    print(io,
        "Recipient(blood=$(r.blood), CPRA=$(r.cpra), ",
        "A=($(r.a1),$(r.a2)), ",
        "B=($(r.b1),$(r.b2)), ",
        "DR=($(r.dr1),$(r.dr2)))"
    )
end


function Base.show(io::IO, ::MIME"text/plain", recipients::AbstractVector{<:Recipient})
    n = length(recipients)
    println(io, "$n-element Vector{Recipient}:")

    n == 0 && return

    limit = get(io, :limit, false)
    rows, _ = displaysize(io)
    max_lines = max(rows - 1, 1)  # available lines after the header

    if !limit || n <= max_lines
        # Show everything
        for r in recipients
            print(io, " ")
            summary(io, r)
            println(io)
        end
        return
    end

    # Long vector display: show head, ellipsis, tail
    avail = max_lines - 1          # reserve one line for " ⋮"
    head = max(1, avail ÷ 2)
    tail = max(1, avail - head)

    if head + tail >= n
        head = min(head, n)
        tail = 0
    end

    for i in 1:head
        print(io, " ")
        summary(io, recipients[i])
        println(io)
    end

    println(io, " ⋮")

    for i in (n-tail+1):n
        print(io, " ")
        summary(io, recipients[i])
        println(io)
    end
end