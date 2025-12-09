# KidneyAllocation.jl

[![Build Status](https://github.com/jojal5/KidneyAllocation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jojal5/KidneyAllocation.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for representing kidney donors and recipients, computing key transplant indices (EPTS and KDRI), and evaluating allocation compatibility and score.

[!NOTE]  
**All donor and recipient examples in this documentation use synthetic data.  
No real clinical or confidential information is included.**

---

## Overview

KidneyAllocation.jl provides:

- Simple, validated structures for **Donor** and **Recipient**
- HLA typing using compact `HLA` integers
- ABO compatibility checks
- **EPTS** (recipient survival index)
- **KDRI** (donor risk index)
- A combined **allocation score** for donor–recipient matching used by Transplant Québec

---

## 🫘 Donor Structure

```julia
struct Donor <: TransplantEntity
    arrival::DateTime
    age::Int64
    blood::ABOGroup
    a1::HLA; a2::HLA
    b1::HLA; b2::HLA
    dr1::HLA; dr2::HLA
    kdri::Float64
end
```

Example:

```julia
donor = Donor(DateTime(2025,1,1), 45, A, 24, 26, 44, 51, 1, 4, 1.27)
```

---

## 👤 Recipient Structure

```julia
struct Recipient <: TransplantEntity
    birth::DateTime
    dialysis::DateTime
    arrival::DateTime
    blood::ABOGroup
    a1::HLA; a2::HLA
    b1::HLA; b2::HLA
    dr1::HLA; dr2::HLA
    CPRA::Int64
    expiration_date::Union{DateTime,Nothing}
end
```

Example:

```julia
recipient = Recipient(DateTime(1985,5,1), DateTime(2015,1,1),
                      DateTime(2024,1,1), A, 24,26, 44,51, 1,4, 80)
```

---

## 📈 EPTS

Compute Estimated Post-Transplant Survival:

```julia
epts = compute_epts(recipient)
```

---

## ⚖️ KDRI

Compute Kidney Donor Risk Index:

```julia
kdri = evaluate_kdri(age=45, height=175, weight=80,
                     ethnicity=6, hypertension=true, diabetes=false,
                     death=4, creatinine=120, hcv=false, dcd=false)
```

---

## 🧮 Allocation Score

A donor–recipient match score:

```julia
score(donor, recipient)
```

Includes components such as wait time, DR mismatch, CPRA, age gap, etc.

---

## 🔧 Example

```julia
if is_abo_compatible(donor.blood, recipient.blood)
    println("Score = ", score(donor, recipient))
else
    println("ABO incompatible")
end
```

