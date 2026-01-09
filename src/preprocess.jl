"""
    parse_hla_int(s::String)

Parse HLA string value to Int64 by extracting only the numerical expression.

## Details

Parse HLA-like values such as "24", "24L", "24Low" -> 24
"""
function parse_hla_int(s::String)
    m = match(r"^\s*(\d+)", s)
    return parse(Int64, m.captures[1])
end