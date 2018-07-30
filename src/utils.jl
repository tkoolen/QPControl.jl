"""
    only(x)

Return the first and only element of `x`. Throws an ArgumentError if
the length of x is not exactly 1.
"""
function only(x)
    length(x) == 1 || throw(ArgumentError("Expected a collection with exactly one element. Got length(x) = $(length(x)) instead."))
    first(x)
end
