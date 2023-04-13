"""
General notes: The aim is to implement a Julia iterator according to the following documentation.
https://docs.julialang.org/en/v1/manual/interfaces/

TODO:
* Implement a length function. We should be able to express the length as a formula, but I haven't thought about it yet... This is important because then the ProgressBar works again.
* Implement any other optional methods?
"""

struct LexiIter
    inner  # the 'inner' iterator that is repeated
    rep::Int  # the number of repetitions
end

import Base.iterate

function iterate(iter::LexiIter)
    res = iterate(iter.inner)
    isnothing(res) && return nothing
    inner_item, inner_state = res
    item = ntuple(_ -> inner_item, iter.rep)
    state = ntuple(_ -> (inner_item, inner_state), iter.rep)
    return item, state
end

function iterate(iter::LexiIter, state::Tuple)
    next_state = _iterate_state(iter, state)
    isnothing(next_state) && return nothing
    next_item = map(s -> first(s), next_state)
    return next_item, next_state
end

function _iterate_state(iter::LexiIter, state::Tuple)
    next_state = [s for s in state]
    i = iter.rep
    res = iterate(iter.inner, last(state[i]))
    while i > 1 && isnothing(res)  # find largest index which doesn't return nothing on advancing
        i -= 1
        res = iterate(iter.inner, last(state[i]))
    end
    i==1 && isnothing(res) && return nothing
    for k in i:iter.rep; next_state[k] = res; end
    return tuple(next_state...)
end