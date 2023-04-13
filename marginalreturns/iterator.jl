function _iterate(state::Vector{Int}, max_vals)
    state = copy(state)
    state[end] += 1
    i = length(state)
    while i > 1 && state[i] > max_vals[i]  # carry over!
        state[i-1] += 1
        i -= 1
    end
    state[1] > max_val[1] && return nothing
    state[i+1:end] .= state[i]
    return state
end

struct LexiIter
    inner
    rep::Int
end

import Base.iterate

function iterate(iter::LexiIter)
    res = iterate(iter.inner)
    isnothing(res) && return nothing
    inner_item, inner_state = res
    item = ntuple(_ -> inner_item, iter.rep)
    state = ntuple(_ -> 1, iter.rep)
    return item, state
end

function iterate(iter::LexiIter, state::Tuple)
    next_state = _iterate_state(iter, state)
    isnothing(next_state) && return nothing
    next_item = map(s -> first(iterate(iter.inner, s-1)), next_state)
    return next_item, next_state
end

function _iterate_state(iter::LexiIter, state::Tuple)
    inner_states = [s for s in state]
    res = iterate(iter.inner, inner_states[end])
    i = iter.rep
    while i > 1 && isnothing(res)  # find smallest index which doesn't return nothing on advancing
        i -= 1
        res = iterate(iter.inner, inner_states[i])
    end
    i==1 && isnothing(res) && return nothing
    _, s = res
    while i <= iter.rep
        inner_states[i] = s
        i += 1
    end
    return tuple(inner_states...)
end