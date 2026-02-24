#!/usr/bin/env julia

"""
Compute the score of the middle column in linear-space alignment.
Returns:
middle_column_scores::Vector{Int}
backtrack::Vector{Int}
"""
function middle_column_score(
    v::String,
    w::String,
    scoring_matrix::Dict{Tuple{Char,Char},Int},
    sigma::Int
    )

    n = length(v)
    m = length(w)
    half = m ÷ 2

    # Two columns only (linear space)
    S = zeros(Int, n+1, 2)
    backtrack = zeros(Int, n+1)

    # Initialize first column
    for i in 0:n
        S[i+1, 1] = -i * sigma
    end

    # Fill columns up to middle
    for j in 1:half
        for i in 0:n
            if i == 0
                S[i+1, 2] = -j * sigma
            else
                scores = (
                    S[i, 1] + scoring_matrix[(v[i], w[j])],  # diagonal
                    S[i+1, 1] - sigma,                      # horizontal
                    S[i, 2] - sigma                         # vertical
                    )

                best = maximum(scores)
                S[i+1, 2] = best
                backtrack[i+1] = findfirst(==(best), scores) - 1
            end
        end

        # Shift column if not final
        if j != half
            S[:,1] .= S[:,2]
        end
    end

    return S[:,2], backtrack
end


"""
Return the middle edge of the alignment graph.
Returns:
((i, j), (i2, j2))
"""
function middle_edge(
    v::String,
    w::String,
    scoring_matrix::Dict{Tuple{Char,Char},Int},
    sigma::Int
    )

    n = length(v)
    m = length(w)
    half = m ÷ 2

    # Source → middle
    source_to_middle, _ =
        middle_column_score(v, w, scoring_matrix, sigma)

    # Middle → sink (reverse strings)
    v_rev = reverse(v)
    w_rev = reverse(w)

    middle_to_sink, backtrack =
        middle_column_score(v_rev, w_rev, scoring_matrix, sigma)

    middle_to_sink = reverse(middle_to_sink)
    backtrack = reverse(backtrack)

    # Component-wise sum
    scores = source_to_middle .+ middle_to_sink

    # Argmax
    max_middle = argmax(scores) - 1   # convert to 0-based

    # Determine next node
    if max_middle == n
        next_node = (max_middle, half + 1)
    else
        direction = backtrack[max_middle + 1]

        candidates = [
            (max_middle + 1, half + 1),  # down
            (max_middle,     half + 1),  # right
            (max_middle + 1, half)       # diag
            ]

        next_node = candidates[direction + 1]
    end

    return (max_middle, half), next_node
end

