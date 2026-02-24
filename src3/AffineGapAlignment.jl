#!/usr/bin/env julia

"""
Global Alignment with Affine Gap Penalty

Returns:
max_score::Int
v_aligned::String
w_aligned::String
"""
function global_alignment_affine_gap_penalty(
    v::String,
    w::String,
    scoring_matrix::Dict{Tuple{Char,Char},Int},
    sigma::Int,
    epsilon::Int
    )

    n = length(v)
    m = length(w)

    # Score matrices
    S_lower  = zeros(Int, n+1, m+1)
    S_middle = zeros(Int, n+1, m+1)
    S_upper  = zeros(Int, n+1, m+1)
    backtrack = zeros(Int, n+1, m+1)

    # Initialize first column
    for i in 2:n+1
        S_lower[i,1]  = -sigma - (i-2)*epsilon
        S_middle[i,1] = -sigma - (i-2)*epsilon
        S_upper[i,1]  = -10*sigma   # negative infinity proxy
    end

    # Initialize first row
    for j in 2:m+1
        S_upper[1,j]  = -sigma - (j-2)*epsilon
        S_middle[1,j] = -sigma - (j-2)*epsilon
        S_lower[1,j]  = -10*sigma
    end

    # Fill matrices
    for i in 2:n+1
        for j in 2:m+1

            S_lower[i,j] = max(
                S_lower[i-1,j] - epsilon,
                S_middle[i-1,j] - sigma
                )

            S_upper[i,j] = max(
                S_upper[i,j-1] - epsilon,
                S_middle[i,j-1] - sigma
                )

            match_score = scoring_matrix[(v[i-1], w[j-1])]

            middle_scores = (
                S_lower[i,j],
                S_middle[i-1,j-1] + match_score,
                S_upper[i,j]
                )

            S_middle[i,j] = maximum(middle_scores)

            # backtrack: 1=lower, 2=diag, 3=upper
            backtrack[i,j] = findfirst(==(S_middle[i,j]), middle_scores)
        end
    end

    # Backtracking
    i, j = n+1, m+1
    max_score = S_middle[i,j]

    v_aligned = collect(v)
    w_aligned = collect(w)

    function insert_indel!(word::Vector{Char}, pos::Int)
        insert!(word, pos, '-')
    end

    while i > 1 && j > 1
        if backtrack[i,j] == 1
            i -= 1
            insert_indel!(w_aligned, j)
            elseif backtrack[i,j] == 3
            j -= 1
            insert_indel!(v_aligned, i)
        else
            i -= 1
            j -= 1
        end
    end

    while i > 1
        i -= 1
        insert_indel!(w_aligned, 1)
    end

    while j > 1
        j -= 1
        insert_indel!(v_aligned, 1)
    end

    return max_score, String(v_aligned), String(w_aligned)
end

