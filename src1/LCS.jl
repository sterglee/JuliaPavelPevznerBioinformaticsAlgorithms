#!/usr/bin/env julia

# -------------------------------
# Longest Common Subsequence Problem
# -------------------------------

function longest_common_subsequence(v::AbstractString, w::AbstractString)::String
    n, m = length(v), length(w)
    # Initialize the dynamic programming matrix
    S = zeros(Int, n+1, m+1)

    # Fill DP table
    for i in 1:n
        for j in 1:m
            if v[i] == w[j]
                S[i+1, j+1] = S[i, j] + 1
            else
                S[i+1, j+1] = max(S[i+1, j], S[i, j+1])
            end
        end
    end

    # Recover the LCS from the DP table
    lcs = IOBuffer()
    i, j = n, m
    while i > 0 && j > 0
        if S[i+1, j+1] == S[i, j+1]
            i -= 1
            elseif S[i+1, j+1] == S[i+1, j]
            j -= 1
        else
            print(lcs, v[i])
            i -= 1
            j -= 1
        end
    end

    # Reverse the result since we built it backwards
    return String(reverse(take!(lcs)))
end

# -------------------------------
# Main script
# -------------------------------
function main()
    # Read input
    lines = readlines("data/stepic_6c.txt")
    dna1, dna2 = lines[1], lines[2]

    # Compute LCS
    lcs = longest_common_subsequence(dna1, dna2)

    # Print and save output
    println(lcs)

end

# Run main
main()
