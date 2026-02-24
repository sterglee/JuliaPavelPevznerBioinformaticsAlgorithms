"""
Returns the overlap alignment of strings v and w.
Match score = 1, Mismatch and Indels = -2.
"""
function overlap_alignment(v::AbstractString, w::AbstractString)
    n, m = length(v), length(w)

    # Initialize the arrays (1-based indexing: n+1 rows, m+1 columns)
    S = zeros(Int, n + 1, m + 1)
    backtrack = zeros(Int, n + 1, m + 1)

    # Initialize the max score and tracking indices
    max_score = -3 * (n + m)
    max_indices = (1, 1)

    # Fill in the Score and Backtrack arrays
    for i in 2:n+1
        for j in 2:m+1
            # Python logic: scores = [Diag, Up, Left]
            # Match = 1, Mismatch = -2, Indel = -2
            match_val = (v[i-1] == w[j-1]) ? 1 : -2

            # 1: Diagonal, 2: Up (gap in w), 3: Left (gap in v)
            scores = [
                S[i-1, j-1] + match_val,
                S[i-1, j] - 2,
                S[i, j-1] - 2
                ]

            val, idx = findmax(scores)
            S[i, j] = val
            backtrack[i, j] = idx

            # Overlap alignment specific:
            # Check for max score ONLY on the last row or last column
            if i == n + 1 || j == m + 1
                if S[i, j] > max_score
                    max_score = S[i, j]
                    max_indices = (i, j)
                end
            end
        end
    end

    curr_i, curr_j = max_indices

    # Initialize aligned strings as substrings up to max indices
    v_aligned = collect(v[1:curr_i-1])
    w_aligned = collect(w[1:curr_j-1])

    # Backtrack to the first row or column
    while curr_i > 1 && curr_j > 1
        dir = backtrack[curr_i, curr_j]

        if dir == 2 # Up (gap in w)
            curr_i -= 1
            insert!(w_aligned, curr_j, '-')
            elseif dir == 3 # Left (gap in v)
            curr_j -= 1
            insert!(v_aligned, curr_i, '-')
        else # Diagonal
            curr_i -= 1
            curr_j -= 1
        end
    end

    # Return result: Score, truncated V string, truncated W string
    # We slice from curr_i/curr_j to remove the "unaligned" heads
    return (
        string(max_score),
        join(v_aligned[curr_i:end]),
        join(w_aligned[curr_j:end])
        )
end

# --- File I/O ---

function main()
    input_file = "data/stepic_7c.txt"
    if !isfile(input_file)
        println("Input file not found.")
        return
    end

    lines = readlines(input_file)
    word1, word2 = strip(lines[1]), strip(lines[2])

    alignment = overlap_alignment(word1, word2)

    # Print to console
    println(join(alignment, "\n"))

    # Save to file
    mkpath("output")
    open("output/Assignment_07C.txt", "w") do f
        write(f, join(alignment, "\n"))
    end
end

main()
