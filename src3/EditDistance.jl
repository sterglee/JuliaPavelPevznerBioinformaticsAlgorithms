"""
Returns the edit distance of strings v and w using dynamic programming.
Based on Assignment 07A from Bioinformatics Algorithms.
"""
function edit_distance(v::AbstractString, w::AbstractString)
    n, m = length(v), length(w)

    # Initialize matrix M (n+1 x m+1)
    # Julia uses 1-based indexing, so M[i, j] maps to M[row, col]
    M = zeros(Int, n + 1, m + 1)

    # Base cases: distance from empty string
    for i in 1:n+1
        M[i, 1] = i - 1
    end
    for j in 1:m+1
        M[1, j] = j - 1
    end

    # Compute each entry of M
    for i in 2:n+1
        for j in 2:m+1
            if v[i-1] == w[j-1]
                # Match: no increase in distance
                M[i, j] = M[i-1, j-1]
            else
                # Mismatch/Insertion/Deletion: take the minimum of the three neighbors + 1
                M[i, j] = min(
                    M[i-1, j] + 1,   # Deletion
                    M[i, j-1] + 1,   # Insertion
                    M[i-1, j-1] + 1  # Substitution
                    )
            end
        end
    end

    # Return the bottom-right entry
    return M[n+1, m+1]
end

# --- Main Execution ---

function main()
    input_file = "data/EditDistance.txt"
    output_file = "output/EditDistance.txt"

    if isfile(input_file)
        # Read the input data
        lines = readlines(input_file)
        word1 = strip(lines[1])
        word2 = strip(lines[2])

        # Get the edit distance
        e_dist = edit_distance(word1, word2)

        # Print result
        println(e_dist)

        # Save result
        mkpath("output")
        write(output_file, string(e_dist))
    else
        println("Input file not found at: ", input_file)
    end
end

# Execute the script
 main()
