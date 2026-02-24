#!/usr/bin/env julia

# Manhattan Tourist Problem
# Computes the longest path from (0,0) to (n,m) using a grid with weighted down and right moves.

using DelimitedFiles

function manhattan_tourist(n::Int, m::Int, down::AbstractArray{<:Integer,2}, right::AbstractArray{<:Integer,2})
    # Initialize the score matrix
    S = zeros(Int, n+1, m+1)

    # Fill first column
    for i in 1:n
        S[i+1, 1] = S[i, 1] + down[i,1]
    end

    # Fill first row
    for j in 1:m
        S[1, j+1] = S[1, j] + right[1,j]
    end

    # Fill the rest of the matrix
    for i in 1:n
        for j in 1:m
            S[i+1, j+1] = max(S[i, j+1] + down[i, j+1], S[i+1, j] + right[i+1, j])
        end
    end

    return S[n+1, m+1]
end

# -------------------------------
# Main script
# -------------------------------
function main()
    # Read input data
    lines = readlines("data/stepic_6b.txt")

    n = parse(Int, lines[1])
    m = parse(Int, lines[2])

    # Split the input into the down and right matrices
    sections = split(join(lines[3:end], "\n"), "\n-\n")

    down_lines = split(sections[1], '\n')
    right_lines = split(sections[2], '\n')

    # Convert to matrices of Int
    down = convert(Matrix{Int}, hcat([parse.(Int, split(l)) for l in down_lines]...)')
        right = convert(Matrix{Int}, hcat([parse.(Int, split(l)) for l in right_lines]...)')

            # Compute longest path
            max_dist = manhattan_tourist(n, m, down, right)

            println(max_dist)

        end

        # Run the main function
        main()
