#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
The associated textbook is Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
Problem Title: Number of Breakpoints Problem
=#

"""
breakpoint_count(permutation::Vector{Int})

Returns the number of breakpoints in a given permutation by prepending 0
and appending n+1, then counting adjacent elements that are not consecutive.
"""
function breakpoint_count(permutation)
    # Prepend 0 and append len+1 to check endpoints
    extended_perm = [0; permutation; length(permutation) + 1]

    # Compare adjacent elements (x, y) and count where x - y != 1
    # In Julia, we can use slicing and the dot operator for element-wise comparison
    return sum(extended_perm[2:end] .- extended_perm[1:end-1] .!= 1)
end

function main()
    # Read the input data
    # We use 'replace' to clean up the parentheses common in these dataset formats
    raw_data = read("data/BreakpointCount.txt", String)
    clean_data = replace(strip(raw_data), r"\(|\)" => "")

    # Convert string to Vector{Int}
    perm = parse.(Int, split(clean_data))

    # Get the number of breakpoints
    num_of_breakpoints = breakpoint_count(perm)

    # Print the answer
    println(num_of_breakpoints)

    # Save the answer
    mkpath("output") # Ensure directory exists
    write("output/BreakpointCount.txt", string(num_of_breakpoints))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
