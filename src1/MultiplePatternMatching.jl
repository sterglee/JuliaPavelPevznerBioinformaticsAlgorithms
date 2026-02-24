#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Problem Title: Multiple Pattern Matching Problem
Assignment #: 09
=#

"""
build_trie(patterns::Vector{String})

Builds a Trie where each node is a Dict mapping characters to the next node ID.
Returns the trie structure and the IDs of nodes that represent the end of a pattern.
"""
function build_trie(patterns)
    # trie[node_id] = Dict(char => next_node_id)
    trie = Dict{Int, Dict{Char, Int}}()
    trie[1] = Dict{Char, Int}()

    # Track which nodes are "leaf" nodes (end of a pattern)
    is_terminal = Set{Int}()

    next_id = 2
    for pattern in patterns
        current_node = 1
        for char in pattern
            if !haskey(trie[current_node], char)
                trie[current_node][char] = next_id
                trie[next_id] = Dict{Char, Int}()
                next_id += 1
            end
            current_node = trie[current_node][char]
        end
        push!(is_terminal, current_node)
    end
    return trie, is_terminal
end

"""
is_prefix_in_trie(text_suffix, trie, is_terminal)

Checks if any pattern in the trie is a prefix of the given text_suffix.
    """
    function is_prefix_in_trie(text_suffix, trie, is_terminal)
        current_node = 1
        for char in text_suffix
            if haskey(trie[current_node], char)
                current_node = trie[current_node][char]
                # If we reach a terminal node, a pattern has been found
                if current_node in is_terminal
                    return true
                end
            else
                # Character not in trie, stop searching this branch
                break
            end
        end
        return false
    end

    function main()
        # Read input
        input_path = "data/stepic_9b.txt"
        if !isfile(input_path)
            println("Error: File not found.")
            return
        end

        lines = readlines(input_path)
        text = strip(lines[1])
        patterns = filter(!isempty, strip.(lines[2:end]))

        # Pre-processing: Build the Trie
        trie, is_terminal = build_trie(patterns)

        # Finding the shortest pattern length to avoid unnecessary checks at the end
        min_len = minimum(length.(patterns))

        # Matching
        matching_indices = Int[]
        for i in 1:(length(text) - min_len + 1)
            # Use a view to avoid memory allocation for the suffix
            if is_prefix_in_trie(@view(text[i:end]), trie, is_terminal)
                push!(matching_indices, i - 1) # Convert to 0-indexed
            end
        end

        # Format and save output
        result = join(matching_indices, " ")
        println(result)

        mkpath("output")
        write("output/Assignment_09B.txt", result)
    end

    if abspath(PROGRAM_FILE) == @__FILE__
        main()
    end

