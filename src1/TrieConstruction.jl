#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Problem Title: Trie Construction Problem
Assignment #: 09
=#

"""
construct_trie(patterns::Vector{String})

Constructs a Trie from a list of patterns. Returns a Dict representing edges:
(parent_node, child_node) => label
"""
function construct_trie(patterns)
    # edges dictionary: (from_node, to_node) => label_char
    edges = Dict{Tuple{Int, Int}, Char}()

    # root node is 1
    new_node_id = 1

    # Nested Dict to represent the tree structure: node_id => {char => next_node_id}
    trie_structure = Dict{Int, Dict{Char, Int}}()
    trie_structure[1] = Dict{Char, Int}()

    for text in patterns
        current_node = 1
        for char in text
            if haskey(trie_structure[current_node], char)
                # Follow existing path
                current_node = trie_structure[current_node][char]
            else
                # Create new node
                new_node_id += 1
                trie_structure[current_node][char] = new_node_id
                trie_structure[new_node_id] = Dict{Char, Int}()

                # Record the edge for output
                edges[(current_node, new_node_id)] = char

                current_node = new_node_id
            end
        end
    end

    return edges
end

function main()
    # Read input data
    input_path = "data/stepic_9a.txt"
    if !isfile(input_path)
        println("Error: Input file not found.")
        return
    end

    patterns = filter(!isempty, strip.(readlines(input_path)))

    # Construct the Trie
    trie_edges_map = construct_trie(patterns)

    # Convert to adjacency format: "parent child symbol"
    # Note: Stepic usually expects 0-indexed nodes, but the problem statement
    # often allows any unique labeling. Here we use 1-based indexing as generated.
    # If 0-indexed is required, use (u-1, v-1).
    adjacency_list = ["$(u-1) $(v-1) $label" for ((u, v), label) in trie_edges_map]

        # Sort for consistent output
        sort!(adjacency_list)

        # Print results
        output_text = join(adjacency_list, "\n")
        println(output_text)

        # Save to file
        mkpath("output")
        write("output/Assignment_09A.txt", output_text)
    end

    if abspath(PROGRAM_FILE) == @__FILE__
        main()
    end

