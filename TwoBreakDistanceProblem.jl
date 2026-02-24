#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Problem Title: 2-Break Distance Problem
Assignment #: 08
=#

"""
two_break_dist(P::Vector{Vector{Int}}, Q::Vector{Vector{Int}})

Returns the 2-Break Distance of Circular Chromosomes P and Q by counting
cycles in the breakpoint graph.
"""
function two_break_dist(P, Q)
    # Construct the breakpoint graph of P and Q.
    # We use a Dict mapping an Int to a Vector of Ints (the connected nodes).
    edges = Dict{Int, Vector{Int}}()

    # Helper function to add edges to the adjacency list
    function add_edge!(u, v)
        push!(get!(edges, u, Int[]), v)
        push!(get!(edges, v, Int[]), u)
    end

    # Build the graph from both sets of chromosomes
    for genome in (P, Q)
        for block in genome
            L = length(block)
            for i in 1:L
                # Julia is 1-indexed. The next element with circular wrap-around:
                u = block[i]
                v = -block[mod1(i + 1, L)]
                add_edge!(u, v)
            end
        end
    end

    # Count the number of cycles in the breakpoint graph using a simple traversal.
    cycles = 0
    while !isempty(edges)
        cycles += 1
        # Get an arbitrary starting node
        start_node = first(keys(edges))
        current = start_node

        while haskey(edges, current)
            # Pick the first available neighbor
            neighbor = popfirst!(edges[current])

            # Clean up the current node if no edges left
            if isempty(edges[current])
                delete!(edges, current)
            end

            # Remove the back-edge from the neighbor to the current node
            if haskey(edges, neighbor)
                # Find the index of 'current' in neighbor's list and remove it
                idx = findfirst(x -> x == current, edges[neighbor])
                deleteat!(edges[neighbor], idx)
                if isempty(edges[neighbor])
                    delete!(edges, neighbor)
                end
            end

            current = neighbor
            # If we returned to start, this specific cycle traversal is done
            if current == start_node
                break
            end
        end
    end

    # Theorem: d(P,Q) = total_blocks - cycles
    total_blocks = sum(length(block) for block in P)
        return total_blocks - cycles
    end

    function main()
        # Read and parse input data
        # The input typically looks like: (+1 +2 +3)(-4 -5)
        lines = readlines("data/stepic_8c.txt")

        function parse_genome(line)
            # Remove outer parens and split by ")("
            line = strip(line)
            blocks_str = split(replace(line, r"^\(|\)$" => ""), ")(")
                                                                   return [parse.(Int, split(b)) for b in blocks_str]
                                                                   end

                                                                   P = parse_genome(lines[1])
                                                                   Q = parse_genome(lines[2])

                                                                   # Calculate Distance
                                                                   dist = two_break_dist(P, Q)

                                                                   # Output results
                                                                   println(dist)
                                                                   mkpath("output")
                                                                   write("output/Assignment_08C.txt", string(dist))
        end

        if abspath(PROGRAM_FILE) == @__FILE__
            main()
        end

