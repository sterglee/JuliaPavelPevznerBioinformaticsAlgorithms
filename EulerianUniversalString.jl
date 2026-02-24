#!/usr/bin/env julia
"""
Binary De Bruijn (Universal String) Generator using Eulerian Cycles
for a given order k.

    Example:
    universal_string_problem(3) -> "00010111"
    """

    # -----------------------------
    # Generate all binary k-mers
    # -----------------------------
    function all_binary_kmers(k::Int)
        if k == 0
            return [""]
        else
            smaller = all_binary_kmers(k-1)
            return [ "0"*s for s in smaller ] ∪ [ "1"*s for s in smaller ]
            end
            end

            # -----------------------------
            # Eulerian cycle in a directed graph
            # -----------------------------
            function eulerian_cycle(graph::Dict{String, Vector{String}})
                g = deepcopy(graph)  # work on a copy to avoid modifying the original
                path = String[]
                stack = [first(keys(g))]  # start at any node

                while !isempty(stack)
                    current = last(stack)
                    if isempty(g[current])
                        push!(path, pop!(stack))
                    else
                        push!(stack, pop!(g[current]))
                    end
                end

                return reverse(path)
            end

            # -----------------------------
            # Generate universal binary string
            # -----------------------------
            function universal_string_problem(k::Int)
                # Build De Bruijn graph: nodes are (k-1)-mers, edges are k-mers
                universal_dict = Dict{String, Vector{String}}()
                for kmer in all_binary_kmers(k)
                    prefix = kmer[1:end-1]
                    suffix = kmer[2:end]
                    if haskey(universal_dict, prefix)
                        push!(universal_dict[prefix], suffix)
                    else
                        universal_dict[prefix] = [suffix]
                    end
                end

                # Find an Eulerian cycle through the graph
                path = eulerian_cycle(universal_dict)

                # Reconstruct the universal string from the Eulerian cycle
                universal_string = join([item[1] for item in path[1:end-1]]) * last(path)[end]

                    return universal_string
                end

                # -----------------------------
                # Example usage
                # -----------------------------
                k = 3
                println("Binary De Bruijn string of order $k:")
                println(universal_string_problem(k))  # Example output: "00010111"

