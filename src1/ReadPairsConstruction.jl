#!/usr/bin/env julia
"""
String Construction from Read-Pairs Problem (Julia version)
"""

using DataStructures

# -----------------------------
# Eulerian path from adjacency dictionary
# -----------------------------
function eulerian_path(edge_dict::Dict{Tuple{String,String},Vector{Tuple{String,String}}})
    # Count in-degree and out-degree for each node
    in_deg = Dict{Tuple{String,String}, Int}()
    out_deg = Dict{Tuple{String,String}, Int}()

    for (u, vs) in edge_dict
        out_deg[u] = length(vs)
        for v in vs
            in_deg[v] = get(in_deg, v, 0) + 1
        end
    end

    nodes = union(keys(edge_dict), keys(in_deg))
    start_node = nothing
    end_node = nothing

    for node in nodes
        outd = get(out_deg, node, 0)
        ind = get(in_deg, node, 0)
        if outd - ind == 1
            start_node = node
            elseif ind - outd == 1
            end_node = node
        end
    end

    # If no specific start, pick any node
    start_node = isnothing(start_node) ? first(keys(edge_dict)) : start_node

    # Copy edges for mutation
    edges = Dict{Tuple{String,String}, Vector{Tuple{String,String}}}()
    for (k,v) in edge_dict
        edges[k] = copy(v)
    end

    # Hierholzer's algorithm
    path = Tuple{String,String}[]
    stack = [start_node]

    while !isempty(stack)
        current = last(stack)
        if haskey(edges, current) && !isempty(edges[current])
            push!(stack, pop!(edges[current]))
        else
            push!(path, pop!(stack))
        end
    end

    return reverse(path)
end

# -----------------------------
# Read input data
# -----------------------------
filename = "data/stepic_5d.txt"
open(filename) do io
    d = parse(Int, readline(io))
    paired_reads = [split(chomp(line), '|') for line in readlines(io)]
        k = length(paired_reads[1][1])

        # Build adjacency dictionary
        edge_dict = Dict{Tuple{String,String}, Vector{Tuple{String,String}}}()
        for pair in paired_reads
            prefix = (pair[1][1:end-1], pair[2][1:end-1])
            suffix = (pair[1][2:end], pair[2][2:end])
            if haskey(edge_dict, prefix)
                push!(edge_dict[prefix], suffix)
            else
                edge_dict[prefix] = [suffix]
            end
        end

        # -----------------------------
        # Find Eulerian path
        # -----------------------------
        path = eulerian_path(edge_dict)

        # -----------------------------
        # Reconstruct string
        # -----------------------------
        first_string  = path[1][1] * join(node[1][end] for node in path[2:end])
            second_string = path[1][2] * join(node[2][end] for node in path[2:end])
            text = first_string[1:k+d] * second_string

            # Print output
            println(text)
            end

