#!/usr/bin/env julia
"""
Eulerian Path Problem
Bioinformatics Algorithms (Part 1) – Coursera
"""

using Base.Iterators: flatten

"""
Generates an Eulerian cycle from the given edges.
Input: Dict{Int, Vector{Int}} mapping each node to its outgoing neighbors.
Output: Vector{Int} representing the Eulerian cycle.
"""
function eulerian_cycle(edge_dict::Dict{Int, Vector{Int}})
    # Copy the dictionary to avoid mutating the original
    edges = deepcopy(edge_dict)

    # Start at any node
    current_node = first(keys(edges))
    path = [current_node]

    # Build the initial cycle
    while true
        next_node = edges[current_node][1]
        push!(path, next_node)

        if length(edges[current_node]) == 1
            delete!(edges, current_node)
        else
            edges[current_node] = edges[current_node][2:end]
        end

        if haskey(edges, path[end])
            current_node = path[end]
        else
            break
        end
    end

    # Expand the cycle until all edges are used
    while !isempty(edges)
        for i in 1:length(path)
            if haskey(edges, path[i])
                current_node = path[i]
                cycle = [current_node]

                while true
                    next_node = edges[current_node][1]
                    push!(cycle, next_node)

                    if length(edges[current_node]) == 1
                        delete!(edges, current_node)
                    else
                        edges[current_node] = edges[current_node][2:end]
                    end

                    if haskey(edges, cycle[end])
                        current_node = cycle[end]
                    else
                        break
                    end
                end

                # Merge the new cycle into the existing path
                path = vcat(path[1:i-1], cycle, path[i+1:end])
                break
            end
        end
    end

    return path
end



# -----------------------------
# Eulerian path function
# -----------------------------
function eulerian_path(edge_dict::Dict{Int, Vector{Int}})
    # Determine unbalanced edges
    out_values = collect(flatten(values(edge_dict)))
    nodes = union(keys(edge_dict), out_values)
    unbalanced_from = nothing
    unbalanced_to = nothing
    for node in nodes
        out_val = count(x -> x == node, out_values)
        in_val = haskey(edge_dict, node) ? length(edge_dict[node]) : 0
        if in_val < out_val
            unbalanced_from = node
            elseif out_val < in_val
            unbalanced_to = node
        end
    end

    # Add an edge connecting the unbalanced edges
    if haskey(edge_dict, unbalanced_from)
        push!(edge_dict[unbalanced_from], unbalanced_to)
    else
        edge_dict[unbalanced_from] = [unbalanced_to]
    end

    # Get the Eulerian cycle including the unbalanced edge
    cycle = eulerian_cycle(edge_dict)

    # Find the location of the added edge
    divide_point = findfirst(i -> cycle[i:i+1] == [unbalanced_from, unbalanced_to], 1:length(cycle)-1)

    # Remove the added edge and rotate the cycle
    return vcat(cycle[divide_point+1:end], cycle[2:divide_point+1])
end

# -----------------------------
# Main Program
# -----------------------------
edges = Dict{Int, Vector{Int}}()
open("data/stepic_5a.txt") do f
    for line in readlines(f)
        parts = split(strip(line), " -> ")
        key = parse(Int, parts[1])
        values = split(parts[2], ",")
        edges[key] = [parse(Int, v) for v in values]
        end
    end

    # Compute Eulerian path
    path = eulerian_path(edges)

    # Output results
    println(join(path, "->"))
