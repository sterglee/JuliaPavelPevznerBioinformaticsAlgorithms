#!/usr/bin/env julia
"""
Contig Generation Problem (Julia version)
"""

using DataStructures

# -----------------------------
# Read input data
# -----------------------------
filename = "data/ContigGeneration.txt"
kmers = open(filename) do io
    [chomp(line) for line in readlines(io)]
    end

    # -----------------------------
    # Build adjacency dictionary
    # -----------------------------
    edges = Dict{String, Vector{String}}()
    for kmer in kmers
        prefix = kmer[1:end-1]
        suffix = kmer[2:end]
        if haskey(edges, prefix)
            push!(edges[prefix], suffix)
        else
            edges[prefix] = [suffix]
        end
    end

    # -----------------------------
    # Determine balanced and unbalanced nodes
    # -----------------------------
    out_values = reduce(vcat, values(edges))
    nodes = union(keys(edges), out_values)
    balanced, unbalanced = Set{String}(), Set{String}()

    for node in nodes
        outdeg = count(x -> x == node, out_values)
        indeg = haskey(edges, node) ? length(edges[node]) : 0
        if indeg == 1 && outdeg == 1
            push!(balanced, node)
        else
            push!(unbalanced, node)
        end
    end

    # -----------------------------
    # Recursive contig generation
    # -----------------------------
    function get_contigs(start::String, path::String, edges::Dict{String, Vector{String}}, balanced::Set{String})
        contigs = String[]
        if haskey(edges, start)
            for next_node in edges[start]
                new_path = path * next_node[end]
                if next_node in balanced
                    append!(contigs, get_contigs(next_node, new_path, edges, balanced))
                else
                    push!(contigs, new_path)
                end
            end
        else
            push!(contigs, path)
        end
        return contigs
    end

    # Generate all contigs starting from unbalanced nodes
    starts = intersect(unbalanced, keys(edges))
    all_contigs = String[]
    for s in starts
        append!(all_contigs, get_contigs(s, s, edges, balanced))
    end

    # Sort and print/save
    all_contigs = sort(all_contigs)
    println(join(all_contigs, "\n"))


