#!/usr/bin/env julia
"""
De Bruijn Graph from a String Problem
Bioinformatics Algorithms (Part 1) – Coursera
"""

# -----------------------------
# Read input
# -----------------------------


dna= "GACCCTTCGGCAGGATTCTAATAACTACTGACATATCAGATTCGGTTGCCTATCTAGCGTGAGCTCATATCCATGCTCTACATCCAGCTTTTAAGCCGAAGGCAGGTCGCTTCTCCGCTTCCCTGAATTCACGGGCCGACGATACGGTGACTAAAGTTTGGCCCGCCACCAAGTTCCCAAGGCTATCACCGAAAGAAAGCGGACGGTT"
k=3

     # -----------------------------
# Build the De Bruijn graph
# -----------------------------
de_bruijn_dict = Dict{String, Set{String}}()

for i in 1:(length(dna)-k+1)
    kmer = dna[i:i+k-1]
    prefix = kmer[1:end-1]
    suffix = kmer[2:end]
    if haskey(de_bruijn_dict, prefix)
        push!(de_bruijn_dict[prefix], suffix)
    else
        de_bruijn_dict[prefix] = Set([suffix])
    end
end

# -----------------------------
# Convert graph to string format
# -----------------------------
de_bruijn_lines = String[]
for (prefix, suffixes) in de_bruijn_dict
    push!(de_bruijn_lines, string(prefix, " -> ", join(collect(suffixes), ",")))
end

# -----------------------------
# Output results
# -----------------------------
println(join(de_bruijn_lines, "\n"))

