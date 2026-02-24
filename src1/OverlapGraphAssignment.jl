#!/usr/bin/env julia
"""
Overlap Graph Problem
Bioinformatics Algorithms (Part 1) – Coursera
"""

# -----------------------------
# Read input
# -----------------------------
dna = strip.(readlines("data/stepic_4b.txt"))

# -----------------------------
# Define overlap check and printing
# -----------------------------
check_overlap(pair) = pair[1][1:end] == pair[2][1:end]  # placeholder
check_overlap(dna1, dna2, k) = dna1[2:end] == dna2[1:end-1]  # overlap of k-1
print_overlap(dna1, dna2) = string(dna1, " -> ", dna2)

# -----------------------------
# Generate all overlapping pairs
# -----------------------------
overlaps = String[]
for i in 1:length(dna)
    for j in 1:length(dna)
        if i != j && dna[i][2:end] == dna[j][1:end-1]  # k-1 overlap
            push!(overlaps, print_overlap(dna[i], dna[j]))
        end
    end
end

# -----------------------------
# Output results
# -----------------------------
println(join(overlaps, "\n"))

