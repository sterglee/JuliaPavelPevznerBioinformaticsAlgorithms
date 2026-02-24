#!/usr/bin/env julia

using DelimitedFiles

# -----------------------------
# Read input
# -----------------------------
lines = readlines("data/stepic_3c.txt")
dna = strip(lines[1])
k = parse(Int, strip(lines[2]))

# First line of profile is nucleotides
profile_nucs = split(strip(lines[3]))
profile_values = [parse.(Float64, split(strip(line))) for line in lines[4:end]]

# Map nucleotide (Char) => index
nuc_loc = Dict{Char, Int}()
for (i, nuc_str) in enumerate(profile_nucs)
    nuc_loc[first(nuc_str)] = i
end

# -----------------------------
# Initialize maximum probability
# -----------------------------
max_prob = -1.0
max_kmer = ""

# -----------------------------
# Loop over all k-mers in the DNA
# -----------------------------
for i in 1:(length(dna)-k+1)
    kmer = dna[i:i+k-1]  # SubString
    prob = 1.0
    for j in 1:k
        prob *= profile_values[j][nuc_loc[kmer[j]]]
    end
    if prob > max_prob
        max_prob = prob
        max_kmer = String(kmer)
    end
end

# -----------------------------
# Print and save the answer
# -----------------------------
println(max_kmer)
