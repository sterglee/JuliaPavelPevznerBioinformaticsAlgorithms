#!/usr/bin/env julia

#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Textbook: Bioinformatics Algorithms: An Active-Learning Approach
Authors: Phillip Compeau & Pavel Pevzner

Problem Title: Minimum Skew Problem
Assignment #: 01
Problem ID: E
=#

# -----------------------------
# Read Input
# -----------------------------
input_path = "data/stepic_1e.txt"
dna = open(input_path, "r") do io
    strip(read(io, String))
end

# -----------------------------
# Compute Skew & Minimum Skew Indices
# -----------------------------
skew_value = 0
min_skew = 1   # Start with Python's initial min_skew = 1
min_ind = Int[]

for (index, nucleotide) in enumerate(dna)
    if nucleotide == 'C'
        skew_value -= 1
    elseif nucleotide == 'G'
        skew_value += 1
    end

    if skew_value == min_skew
        push!(min_ind, index)       # Julia is 1-based, matches Python's index+1
    elseif skew_value < min_skew
        min_skew = skew_value
        min_ind = [index]
    end
end

# -----------------------------
# Print & Save Output
# -----------------------------
result = join(min_ind, " ")
println(result)

output_path = "output/Assignment_01E.txt"
open(output_path, "w") do io
    write(io, result)
end

