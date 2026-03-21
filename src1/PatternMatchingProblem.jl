#!/usr/bin/env julia

#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Textbook: Bioinformatics Algorithms: An Active-Learning Approach
Authors: Phillip Compeau & Pavel Pevzner

Problem Title: Pattern Matching Problem
Assignment #: 01
Problem ID: C
=#

# -----------------------------
# Read Input
# -----------------------------
input_path = "data/stepic_1c.txt"
pattern, text = open(input_path, "r") do io
    strip.(readlines(io))
end

# -----------------------------
# Pattern Matching
# -----------------------------
pattern_len = length(pattern)
text_len = length(text)

pattern_loc = Int[]

@inbounds for i in 1:(text_len - pattern_len + 1)
    if text[i:i + pattern_len - 1] == pattern
        push!(pattern_loc, i - 1)  # 0-based indexing to match Python
    end
end

# Convert to string for output
result = join(pattern_loc, " ")
println(result)

# -----------------------------
# Write Output
# -----------------------------
output_path = "output/Assignment_01C.txt"
open(output_path, "w") do io
    write(io, result)
end

