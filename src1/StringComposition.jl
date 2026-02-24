#!/usr/bin/env julia
"""
String Composition Problem
Bioinformatics Algorithms (Part 1) – Coursera
"""

# -----------------------------
# Read input
# -----------------------------
lines = readlines("data/stepic_4a.txt")
k = parse(Int, strip(lines[1]))
text = strip(lines[2])

# -----------------------------
# Generate all k-mers and sort lexicographically
# -----------------------------
composition = sort([text[i:i+k-1] for i in 1:(length(text)-k+1)])

    # -----------------------------
    # Output results
    # -----------------------------
    println(join(composition, "\n"))

