#!/usr/bin/env julia

#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Textbook: Bioinformatics Algorithms: An Active-Learning Approach
Authors: Phillip Compeau & Pavel Pevzner

Problem Title: Reverse Complement Problem
=#

# -----------------------------
# Reverse Complement Function
# -----------------------------

const DNA_COMPLEMENT = Dict(
    'A' => 'T',
    'T' => 'A',
    'C' => 'G',
    'G' => 'C'
)

function ReverseComplementDNA(seq::AbstractString)
    n = length(seq)
    result = Vector{Char}(undef, n)

    @inbounds for i in 1:n
        result[n - i + 1] = DNA_COMPLEMENT[seq[i]]
    end

    return String(result)
end


# -----------------------------
# Main
# -----------------------------

input_path = "data/ReverseComplement.txt"
output_path = "output/ReverseComplement.txt"

dna = open(input_path, "r") do io
    strip(read(io, String))
end

reverse_complement = ReverseComplementDNA(dna)

# Write result
open(output_path, "w") do io
    write(io, reverse_complement)
end

println(reverse_complement)
