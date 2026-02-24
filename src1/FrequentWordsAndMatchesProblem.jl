#!/usr/bin/env julia

#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Textbook: Bioinformatics Algorithms: An Active-Learning Approach
Authors: Phillip Compeau & Pavel Pevzner

Problem Title: Frequent Words with Mismatches Problem
=#

using Combinatorics  # For combinations()

# -----------------------------
# Generate Mismatches
# -----------------------------

const NUCLEOTIDES = ['A','C','G','T']

"""
    create_mismatches(swap_list::Vector{Tuple{String,Vector{Int}}})

Recursive function to generate k-mers with mismatches at specified indices.
"""
function create_mismatches(swap_list)
    mismatch_list = []

    # Multiple indices remaining
    if length(swap_list[1][2]) > 1
        for (kmer, indicies) in swap_list
            index = indicies[1]
            for nuc in NUCLEOTIDES
                if nuc != kmer[index]
                    push!(mismatch_list, (replace(kmer, index => nuc), indicies[2:end]))
                end
            end
        end
        return create_mismatches(mismatch_list)
    else
        # Last mismatch
        for (kmer, indicies) in swap_list
            index = indicies[1]
            for nuc in NUCLEOTIDES
                if nuc != kmer[index]
                    push!(mismatch_list, replace(kmer, index => nuc))
                end
            end
        end
        return mismatch_list
    end
end

"""
    mismatch_list(kmer::String, d::Int)

Returns all k-mers that differ by at most d mismatches.
"""
function mismatch_list(kmer::String, d::Int)
    all_mismatches = [kmer]
    for i in 1:d
        combos = collect(combinations(1:length(kmer), i))
        swap_list = [(kmer, combo) for combo in combos]
        append!(all_mismatches, create_mismatches(swap_list))
    end
    return all_mismatches
end

# -----------------------------
# Main
# -----------------------------
dna = "GTACAGACGGACAAAAATTGTACAGACGAGTGTTGCCCAAGTTAGAACATTCAACTGTACAGACGATTCAACTAAGTTAGAACAGTGTTGCCCATTCAACTAAGTTAGAACGACAAAAATTATTCAACTAAGTTAGAACATTCAACTGTACAGACGAGTGTTGCCCAGTGTTGCCCAGTGTTGCCCAGTGTTGCCCATTCAACTAGTGTTGCCCATTCAACTAGTGTTGCCCGACAAAAATTAGTGTTGCCCAGTGTTGCCCAGTGTTGCCCAAGTTAGAACAAGTTAGAACGACAAAAATTATTCAACTAGTGTTGCCCAGTGTTGCCCGTACAGACGGTACAGACGAAGTTAGAACAAGTTAGAACAAGTTAGAACGACAAAAATTAAGTTAGAACGACAAAAATTAAGTTAGAACGTACAGACGAAGTTAGAACAGTGTTGCCCAGTGTTGCCCGTACAGACGGTACAGACGAGTGTTGCCCGTACAGACGATTCAACTGACAAAAATTATTCAACTAGTGTTGCCCGTACAGACGGTACAGACGATTCAACTATTCAACTAGTGTTGCCCGTACAGACGGTACAGACGATTCAACTATTCAACTAAGTTAGAACAGTGTTGCCCGTACAGACGGTACAGACGATTCAACTATTCAACTAGTGTTGCCCATTCAACTAGTGTTGCCCGTACAGACGATTCAACTAGTGTTGCCCATTCAACTGTACAGACGGTACAGACGATTCAACTAGTGTTGCCCGTACAGACGATTCAACTAGTGTTGCCCGACAAAAATTAAGTTAGAACGACAAAAATTATTCAACTATTCAACTGACAAAAATTATTCAACTAAGTTAGAACATTCAACTGACAAAAATTAAGTTAGAACGACAAAAATTATTCAACTAAGTTAGAACAGTGTTGCCCATTCAACTATTCAACTATTCAACTGTACAGACG"

d = 0
k=5
# Count k-mers with up to d mismatches
mismatch_dict = Dict{String, Int}()

n = length(dna)
@inbounds for i in 1:(n - k + 1)
    kmer = dna[i:i + k - 1]
    kmer
    for m in mismatch_list(kmer, d)
        mismatch_dict[m] = get(mismatch_dict, m, 0) + 1
    end
end

# Find max count
max_val = maximum(values(mismatch_dict))

# Collect most frequent k-mers
kmers = [kmer for (kmer, count) in mismatch_dict if count == max_val]

# -----------------------------
# Print & Save Output
# -----------------------------
result = join(kmers, " ")
println(result)


