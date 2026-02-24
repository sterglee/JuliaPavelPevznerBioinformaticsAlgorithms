#!/usr/bin/env julia
"""
A solution to the (k,d)-Motif Enumeration assignment (Stepic 3A)
Converted from Python to Julia.
"""

using Combinatorics  # For generating mismatches
using Printf

# -----------------------------
# Function to generate all (k,d)-mismatches
# -----------------------------
function mismatch_list(kmer::String, d::Int)
    """
    Returns all k-mers within Hamming distance d of the input kmer.
    """
    nucleotides = ['A','C','G','T']
    k = length(kmer)
    if d == 0
        return Set([kmer])
    elseif k == 0
        return Set([""])
    end

    # Recursive generation of mismatches
    first_char = kmer[1]
    suffix = kmer[2:end]

    suffix_mismatches = mismatch_list(suffix, d)
    result = Set{String}()
    for sm in suffix_mismatches
        # If no mismatch used, we can either keep the first_char or change it
        if count(x -> x != first_char, sm) < d
            for n in nucleotides
                push!(result, string(n) * sm)
            end
        else
            push!(result, string(first_char) * sm)
        end
    end

    return result
end

# -----------------------------
# Read input
# -----------------------------
       dna_list = ["TCCGTCTTGCGGTAGCGCACCTCTG", "GGAGCCATATGGCAGGTTATGACAA", "CTGCATATGGTTCTACCGTCGGTAG",  "TATGAGGAAGACAATGTCTTCACGG",
"CGCGTACCCTTCCTAGGGTTGGTAG", "AACCAGGCAGAAGTCTGTCTTACTT"]
  k=3
  d=1
    # Generate sets of (k,d)-motifs for each DNA sequence
    motif_sets = []
    for dna in dna_list
        motifs = Set{String}()
        for i in 1:(length(dna)-k+1)
            kmer = dna[i:i+k-1]
            union!(motifs, mismatch_list(kmer, d))
        end
        push!(motif_sets, motifs)
    end

    # Intersect all sets to get common motifs
    common_motifs = reduce(intersect, motif_sets)

    # Sort motifs
    sorted_motifs = sort(collect(common_motifs))

    # Output
    println(join(sorted_motifs, " "))


