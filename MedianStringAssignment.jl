#!/usr/bin/env julia
"""
Median String Problem (Stepic 3B) – Pure Base Julia version
"""

using Printf

# -----------------------------
# Hamming distance
# -----------------------------
function hamming_distance(s1::AbstractString, s2::AbstractString)
    @assert length(s1) == length(s2)
    sum(c1 != c2 for (c1, c2) in zip(s1, s2))
    end

    # -----------------------------
    # Motif score
    # -----------------------------
    function motif_score(pattern::String, motif::String)
        k = length(pattern)
        scores = [hamming_distance(pattern, motif[i:i+k-1]) for i in 1:(length(motif)-k+1)]
            return minimum(scores)
        end

        # -----------------------------
        # Generate all k-mers recursively
        # -----------------------------
        function generate_kmers(k::Int, alphabet::Vector{Char})
            if k == 1
                return [string(c) for c in alphabet]
                else
                    shorter = generate_kmers(k-1, alphabet)
                    kmers = String[]
                    for s in shorter
                        for c in alphabet
                            push!(kmers, s * string(c))
                        end
                    end
                    return kmers
                end
            end

            # -----------------------------
            # Main program
            # -----------------------------


        # -----------------------------
        # Read input
        # -----------------------------
            k = 10
            dna_list = ["GTGCGTATTAAAATAGGAACGGATCATCTGATTAAACAGGGG",
                        "GTCCGAAGATCGCTAATAGCCAAGCAGGAGAGTCTACTGGAC",
                        "GAGCTACAGGGGTCAGCTCTGCGGGTAGTGCTCAGATAACGA",
                        "GGGATCCAGGAGCAACTGTTCTCACTCCTCAGTGTCCACTCT",
                        "ACGCATCCCGGGCAGGTGCACCATAGTGCGAGAAGCAATCAT",
                        "CGGATCATTTACAAGTCCCAAGTCAGTTTCCAGGCGTCGTAG",
                        "CGATTGGATTGGTACCCTGGCCCTGGTAACGAGGGTCAGGAG",
                        "CAGGCGGCTGGATGTGTATACTATGGAAAACCAGGAAGCAGT",
                        "GATGTTCCTAGACAGGCGTTTTACAAGTTCCGTTAATCTGTG",
                        "ACCGTTTGTTGAAGATTGAGATTGACTAGACAGGAGTCAACC"]

            # Initialize best pattern
            best_pattern = (k*length(dna_list)+1, "")

            nucleotides = ['A','C','G','T']

            # Generate all k-mers
            all_patterns = generate_kmers(k, nucleotides)

            # Loop over all k-mers
            for pattern in all_patterns
                current_score = sum(motif_score(pattern, dna) for dna in dna_list)
                    if current_score < best_pattern[1]  # Update if better
                        best_pattern = (current_score, pattern)
                    end
                end

                # Output
                println(best_pattern[2])



