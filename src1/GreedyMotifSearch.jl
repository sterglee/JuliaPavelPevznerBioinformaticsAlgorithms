#!/usr/bin/env julia
"""
Greedy Motif Search
Bioinformatics Algorithms (Part 1) – Coursera
"""

using Printf, Statistics

# -----------------------------
# Hamming distance
# -----------------------------
function hamming_distance(s1::AbstractString, s2::AbstractString)
    sum(c1 != c2 for (c1, c2) in zip(s1, s2))
    end

    # -----------------------------
    # Score of motifs
    # -----------------------------
    function score(motifs::Vector{String})
        k = length(motifs[1])
        t = length(motifs)
        total_score = 0
        for i in 1:k
            column = [motifs[j][i] for j in 1:t]
                counts = Dict(nuc => count(==(nuc), column) for nuc in "ACGT")
                    total_score += t - maximum(values(counts))
                end
                return total_score
            end

            # -----------------------------
            # Profile matrix
            # -----------------------------
            function profile(motifs::Vector{String})
                k = length(motifs[1])
                t = length(motifs)
                prof = Vector{Vector{Float64}}(undef, k)
                for i in 1:k
                    column = [motifs[j][i] for j in 1:t]
                        prof[i] = [count(==(nuc), column)/t for nuc in "ACGT"]
                        end
                        return prof
                    end

                    # -----------------------------
                    # Profile-most probable k-mer
                    # -----------------------------
                    function profile_most_probable_kmer(dna::String, k::Int, prof::Vector{Vector{Float64}})
                        nuc_index = Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4)
                        max_prob = -1.0
                        best_kmer = ""
                        for i in 1:(length(dna)-k+1)
                            kmer = dna[i:i+k-1]
                            prob = 1.0
                            for j in 1:k
                                prob *= prof[j][nuc_index[kmer[j]]]
                            end
                            if prob > max_prob
                                max_prob = prob
                                best_kmer = kmer
                            end
                        end
                        return best_kmer
                    end

                    # -----------------------------
                    # Greedy Motif Search
                    # -----------------------------
                    function greedy_motif_search(dna_list::Vector{String}, k::Int, t::Int)
                        best_motifs = [dna[1:k] for dna in dna_list]
                            n = length(dna_list[1])

                            for i in 1:(n - k + 1)
                                motifs = [dna_list[1][i:i+k-1]]
                                for j in 2:t
                                    prof = profile(motifs)
                                    next_motif = profile_most_probable_kmer(dna_list[j], k, prof)
                                    push!(motifs, next_motif)
                                end
                                if score(motifs) < score(best_motifs)
                                    best_motifs = copy(motifs)
                                end
                            end

                            return best_motifs
                        end

                        # -----------------------------
                        # Main Program
                        # -----------------------------
                        lines = readlines("data/stepic_3d.txt")
                        k, t = parse.(Int, split(String(lines[1])))
                        dna_list = [String(strip(line)) for line in lines[2:end]]

                            best_motifs = greedy_motif_search(dna_list, k, t)

                            println(join(best_motifs, "\n"))

