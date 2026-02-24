#!/usr/bin/env julia
"""
Randomized Motif Search
Bioinformatics Algorithms (Part 1) – Coursera
"""

using Random, Statistics, Printf

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
            # Profile with pseudocounts
            # -----------------------------
            function profile_with_pseudocounts(motifs::Vector{String})
                k = length(motifs[1])
                t = length(motifs)
                prof = Vector{Vector{Float64}}(undef, k)
                for i in 1:k
                    column = [motifs[j][i] for j in 1:t]
                        prof[i] = [(count(==(nuc), column) + 1)/(t + 4) for nuc in "ACGT"]
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
                    # Motifs from profile
                    # -----------------------------
                    function motifs_from_profile(profile::Vector{Vector{Float64}}, dna_list::Vector{String}, k::Int)
                        [profile_most_probable_kmer(seq, k, profile) for seq in dna_list]
                        end

                        # -----------------------------
                        # Randomized Motif Search
                        # -----------------------------
                        function randomized_motif_search(dna_list::Vector{String}, k::Int, t::Int)
                            # Randomly select initial k-mers
                            rand_starts = [rand(1:(length(dna_list[1])-k+1)) for _ in 1:t]
                                motifs = [dna_list[i][rand_starts[i]:rand_starts[i]+k-1] for i in 1:t]

                                    best_motifs = copy(motifs)
                                    best_score = score(best_motifs)

                                    while true
                                        prof = profile_with_pseudocounts(motifs)
                                        motifs = motifs_from_profile(prof, dna_list, k)
                                        current_score = score(motifs)
                                        if current_score < best_score
                                            best_motifs = copy(motifs)
                                            best_score = current_score
                                        else
                                            return best_motifs, best_score
                                        end
                                    end
                                end

                                # -----------------------------
                                # Main Program
                                # -----------------------------
                                lines = readlines("data/stepic_3f.txt")
                                k, t = parse.(Int, split(lines[1]))
                                dna_list = [String(strip(line)) for line in lines[2:end]]

                                    best_motifs, best_score = [], k*t  # initialize

                                    # Repeat the randomized search 1000 times
                                    for _ in 1:1000
                                        motifs, s = randomized_motif_search(dna_list, k, t)
                                        if s < best_score
                                            best_motifs = copy(motifs)
                                            best_score = s
                                        end
                                    end

                                    println(join(best_motifs, "\n"))

