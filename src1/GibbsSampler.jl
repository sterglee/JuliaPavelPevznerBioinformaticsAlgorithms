#!/usr/bin/env julia
#Gibbs Sampler
#Bioinformatics Algorithms (Part 1) – Coursera

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
# Gibbs Sampler
# -----------------------------
function gibbs_sampler(dna_list::Vector{String}, k::Int, t::Int, N::Int)
    # Randomly generate k-mers from each sequence
    rand_starts = [rand(1:(length(dna_list[1])-k+1)) for _ in 1:t]
    motifs = [dna_list[i][rand_starts[i]:rand_starts[i]+k-1] for i in 1:t]

    best_motifs = copy(motifs)
    best_score = score(best_motifs)

    for _ in 1:N
        r = rand(1:t)
        # Remove r-th motif
        motifs_except_r = [motifs[i] for i in 1:t if i != r]
        current_profile = profile_with_pseudocounts(motifs_except_r)
        # Sample new motif for position r
        motifs[r] = profile_most_probable_kmer(dna_list[r], k, current_profile)
        current_score = score(motifs)
        if current_score < best_score
            best_motifs = copy(motifs)
            best_score = current_score
        end
    end
    return best_motifs, best_score
end

# -----------------------------
# Main Program
# -----------------------------
lines = readlines("data/stepic_3g.txt")
k, t, N = parse.(Int, split(lines[1]))
dna_list = [String(strip(line)) for line in lines[2:end]]

best_motifs, best_score = [], k*t  # Initialize best score

# Repeat Gibbs sampler 20 times
for _ in 1:20
    motifs, s = gibbs_sampler(dna_list, k, t, N)
    if s < best_score
        best_motifs = copy(motifs)
        best_score = s
    end
end

# Output results
println(join(best_motifs, "\n"))
