#!/usr/bin/env julia
"""
Convolution Cyclopeptide Sequencing )
"""

using Printf, Statistics

# -----------------------------
# Append proteins to current peptides
# -----------------------------
function append_protein(add_list::Vector{Vector{Int}}, protein_alphabet::Vector{Int})
    newlist = Vector{Vector{Int}}()
    for peptide in add_list
        for p in protein_alphabet
            push!(newlist, vcat(peptide, p))
        end
    end
    return newlist
end

# -----------------------------
# Circular spectrum
# -----------------------------
function spectrum(peptide::Vector{Int})
    n = length(peptide)
    doubled = vcat(peptide, peptide)
    spec = [0, sum(peptide)]
    for i in 1:(n-1)
        for j in 1:n
            subpep = doubled[j:j+i-1]
            push!(spec, sum(subpep))
        end
    end
    sort!(spec)
    return spec
end

# -----------------------------
# Spectrum score
# -----------------------------
function spectrum_score(peptide::Vector{Int}, exp_spec::Vector{Int})
    pep_spec = spectrum(peptide)
    if maximum(pep_spec) > maximum(exp_spec)
        return -1
    end
    score = sum(min(count(x -> x==p, pep_spec), count(x -> x==p, exp_spec)) for p in unique(pep_spec))
        return score
    end

    # -----------------------------
    # Main function
    # -----------------------------
    function convolution_cyclopeptide(M::Int, N::Int, exp_spec::Vector{Int})
        # Step 1: compute the convolution
        convolution = [i - j for i in exp_spec, j in exp_spec if i - j > 0]

            # Step 2: count occurrences of each mass in the range 57..200
            convo_dict = Dict{Int, Int}()
            for c in convolution
                if 57 <= c <= 200
                    convo_dict[c] = get(convo_dict, c, 0) + 1
                end
            end

            # Step 3: build protein alphabet of top M frequent elements (with ties)
            sorted_counts = sort(collect(convo_dict), by=x->x[2], rev=true)
            protein_alphabet = Int[]
            i = 1
            while length(protein_alphabet) < M && i <= length(sorted_counts)
                count_val = sorted_counts[i][2]
                tied = [k for (k,v) in sorted_counts[i:end] if v == count_val]
                protein_alphabet = vcat(protein_alphabet, tied)
                i += length(tied)
                end
                protein_alphabet = protein_alphabet[1:min(M, length(protein_alphabet))]

                # Step 4: initialize leader
                overall_leader = Vector{Int}()
                overall_score = -1

                # Step 5: initialize peptides
                seq = [[p] for p in protein_alphabet if spectrum_score([p], exp_spec) != -1]

                    # Step 6: iterative expansion
                    while !isempty(seq)
                        # Score all current peptides
                        scored = [(spectrum_score(peptide, exp_spec), peptide) for peptide in seq]
                            scored = filter(x -> x[1] != -1, scored)

                            if isempty(scored)
                                break
                            end

                            # Select top N leaders (including ties)
                            sorted_scores = sort(scored, by=x->x[1], rev=true)
                            leaders = Vector{Vector{Int}}()
                            topN_score = sorted_scores[min(N, length(sorted_scores))][1]
                            for (sc, pep) in sorted_scores
                                if sc >= topN_score
                                    push!(leaders, pep)
                                else
                                    break
                                end
                            end

                            # Update overall leader
                            for (sc, pep) in scored
                                if sc > overall_score
                                    overall_score = sc
                                    overall_leader = pep
                                end
                            end

                            # Expand leaders
                            seq = append_protein(leaders, protein_alphabet)
                            seq = [p for p in seq if spectrum_score(p, exp_spec) != -1]
                            end

                            # Convert leader peptide to mass string
                            leader_peptide_str = join(overall_leader, "-")
                            println(leader_peptide_str)
                            return leader_peptide_str
                            end

                            # -----------------------------
                            # Example usage
                            # -----------------------------
                            # Read data
                            M, N, exp_spec = 20, 60, sort([0, 71, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186])  # replace with file read if needed

                            # Run the algorithm
                            convolution_cyclopeptide(M, N, exp_spec)

