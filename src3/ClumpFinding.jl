function clump_finding(genome_string, k, L, t)
    # Initialize list of clumps
    clump_list = String[]
    n = length(genome_string)
    
    # Loop over each L-length interval
    # range(0, n-L+1) in Python becomes 1:(n-L+1) in Julia
    for u in 1:(n - L + 1)
        tmp_string = genome_string[u : u + L - 1]
    
        # Initialize list of k-mers checked in this window
        pattern_list = String[]
    
        # Loop over each k-mer within the window
        for i in 1:(L - k + 1)
            # Get current k-mer
            pattern = tmp_string[i : i + k - 1]

            # Check if this is a NEW pattern in this window
            if pattern ∉ pattern_list
                push!(pattern_list, pattern)

                # Count frequency within the window
                count = 0
                for j in i:(L - k + 1)
                    if pattern == tmp_string[j : j + k - 1]
                        count += 1
                    end
                end

                # If frequent enough, add to global clump list
                if count >= t
                    if pattern ∉ clump_list
                        push!(clump_list, pattern)
                    end
                end
            end
        end
    end

    return clump_list
end

# Test Case
dna_string = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC"
println("Sample test: ", clump_finding(dna_string, 5, 75, 4))

