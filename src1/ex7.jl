
function CheckClumpLength(indicies, t, L)
      for i in 1:(length(indicies) - t + 1)
        if indicies[t + i - 1] - indicies[i] <= L
            return true
        end
    end
    return false
end

k=0
dna=""
open("data/stepic_1d.txt") do input_data
    lines = readlines(input_data)
    global dna = strip(lines[1])
    global k, L, t = parse.(Int, split(strip(lines[2])))
end


# Find all k-mers, count their appearances, and store their indices. 
kmer_dict = Dict{String, Tuple{Int, Vector{Int}}}()
for i in 1:(length(dna) - k + 1)
    kmer = dna[i:i+k-1]
    if haskey(kmer_dict, kmer)
        kmer_dict[kmer] = (kmer_dict[kmer][1] + 1, push!(kmer_dict[kmer][2], i))
    else
        kmer_dict[kmer] = (1, [i])
    end
end

# The candidate k-mers that appear at least t times, along with the indices where they appear.
kmer_candidates = [(kmer[1], kmer[2][2]) for kmer in kmer_dict if kmer[2][1] >= t]

# Check that at least t candidate k-mers fall within a clump of size L.
kmer_clumps = String[]
for candidate in kmer_candidates
    if CheckClumpLength(candidate[2], t, L)
        push!(kmer_clumps, candidate[1])
    end
end

# Print and save the solution.
println(join(kmer_clumps))
open("output/Assignment_01D.txt", "w") do output_data
    write(output_data, join(kmer_clumps))
end
