using Combinatorics
using DataStructures

# --- Change: Use AbstractString to avoid MethodError ---
function get_mismatches(kmer::AbstractString, d::Int)
    nucleotides = ['A', 'C', 'G', 'T']
    k = length(kmer)
    mismatch_set = Set{String}()
    push!(mismatch_set, String(kmer))
    
    if d == 0
        return mismatch_set
    end

    # Iterate through 1 to d possible mismatch counts
    for i in 1:d
        # combinations(1:k, i) gives all sets of positions to mutate
        for indices in combinations(1:k, i)
            # Iterators.product generates all possible nucleotide combinations for those positions
            for replacements in Iterators.product(fill(nucleotides, i)...)
                new_kmer = collect(kmer)
                
                # Apply replacements
                for (idx, char) in zip(indices, replacements)
                    new_kmer[idx] = char
                end
                
                push!(mismatch_set, String(new_kmer))
            end
        end
    end
    return mismatch_set
end

function solve()
    input_path = "data/stepic_1g.txt"
    if !isfile(input_path)
        println("Error: File not found at $input_path")
        return
    end

    # 1. Read input
    lines = readlines(input_path)
    if length(lines) < 2
        println("Error: Input file must have DNA on line 1 and k d on line 2.")
        return
    end

    dna = strip(lines[1])
    params = parse.(Int, split(strip(lines[2])))
    k, d = params[1], params[2]

    # 2. Count mismatches
    # DefaultDict simplifies the counting logic
    mismatch_dict = DefaultDict{String, Int}(0)

    for i in 1:(length(dna) - k + 1)
        # Slicing creates a SubString
        current_window = dna[i:i+k-1]
        
        # Get all neighbors and increment their counts
        for neighbor in get_mismatches(current_window, d)
            mismatch_dict[neighbor] += 1
        end
    end

    # 3. Find and print maximums
    if isempty(mismatch_dict) return end
    
    max_count = maximum(values(mismatch_dict))
    results = [kmer for (kmer, count) in mismatch_dict if count == max_count]

    # Stepik usually expects space-separated results
    println(join(results, " "))
    
    # Optional: Save to output
    # mkpath("output")
    # write("output/Assignment_01G.txt", join(results, " "))
end

solve()

