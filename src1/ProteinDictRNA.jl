using DelimitedFiles
using DataStructures
using Combinatorics

# -----------------------------
# 1. Helper Functions
# -----------------------------

const DNA_COMPLEMENT = Dict('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C')

function ReverseComplementDNA(seq::AbstractString)
    n = length(seq)
    result = Vector{Char}(undef, n)
    for i in 1:n
        result[n - i + 1] = DNA_COMPLEMENT[seq[i]]
    end
    return String(result)
end

function get_mismatches(kmer::AbstractString, d::Int)
    nucleotides = ['A', 'C', 'G', 'T']
    k = length(kmer)
    mismatch_set = Set{String}()
    push!(mismatch_set, String(kmer))
    
    for i in 1:d
        for indices in combinations(1:k, i)
            for replacements in Iterators.product(fill(nucleotides, i)...)
                new_kmer = collect(kmer)
                for (idx, char) in zip(indices, replacements)
                    new_kmer[idx] = char
                end
                push!(mismatch_set, String(new_kmer))
            end
        end
    end
    return mismatch_set
end

function ProteinDictRNA()
    # Simplified standard genetic code
    return Dict(
        "UUU"=>"F", "UUC"=>"F", "UUA"=>"L", "UUG"=>"L", "UCU"=>"S", "UCC"=>"S", "UCA"=>"S", "UCG"=>"S",
        "UAU"=>"Y", "UAC"=>"Y", "UAA"=>"Stop", "UAG"=>"Stop", "UGU"=>"C", "UGC"=>"C", "UGA"=>"Stop", "UGG"=>"W",
        "CUU"=>"L", "CUC"=>"L", "CUA"=>"L", "CUG"=>"L", "CCU"=>"P", "CCC"=>"P", "CCA"=>"P", "CCG"=>"P",
        "CAU"=>"H", "CAC"=>"H", "CAA"=>"Q", "CAG"=>"Q", "CGU"=>"R", "CGC"=>"R", "CGA"=>"R", "CGG"=>"R",
        "AUU"=>"I", "AUC"=>"I", "AUA"=>"I", "AUG"=>"M", "ACU"=>"T", "ACC"=>"T", "ACA"=>"T", "ACG"=>"T",
        "AAU"=>"N", "AAC"=>"N", "AAA"=>"K", "AAG"=>"K", "AGU"=>"S", "AGC"=>"S", "AGA"=>"R", "AGG"=>"R",
        "GUU"=>"V", "GUC"=>"V", "GUA"=>"V", "GUG"=>"V", "GCU"=>"A", "GCC"=>"A", "GCA"=>"A", "GCG"=>"A",
        "GAU"=>"D", "GAC"=>"D", "GAA"=>"E", "GAG"=>"E", "GGU"=>"G", "GGC"=>"G", "GGA"=>"G", "GGG"=>"G"
    )
end

# -----------------------------
# 2. Main Problems
# -----------------------------

# --- Problem 1H: Frequent Words with Mismatches and Reverse Complements ---
function solve_1h()
    input_path = "data/stepic_1h.txt"
    if !isfile(input_path) return end

    # readdlm returns a matrix; we take the first column
    data = readdlm(input_path, '\n', String)
    dna = strip(data[1])
    k_d = parse.(Int, split(data[2]))
    k, d = k_d[1], k_d[2]

    mismatch_dict = DefaultDict{String, Int}(0)

    for i in 1:(length(dna) - k + 1)
        kmer = dna[i:i+k-1]
        rev_kmer = ReverseComplementDNA(kmer)
        
        # Combine neighbors of original and neighbors of reverse complement
        # We use a Set to ensure we don't double count if kmer == rev_kmer
        for neighbor in get_mismatches(kmer, d)
            mismatch_dict[neighbor] += 1
        end
        for neighbor in get_mismatches(rev_kmer, d)
            mismatch_dict[neighbor] += 1
        end
    end

    max_val = maximum(values(mismatch_dict))
    kmers = [k for (k, v) in mismatch_dict if v == max_val]
    
    result = join(kmers, " ")
    println("1H Result: ", result)
    mkpath("output")
    write("output/Assignment_01H.txt", result)
end

# --- Problem 2A: Protein Translation ---
function solve_2a()
    input_path = "data/stepic_2a.txt"
    if !isfile(input_path) return end
    
    rna = strip(read(input_path, String))
    rna_dict = ProteinDictRNA()
    protein = ""
    
    for i in 1:3:length(rna)-2
        codon = rna[i:i+2]
        amino_acid = rna_dict[codon]
        amino_acid == "Stop" && break
        protein *= amino_acid
    end
    
    println("2A Result: ", protein)
    write("output/Assignment_02A.txt", protein)
end

# Run the solvers
solve_1h()
solve_2a()

