# A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
# The associated textbook is Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
# The course is run on Coursera and the assignments and textbook are hosted on Stepic.

# Problem Title: Creating a Distance Matrix
# Assignment #: 01
# Problem ID: A 
# URL: https://beta.stepic.org/Bioinformatics-Algorithms-2/Hidden-Messages-in-the-Replication-Origin-2/#step-4

using DelimitedFiles

input_data = readlines("data/stepic_1a.txt")
dna, k = strip.(input_data)
k = parse(Int, k)

kmer_dict = Dict{String, Int}()

for i in 1:(length(dna) - k + 1)
    kmer = dna[i:i+k-1]
    if haskey(kmer_dict, kmer)
        kmer_dict[kmer] += 1
    else
        kmer_dict[kmer] = 1
    end
end

max_count = maximum(values(kmer_dict))
kmers = [kmer for (kmer, count) in kmer_dict if count == max_count]

println(join(kmers, " "))
write("output/Assignment_01A.txt", join(kmers, " "))