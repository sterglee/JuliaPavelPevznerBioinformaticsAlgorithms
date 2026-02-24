# Function to calculate Hamming Distance between two strings
function HammingDistance(String1, String2)
    if length(String1) != length(String2)
        return "Error: Lengths do not match."
    end
    
    distance = 0
    # Julia strings are 1-indexed
    for i in 1:length(String1)
        if String1[i] != String2[i]
            distance += 1
        end
    end
    return distance
end

# Function for Approximate Pattern Matching
function ApproximateMatching(Text, Pattern, d)
    t, p = length(Text), length(Pattern)
    indices = Int[] # Define an empty array of Integers
    
    # Loop from 1 to the last possible starting position
    for i in 1:(t - p + 1)
        # Use slicing [i:i+p-1] to get the substring
        distance = HammingDistance(Pattern, Text[i:i + p - 1])
        if distance <= d
            # Subtract 1 to maintain 0-based indexing for the output
            push!(indices, i - 1)
        end
    end
    return indices
end

# --- Execution ---
patt = "ATTCTGGA"
text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC"

result = ApproximateMatching(text, patt, 3)
println("Sample test: Starting indices: $result")

