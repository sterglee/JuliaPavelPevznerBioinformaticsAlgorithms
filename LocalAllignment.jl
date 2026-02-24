"""
Bioinformatics Algorithms: Local Alignment (Smith-Waterman)
"""

# Helper to provide the PAM250 scoring matrix
# You can expand this dictionary with the full PAM250 values
function get_pam250_dict()
    # This is a subset for demonstration. 
    # In a real scenario, you'd parse the full 20x20 matrix.
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    # Example values (Match = 2, Mismatch = -1 for simplicity in this snippet)
    # Replace this logic with your actual PAM250 data loading
    pam = Dict{Tuple{Char, Char}, Int}()
    
    # Example: If you have the matrix as a string or file, parse it here.
    # For now, let's assume a simple match/mismatch if the key is missing:
    return pam
end

# Robust scoring function to handle missing keys and symmetry
function get_score(matrix, char_v::Char, char_w::Char)
    if haskey(matrix, (char_v, char_w))
        return matrix[(char_v, char_w)]
    elseif haskey(matrix, (char_w, char_v))
        return matrix[(char_w, char_v)]
    else
        # Standard PAM250 default for mismatches if not in dict
        return char_v == char_w ? 5 : -4 
    end
end

function local_alignment(v::AbstractString, w::AbstractString, scoring_matrix, sigma::Int)
    n, m = length(v), length(w)
    
    # S[i, j] stores scores, backtrack[i, j] stores move indices
    S = zeros(Int, n + 1, m + 1)
    backtrack = zeros(Int, n + 1, m + 1)

    # 1. Fill the Scoring Matrix
    for i in 2:n+1
        for j in 2:m+1
            # 1: Up, 2: Left, 3: Diagonal, 4: Restart (0)
            match_score = get_score(scoring_matrix, v[i-1], w[j-1])
            
            scores = [
                S[i-1, j] - sigma, 
                S[i, j-1] - sigma, 
                S[i-1, j-1] + match_score, 
                0
            ]
            
            val, idx = findmax(scores)
            S[i, j] = val
            backtrack[i, j] = idx
        end
    end

    # 2. Find the peak (Starting point for local alignment)
    max_score, max_idx = findmax(S)
    curr_i, curr_j = Tuple(max_idx)

    # 3. Backtrack to find the alignment
    res_v = Char[]
    res_w = Char[]

    while curr_i > 1 && curr_j > 1 && backtrack[curr_i, curr_j] != 4
        direction = backtrack[curr_i, curr_j]
        
        if direction == 1      # Up (Gap in W)
            pushfirst!(res_v, v[curr_i-1])
            pushfirst!(res_w, '-')
            curr_i -= 1
        elseif direction == 2  # Left (Gap in V)
            pushfirst!(res_v, '-')
            pushfirst!(res_w, w[curr_j-1])
            curr_j -= 1
        elseif direction == 3  # Diagonal (Match/Mismatch)
            pushfirst!(res_v, v[curr_i-1])
            pushfirst!(res_w, w[curr_j-1])
            curr_i -= 1
            curr_j -= 1
        end
    end

    return (string(max_score), join(res_v), join(res_w))
end

# --- Main Logic ---

# 1. Initialize the scoring matrix
PAM250_DICT = get_pam250_dict()

# 2. Define or Load your sequences
 word1 = "AACTAAGGTT"
 word2 = "ATAAGGT"

# 3. Run the alignment
 alignment = local_alignment(word1, word2, PAM250_DICT, 5)
 println(join(alignment, "\n"))