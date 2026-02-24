"""
Hirschberg's Algorithm: Space-Efficient Global Alignment.
Constructs the optimal alignment using O(min(n, m)) space for scoring.
    """

    # --- 1. Scoring Matrix Parser ---
    function parse_blosum(filepath::String)
        if !isfile(filepath)
            error("Scoring matrix file not found at: $filepath")
        end
        lines = filter(l -> !startswith(l, "#") && !isempty(strip(l)), readlines(filepath))
        amino_acids = [c[1] for c in split(lines[1]) if length(c) == 1]
            matrix_dict = Dict{Tuple{Char, Char}, Int}()

            for i in 2:length(lines)
                parts = split(lines[i])
                row_char = uppercase(parts[1][1])
                col_idx = 1
                for j in 2:length(parts)
                    val = tryparse(Int, parts[j])
                    if val !== nothing && col_idx <= length(amino_acids)
                        col_char = uppercase(amino_acids[col_idx])
                        matrix_dict[(row_char, col_char)] = val
                        col_idx += 1
                    end
                end
            end
            return matrix_dict
        end

        # --- 2. Main Alignment Logic ---
        function space_efficient_global_alignment(v::AbstractString, w::AbstractString, scoring_matrix::Dict, sigma::Int)

            # Helper for Linear Space scoring
            function middle_column_score(seq1::AbstractString, seq2::AbstractString)
                n = length(seq1)
                prev = [-(i-1) * sigma for i in 1:n+1]
                    curr = zeros(Int, n + 1)

                    for j in 1:length(seq2)
                        curr[1] = -j * sigma
                        for i in 2:n+1
                            char1, char2 = seq1[i-1], seq2[j]
                            match_score = get(scoring_matrix, (uppercase(char1), uppercase(char2)), -4)

                            match = prev[i-1] + match_score
                            deletion = curr[i-1] - sigma
                            insertion = prev[i] - sigma
                            curr[i] = max(match, deletion, insertion)
                        end
                        prev = copy(curr)
                    end
                    return prev
                end

                # Needleman-Wunsch for small sub-problems (Base Case)
                function small_nw(v_s::AbstractString, w_s::AbstractString)
                    n, m = length(v_s), length(w_s)
                    S = zeros(Int, n+1, m+1)
                    for i in 2:n+1 S[i,1] = S[i-1,1] - sigma end
                    for j in 2:m+1 S[1,j] = S[1,j-1] - sigma end

                    for i in 2:n+1, j in 2:m+1
                        score = get(scoring_matrix, (uppercase(v_s[i-1]), uppercase(w_s[j-1])), -4)
                        S[i,j] = max(S[i-1,j-1] + score, S[i-1,j] - sigma, S[i,j-1] - sigma)
                    end

                    res_v, res_w = Char[], Char[]
                    i, j = n+1, m+1
                    while i > 1 || j > 1
                        if i > 1 && j > 1 && S[i,j] == S[i-1,j-1] + get(scoring_matrix, (uppercase(v_s[i-1]), uppercase(w_s[j-1])), -4)
                            pushfirst!(res_v, v_s[i-1]); pushfirst!(res_w, w_s[j-1]); i-=1; j-=1
                            elseif i > 1 && S[i,j] == S[i-1,j] - sigma
                            pushfirst!(res_v, v_s[i-1]); pushfirst!(res_w, '-'); i-=1
                        else
                            pushfirst!(res_v, '-'); pushfirst!(res_w, w_s[j-1]); j-=1
                        end
                    end
                    return [join(res_v), join(res_w)]
                end

                # Find Middle Edge
                function get_middle_edge(top, bottom, left, right)
                    mid_col = floor(Int, (left + right) / 2)

                    pref = middle_column_score(v[top+1:bottom], w[left+1:mid_col])
                    suff = middle_column_score(reverse(v[top+1:bottom]), reverse(w[mid_col+2:right]))

                    combined = pref .+ reverse(suff)
                    row_idx = argmax(combined) - 1
                    mid_node = (top + row_idx, mid_col)

                    # Determine next node (Edge direction)
                    if mid_node[1] < bottom && mid_node[2] < right
                        diag = pref[row_idx+1] + get(scoring_matrix, (uppercase(v[mid_node[1]+1]), uppercase(w[mid_node[2]+1])), -4) + reverse(suff)[row_idx+2]
                        down = pref[row_idx+2] - sigma + reverse(suff)[row_idx+2]

                        if combined[row_idx+1] == diag
                            next_node = (mid_node[1] + 1, mid_node[2] + 1)
                            elseif combined[row_idx+1] == down
                            next_node = (mid_node[1] + 1, mid_node[2])
                        else
                            next_node = (mid_node[1], mid_node[2] + 1)
                        end
                        elseif mid_node[1] < bottom
                        next_node = (mid_node[1] + 1, mid_node[2])
                    else
                        next_node = (mid_node[1], mid_node[2] + 1)
                    end
                    return mid_node, next_node
                end

                # Recursive D&C
                function linear_space_alignment2(top, bottom, left, right)
                    if left == right
                        return [v[top+1:bottom], "-"^(bottom - top)]
                        elseif top == bottom
                        return ["-"^(right - left), w[left+1:right]]
                        elseif bottom - top == 1 || right - left == 1
                        return small_nw(v[top+1:bottom], w[left+1:right])
                    else
                        mid, nxt = get_middle_edge(top, bottom, left, right)

                        char_v = nxt[1] > mid[1] ? v[nxt[1]] : '-'
                        char_w = nxt[2] > mid[2] ? w[nxt[2]] : '-'

                        A = linear_space_alignment2(top, mid[1], left, mid[2])
                        B = linear_space_alignment2(nxt[1], bottom, nxt[2], right)

                        return [A[1] * char_v * B[1], A[2] * char_w * B[2]]
                    end
                end

                # Execution
                res = linear_space_alignment2(0, length(v), 0, length(w))

                # Calculate Score
                final_score = 0
                for (a, b) in zip(res[1], res[2])
                    if a == '-' || b == '-'
                        final_score -= sigma
                    else
                        final_score += get(scoring_matrix, (uppercase(a), uppercase(b)), -4)
                    end
                end
                return [string(final_score), res[1], res[2]]
            end

            # --- 3. Main Routine ---
            function main()
                # Ensure these files exist in your working directory
                input_path = "data/stepic_7f.txt"
                matrix_path = "scripts/BLOSUM62.txt"

                if !isfile(input_path) || !isfile(matrix_path)
                    println("Error: Required files not found.")
                    return
                end

                lines = readlines(input_path)
                # Using String() ensures we aren't passing SubStrings if dispatch is picky
                v = String(strip(lines[1]))
                w = String(strip(lines[2]))

                println("Loading BLOSUM62...")
                blosum62 = parse_blosum(matrix_path)

                println("Running Hirschberg Alignment...")
                @time result = space_efficient_global_alignment(v, w, blosum62, 5)

                # Print and Save
                println("\nScore: ", result[1])
                mkpath("output")
                write("output/Assignment_07F.txt", join(result, "\n"))
                println("Results written to output/Assignment_07F.txt")
            end

            # Start the program
            main()

