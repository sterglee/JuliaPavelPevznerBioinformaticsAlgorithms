"""
Returns the Multiple Longest Common Subsequence and alignment for three strings.
    """
    function multiple_alignment_3(v::AbstractString, w::AbstractString, u::AbstractString)
        n, m, l = length(v), length(w), length(u)

        # Initialize 3D matrices for scores and backtracking
        # S[i, j, k] stores the LCS length
        S = zeros(Int, n + 1, m + 1, l + 1)
        # backtrack stores the move index (0 to 6)
        backtrack = zeros(Int, n + 1, m + 1, l + 1)

        # Fill the 3D Scoring Matrix
        for i in 2:n+1
            for j in 2:m+1
                for k in 2:l+1
                    # Check for a 3-way match
                    match = (v[i-1] == w[j-1] == u[k-1]) ? 1 : 0

                    # Possible moves (following the Python logic indices)
                    # 0: Diagonal (i-1, j-1, k-1)
                    # 1: (i-1, j, k)
                    # 2: (i, j-1, k)
                    # 3: (i, j, k-1)
                    # 4: (i-1, j-1, k)
                    # 5: (i-1, j, k-1)
                    # 6: (i, j-1, k-1)
                    scores = [
                        S[i-1, j-1, k-1] + match,
                        S[i-1, j, k],
                        S[i, j-1, k],
                        S[i, j, k-1],
                        S[i-1, j-1, k],
                        S[i-1, j, k-1],
                        S[i, j-1, k-1]
                        ]

                    val, idx = findmax(scores)
                    S[i, j, k] = val
                    backtrack[i, j, k] = idx - 1 # Store 0-indexed to match Python logic
                end
            end
        end

        # Backtracking
        curr_i, curr_j, curr_k = n + 1, m + 1, l + 1
        max_score = S[curr_i, curr_j, curr_k]

        res_v, res_w, res_u = Char[], Char[], Char[]

        while curr_i > 1 && curr_j > 1 && curr_k > 1
            dir = backtrack[curr_i, curr_j, curr_k]

            if dir == 1
                pushfirst!(res_v, v[curr_i-1]); pushfirst!(res_w, '-'); pushfirst!(res_u, '-')
                curr_i -= 1
                elseif dir == 2
                pushfirst!(res_v, '-'); pushfirst!(res_w, w[curr_j-1]); pushfirst!(res_u, '-')
                curr_j -= 1
                elseif dir == 3
                pushfirst!(res_v, '-'); pushfirst!(res_w, '-'); pushfirst!(res_u, u[curr_k-1])
                curr_k -= 1
                elseif dir == 4
                pushfirst!(res_v, v[curr_i-1]); pushfirst!(res_w, w[curr_j-1]); pushfirst!(res_u, '-')
                curr_i -= 1; curr_j -= 1
                elseif dir == 5
                pushfirst!(res_v, v[curr_i-1]); pushfirst!(res_w, '-'); pushfirst!(res_u, u[curr_k-1])
                curr_i -= 1; curr_k -= 1
                elseif dir == 6
                pushfirst!(res_v, '-'); pushfirst!(res_w, w[curr_j-1]); pushfirst!(res_u, u[curr_k-1])
                curr_j -= 1; curr_k -= 1
            else # dir == 0 (Diagonal)
                pushfirst!(res_v, v[curr_i-1]); pushfirst!(res_w, w[curr_j-1]); pushfirst!(res_u, u[curr_k-1])
                curr_i -= 1; curr_j -= 1; curr_k -= 1
            end
        end

        # Handle remaining prefixes if we hit one wall before (1,1,1)
        while curr_i > 1
            pushfirst!(res_v, v[curr_i-1]); pushfirst!(res_w, '-'); pushfirst!(res_u, '-')
            curr_i -= 1
        end
        while curr_j > 1
            pushfirst!(res_v, '-'); pushfirst!(res_w, w[curr_j-1]); pushfirst!(res_u, '-')
            curr_j -= 1
        end
        while curr_k > 1
            pushfirst!(res_v, '-'); pushfirst!(res_w, '-'); pushfirst!(res_u, u[curr_k-1])
            curr_k -= 1
        end

        return (string(max_score), join(res_v), join(res_w), join(res_u))
    end

    # --- Main Logic ---

    function main()
        input_file = "data/stepic_7g.txt"
        if !isfile(input_file)
            println("Input file not found.")
            return
        end

        lines = readlines(input_file)
        word1, word2, word3 = strip(lines[1]), strip(lines[2]), strip(lines[3])

        alignment = multiple_alignment_3(word1, word2, word3)

        println(join(alignment, "\n"))

        mkpath("output")
        write("output/Assignment_07G.txt", join(alignment, "\n"))
    end

    main()
