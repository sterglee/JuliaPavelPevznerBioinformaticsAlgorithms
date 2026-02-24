using Statistics

"""
dpChange(money::Int, coins::Vector{Int}) -> Int

Solve the Change Problem: find the minimum number of coins needed to make `money` using denominations in `coins`.
"""
function dpChange(money::Int, coins::Vector{Int})::Int
    minNumCoins = getMinNumCoins(money, coins)
    return minNumCoins[money + 1]
end

function getMinNumCoins(money::Int, coins::Vector{Int})::Vector{Int}
    minNumCoins = fill(typemax(Int), money + 1)
    minNumCoins[1] = 0  # minNumCoins(0) in Scala is index 1 in Julia
    for m in 1:money
        for coin in coins
            if m >= coin && minNumCoins[m - coin + 1] + 1 < minNumCoins[m + 1]
                minNumCoins[m + 1] = minNumCoins[m - coin + 1] + 1
            end
        end
    end
    return minNumCoins
end

"""
longestPath(n::Int, m::Int, down::Vector{Vector{Int}}, right::Vector{Vector{Int}}) -> Int

Find the length of the longest path in the Manhattan Tourist Problem.
"""
function longestPath(n::Int, m::Int, down::Vector{Vector{Int}}, right::Vector{Vector{Int}})::Int
    s = zeros(Int, n+1, m+1)

    for i in 2:n+1
        s[i, 1] = s[i-1, 1] + down[i-1-1][1-1]  # down(i-1)(0) in Scala is down[i-2][0] in Julia
    end
    for j in 2:m+1
        s[1, j] = s[1, j-1] + right[1-1][j-1-1]
    end

    for i in 2:n+1
        for j in 2:m+1
            s[i, j] = max(s[i-1, j] + down[i-1-1][j-1], s[i, j-1] + right[i-1][j-1-1])
        end
    end

    return s[n+1, m+1]
end
