# sylvester.jl
#
# Bartels–Stewartov algoritam za AX + XB = C:
#
#   1. A' = Uᵀ A U  (donja forma),  B' = Vᵀ B V  (gornja forma)
#   2. C' = Uᵀ C V
#   3. riješi A'X' + X'B' = C' blok po blok
#   4. X  = U X' Vᵀ

function block_ranges(T; lower = false)
    n = size(T, 1)
    ranges = UnitRange{Int}[]
    i = 1
    while i <= n
        coupled = i < n && (lower ? T[i, i+1] : T[i+1, i]) != 0
        if coupled
            push!(ranges, i:i+1); i += 2
        else
            push!(ranges, i:i);   i += 1
        end
    end
    return ranges
end


function solve_quasitriangular_sylvester!(Al, Bu, Y)
    rA = block_ranges(Al, lower = true)
    rB = block_ranges(Bu)

    for jl in rB, ik in rA
        above = 1:first(ik) - 1
        left = 1:first(jl) - 1
        D = Y[ik, jl] - Al[ik, above] * Y[above, jl] - Y[ik, left] * Bu[left, jl]
        Y[ik, jl] = solve_small_sylvester(Al[ik, ik], Bu[jl, jl], D)
    end
    return Y
end

function bartels_stewart(A, B, C; builtin_schur = false, tol = eps())
    m, n = size(C)
    size(A) == (m, m) || throw(DimensionMismatch("A mora biti $(m)×$(m)"))
    size(B) == (n, n) || throw(DimensionMismatch("B mora biti $(n)×$(n)"))


    TAt, U = schur_form(transpose(A); builtin = builtin_schur, tol = tol)
    Al = Matrix(transpose(TAt))
    Bu, V = schur_form(B; builtin = builtin_schur, tol = tol)

    Y = U' * C * V

    solve_quasitriangular_sylvester!(Al, Bu, Y)

    return U * Y * V'
end