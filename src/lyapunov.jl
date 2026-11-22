
function solve_quasitriangular_lyapunov!(Ap, Y)
    r = block_ranges(Ap)
    At = Matrix(transpose(Ap))       

    for l in eachindex(r), k in l:length(r)   
        ik, jl = r[k], r[l]
        above = 1:first(ik)-1
        left  = 1:first(jl)-1

        D = Y[ik, jl] - At[ik, above] * Y[above, jl] - Y[ik, left] * Ap[left, jl]

        if k == l
            Y[ik, ik] = solve_small_lyapunov(Ap[ik, ik], (D + D') / 2)
        else
            Y[ik, jl] = solve_small_sylvester(At[ik, ik], Ap[jl, jl], D)
            Y[jl, ik] = transpose(Y[ik, jl])   
        end
    end
    return Y
end
