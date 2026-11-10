# smallsolve.jl


function solve_small_sylvester(Akk, Bll, D)
    p, q = size(D)
    M = kron(I(q), Akk) + kron(transpose(Bll), I(p))
    x = M \ vec(D)
    return reshape(x, p, q)
end

function solve_small_lyapunov(Akk, D)
    if size(Akk, 1) == 1
        return D ./ (2 .* Akk)
    end

    a = Akk
    M = [2a[1,1]   2a[2,1]          0
         a[1,2]    a[1,1] + a[2,2]  a[2,1]
         0         2a[1,2]          2a[2,2]]
    x = M \ [D[1,1], D[2,1], D[2,2]]

    return [x[1] x[2]
            x[2] x[3]]
end
