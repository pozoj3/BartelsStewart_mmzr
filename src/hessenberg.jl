function householder_vector(x::AbstractVector)
    v = float.(copy(x))
    α = norm(v)
    α == 0 && return v
    v[1] += copysign(α, v[1])
    v ./= norm(v)
    return v
end

reflect_left!(M, v)  = (M .-= 2 .* v * (v' * M); M)
reflect_right!(M, v) = (M .-= 2 .* (M * v) * v'; M)

function hessenberg_householder(A::AbstractMatrix)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("matrica mora biti kvadratna"))

    H = float.(Matrix(A))
    Q = Matrix{eltype(H)}(I, n, n)

    for k in 1:n-2
        v = householder_vector(H[k+1:n, k])

        reflect_left!(view(H, k+1:n, k:n), v)
        reflect_right!(view(H, :, k+1:n), v)
        reflect_right!(view(Q, :, k+1:n), v)

        H[k+2:n, k] .= 0    
    end
    return H, Q
end