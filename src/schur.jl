# Realna Schurova forma QR algoritmom (Francisov dvostruki korak s
# implicitnim pomakom), polazeći od Hessenbergove forme:
#
#     T = Uᵀ A U,   U ortogonalna,
#
# T je gornja kvazi-trokutasta: na dijagonali su blokovi 1×1 i 2×2,
# ispod njih same nule.

# POMOCNA FUNKCIJA ZA FRANCIS STEP, PRIMA MATRICE H i Q, indexe hi,lo i vrijednosti s i t (ZBROJ I UMNOŽAK KORIJENA)
function francis_step!(H,Q,hi,lo,s,t)

    # (H - mi1*I)(H - m2*I) = H^2 - s*H + t*I
    h11, h12 = H[lo, lo], H[lo, lo+1]
    h21, h22 = H[lo+1, lo], H[lo+1, lo+1]
    x = h11^2 + h12 * h21 - s * h11 + t
    y = h21 * (h11 + h22 - s)
    z = h21 * H[lo+2, lo+1]


    for k in lo-1:hi-3
        v = householder_vector([x, y, z])
        rows = k+1:k+3
        reflect_left!(view(H, rows, max(lo, k):N), v)
        reflect_right!(view(H, 1:min(k+4, hi), rows), v)
        reflect_right!(view(Q, :, rows), v)

        x = H[k+2, k+1]
        y = H[k+3, k+1]
        if k < hi - 3
            z = H[k+4, k+1]
        end
    end

    # REFLEKSIJA REDA 2
    v = householder_vector([x, y])
    rows = hi-1:hi
    reflect_left!(view(H, rows, hi-2:N), v)
    reflect_right!(view(H, 1:hi, rows), v)
    reflect_right!(view(Q, :, rows), v)
    return H
end

# Vraća gornju kvazi-trokutastu matricu `T` (dijagonalni blokovi reda 1 ili 2)
# i ortogonalnu matricu `U` takve da je `T = Uᵀ A U`.

function real_schur(A::AbstractMatrix; tol = eps(), maxiter = 100)
    H, Q = hessenberg_householder(A)
    n = size(H, 1)
    n == 0 && return H, Q
    normH = opnorm(H, Inf)

    hi = n          # DONJI DESNI KUT NERIJEŠENOG DIJELA
    iter = 0        # BROJ QR KORAKA ZA TRENUTNI BLOK
    while hi >= 1
        # TRAŽIMO LO, POČETAK NOVOG 'PROZORA'
        lo = hi
        while lo > 1 && abs(H[lo, lo-1]) > tol * normH
            lo -= 1
        end
        lo > 1 && (H[lo, lo-1] = 0)

        if lo == hi                  # izdvojen blok 1×1 (realna sv. vrijednost)
            hi -= 1; iter = 0
        elseif lo == hi - 1          # izdvojen blok 2×2
            hi -= 2; iter = 0
        else                         # prozor reda ≥ 3: jedan QR korak
            iter += 1
            iter > maxiter && error("QR algoritam nije konvergirao")

            if iter % 10 == 0
                # "izvanredni" pomak ako konvergencija zapne
                w = abs(H[hi, hi-1]) + abs(H[hi-1, hi-2])
                a = H[hi, hi] + 0.75w
                s, t = 2a, a^2 + 0.4375w^2
            else
                # standardni pomaci: sv. vrijednosti donjeg desnog bloka 2×2
                s = H[hi-1, hi-1] + H[hi, hi]
                t = H[hi-1, hi-1] * H[hi, hi] - H[hi-1, hi] * H[hi, hi-1]
            end
            francis_step!(H, Q, lo, hi, s, t)
        end
        
    end
    # ISPOD PRVE PODDIJAGONALE SU NULE DO NA ZAOKRUŽIVANJE
    for j in 1:n, i in j+2:n
        H[i, j] = 0
    end
    return H, Q
end

function schur_form(A; builtin = false, tol = eps())
    if builtin
        F = schur(float.(Matrix(A)))
        return F.T, F.Z
    end
    return real_schur(A; tol = tol)
end

