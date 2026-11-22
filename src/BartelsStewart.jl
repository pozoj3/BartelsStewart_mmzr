"""
    BartelsStewart

Implementacija Bartels–Stewartova algoritma (Comm. ACM 15(9), 1972.) za

- `AX + XB = C`             → `bartels_stewart(A, B, C)`
- `AᵀX + XA = C`, `C = Cᵀ`  → `lyapunov_bs(A, C)`
"""
module BartelsStewart

using LinearAlgebra

include("hessenberg.jl")   # Householder → gornja Hessenbergova forma
include("schur.jl")        # QR algoritam → realna Schurova forma
include("smallsolve.jl")   # mali sustavi (blokovi reda 1 ili 2)
include("sylvester.jl")    # AX + XB = C
#include("lyapunov.jl")     # AᵀX + XA = C

export bartels_stewart, lyapunov_bs

end 