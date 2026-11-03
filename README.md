# BartelsStewart_mmzr
Projekt za kolegij Matematičko modeliranje i znanstveno računanje - Bartels-Stewartov algoritam za rješavanje Syvesterove jednadžbe

# Kostur projekta:
BartelsStewart/
├── Project.toml
├── src/
│   ├── BartelsStewart.jl   # glavni modul: include-a ostale datoteke, exporta funkcije
│   ├── hessenberg.jl       # Householderovo svođenje na Hessenbergovu formu (HSHLDR + BCKMLT)
│   ├── schur.jl            # QR algoritam → realna Schurova forma (SCHUR)
│   ├── smallsolve.jl       # mali sustavi 1×1, 2×2, 4×4 (SYSSLV)
│   ├── sylvester.jl        # AX + XB = C: transformacija + blok supstitucija (AXPXB + SHRSLV)
│   └── lyapunov.jl         # AᵀX + XA = C, simetrični slučaj (ATXPXA + SYMSLV)
└── test/
    └── runtests.jl         # usporedba s LinearAlgebra.sylvester, provjera reziduala
