# BartelsStewart_mmzr
Projekt za kolegij Matematičko modeliranje i znanstveno računanje - Bartels-Stewartov algoritam za rješavanje Syvesterove jednadžbe
TEST ZA PUSH
## Kostur projekta:
```text
BartelsStewart/
├── Project.toml
├── README.md
├── src/
│   ├── BartelsStewart.jl   # glavni modul
│   ├── hessenberg.jl       # Householder → gornja Hessenbergova forma
│   ├── schur.jl            # QR algoritam → realna Schurova forma
│   ├── smallsolve.jl       # mali sustavi reda 1, 2, 3 i 4
│   ├── sylvester.jl        # AX + XB = C
│   └── lyapunov.jl         # AᵀX + XA = C
└── test/
    └── runtests.jl         # testovi
```

