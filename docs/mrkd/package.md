# About MomentTools.jl



## Dependencies

A version of Julia >= 1.3 should be used.

It depends on the following packages, which are installed automatically:

  - `Combinatorics` for the computation of polar ideals.
  - `Dualization` for the use of dual optimization solvers on the moment optimization problems.
  - `DynamicPolynomials` 
  - `JuMP`
  - `LinearAlgebra`
  - `MultivariateSeries` (v>=1.1.2) for the representation of moment sequences
    and the decomposition of series.

To solve the moment optimization problems, SDP optimizers have also to be installed, such as `MosekTools`, `CSD`, ...
    
    
## Development

The git project of the package is
    [https://github.com/AlgebraicGeometricModeling/MomentTools.jl](https://github.com/AlgebraicGeometricModeling/MomentTools.jl).
    
The main developpers are (so far)

  - Lorenzo Baldi
  - Bernard Mourrain

The development is done in relation with the activity of the network [POEMA](http://poema-network.eu/).

## See also 

Other Julia packages are currently developed for polynomial and moment optimization:

  - [`MomentOpt`](https://github.com/lanl-ansi/MomentOpt.jl)
  - [`SumOfSquares`](https://github.com/JuliaOpt/SumOfSquares.jl)
  - [`PolyJuMP`](https://github.com/JuliaOpt/PolyJuMP.jl)
