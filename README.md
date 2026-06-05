# AlgEt

Description
--

A package for Magma to make computation with étale algebras.
An étale algebra is described by a direct product of number fields.
This guarantees that every operation on elements is very fast.

For the theory on which this code is based, see the `References` section at the bottom.

Please send comments and bug reports to `stefano.marseglia89@gmail.com`.

Details
--

We introduce new types `AlgEtQ`, `AlgEtQElt`, `AlgEtQOrd` and `AlgEtQIdl` which correspond to étale algebras over the rationals, their elements, orders and ideals, respectively.
We have functionalities to compute isomorphism classes of invertible and non-invertible ideals for orders.
For an introduction to the underlying theory, complete descriptions of the intrinsics and examples, we refer to the [`User Manual for AlgEtQ`](https://github.com/stmar89/AlgEt/blob/main/doc/UserManualAlgEtQ.md).
To use them, use the magma command `AttachSpec("spec")`, after opening magma in the folder where you have downloaded the repo.

We also introduce the type `AlgEtQMod`, for modules over `AlgEtQ`.
For complete descriptions and more details we refer to the [`List of commands for AlgEtQMod`](https://github.com/stmar89/AlgEt/blob/main/doc/ListOfCommandsAlgEtQMod.md)
To use them, use the magma command `AttachSpec("specMod")`, after opening magma in the folder where you have downloaded the repo.

We also introduce the type `AlgEtQMtrx`, for matrices over `AlgEtQ`.
For complete descriptions and more details we refer to the [`List of commands for AlgEtQMtrx`](https://github.com/stmar89/AlgEt/blob/main/doc/ListOfCommandsAlgEtQMtrx.md)
To use them, use the magma command `AttachSpec("specMtrx")`, after opening magma in the folder where you have downloaded the repo.

Functionalities for p-adic etale algebras, including how to build completion at rational primes, are developed by Casper Putz. Available [`here`](https://github.com/CPutz/etale-algebra-family).

In the folder [`examples`](https://github.com/stmar89/AlgEt/blob/main/examples), you will find files containing the code to reproduce the examples from the papers in the references below, which should be of help to get a quick start on the functionalities.

<!---
In the file [`examples.txt`](https://github.com/stmar89/PolsAbVarFpCanLift/blob/main/doc/examples.txt) there is the code to see how to use the main functions of the package.
-->

Changelog
--
- v1.0.0 : October 14, 2025 - version submitted for inclusion in Magma.
- v1.0.1 : June 5, 2026.
  - Bugfixes:
    - Added missing declaration of the attribute `PlacesAboveRationalPrime` for `AlgEtQ`.
    - Implemented workaround in `AbsoluteCoordinates` by using `Eltseq` instead of `Flat` to avoid memory leak (which will probably be fixed in the near future).
  - Breaking changes and additional functionalities:
    - The intrinsic `EtaleAlgebra` now can take as input polynomials with rational or integral coefficients.
    - The precision parameter for various intrinsics is now consistently called `Precision`.
    - Privatized internal intrinsics to compute intermediate ideal with underscore prefix (`_MinimalIntermediateIdeals`, etc.).
    - New public `IntermediateIdeals` intrinsic with boolean parameters `Maximal`, `Minimal`, `PrescribedMultiplicatorRing` returning only maximal, minimal and filtered by multiplicator ring intermediate ideals.
    - The intrinsic `Map` to acess the representatives of the abstract representations of the ICM and the weak equivalence class monoid is now called `RepresentativeMap`.

References
--

Stefano Marseglia,<br>
*Computing the ideal class monoid of an order*,<br>
J. Lond. Math. Soc. 101 (2020), no. 3, 984-1007, [`DOI`](https://doi.org/10.1112/jlms.12294) [`arXiv`](https://arxiv.org/abs/1805.09671)

Stefano Marseglia,<br>
*Cohen-Macaulay type of orders, generators and ideal classes*, Journal of Algebra 658 (2024), 247-276. [`DOI`](https://doi.org/10.1016/j.jalgebra.2024.05.051) [`arXiv`](https://arxiv.org/abs/2206.03758)

Stefano Marseglia,<br>
*Modules over orders, conjugacy classes of integral matrices, and abelian varieties over finite fields*, Research in Number Theory 11 (2025) no. 1, paper No. 27. 
Part of the proceedings of the Sixteenth Algorithmic Number Theory Symposium (ANTS XVI).
[`DOI`](https://doi.org/10.1007/s40993-024-00584-9) [`arXiv`](https://arxiv.org/abs/2208.05409)

Stefano Marseglia,<br>
*Local isomorphism classes of fractional ideals of orders in étale algebras*, Journal of Algebra 673 (2025), 77-102. [`DOI`](https://doi.org/10.1016/j.jalgebra.2025.02.030) [`arXiv`](https://arxiv.org/abs/2311.18571)
