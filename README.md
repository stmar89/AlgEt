# AlgEt

Description
--

A package for Magma to make computation with étale algebras.
An étale algebra is described by a direct product of number fields
This guarantees that every operation on elements is very fast.

For the theory on which this code is based, see the `References` section at the bottom.

Please send comments and bug reports to `stefano.marseglia89@gmail.com`.

Details
--

We introduce new type `AlgEtQ`, `AlgEtQElt`, `AlgEtQOrd` and `AlgEtQIdl` which correspond to étale algebras over the rationals, their elements, orders and ideals, respectively.
We have functionalities to compute isomorphism classes of invertible and non-invertible ideals for orders.
For complete descriptions and more details we refer to the [`List of commands for AlgEtQ`](https://github.com/stmar89/AlgEt/blob/main/doc/ListOfCommandsAlgEtQ.md).
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

References
--
Check my webpage for more up-to-date bibliography info.

Stefano Marseglia,<br>
*Computing the ideal class monoid of an order*,<br>
J. Lond. Math. Soc. 101 (2020), no. 3, 984-1007, [`DOI`](https://doi.org/10.1112/jlms.12294)

Stefano Marseglia,<br>
*Cohen-Macaulay type of orders, generators and ideal classes*, [`arXiv`](https://arxiv.org/abs/2206.03758)

Stefano Marseglia,<br>
*Modules over orders, conjugacy classes of integral matrices, and abelian varieties over finite fields*, [`arXiv`](https://arxiv.org/abs/2208.05409)

Stefano Marseglia,<br>
*Local isomorphism classes of fractional ideals of orders in étale algebras*, [`arXiv`](https://arxiv.org/abs/2311.18571)
