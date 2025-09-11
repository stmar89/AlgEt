# AlgEt

Description
--

A package for Magma to make computation with étale algebras.
An étale algebra is described by a direct product of number fields
This guarantees that every operation on elements is very fast.

For the theory on which this code is based, see the `References` section at the bottom.

Please send comments and bug reports to `stefano.marseglia89@gmail.com`.

We introduce new type `AlgEtQ`, `AlgEtQElt`, `AlgEtQOrd` and `AlgEtQIdl` which correspond to étale algebras over the rationals, their elements, orders and ideals, respectively.
There are functionalities for ideal operations, factorizations, trace duality, etc.
There are intrinsics to compute isomorphism classes of invertible (PicardGroup) and non-necessarily-invertible ideals (IdealClassMonoid) for orders.


Instructions
--

The folder `doc` contains a Makefile and magma script to automatically constuct the documentation (from the description of each intrinsic). The output is a MarkDown file `ListOfCommandsAlgEtQ.md` in the same folder.

Each .m file in the folder AlgEtQ contains several quick tests.
These can be combined all together and run using the Makefile in the folder `dev/fast_tests_AlgEtQ_make/`. This test-suit should take only 3-5 minutes.
In the file `dev/slow_tests_AlgEtQ.m`, there is also a series of much slower tests (a few hours), included mainly to compare the timings with previous versions.

To use them, use the magma command `AttachSpec("spec")`, after opening magma in the folder where you have downloaded the repo.

In the folder `examples`, you will find files containing the code to reproduce the examples from the papers in the references below, which should be of help to get a quick start on the functionalities.

References
--

Stefano Marseglia,<br>
*Computing the ideal class monoid of an order*,<br>
J. Lond. Math. Soc. 101 (2020), no. 3, 984-1007, [`DOI`](https://doi.org/10.1112/jlms.12294)

Stefano Marseglia,<br>
*Cohen-Macaulay type of orders, generators and ideal classes*, [`arXiv`](https://arxiv.org/abs/2206.03758)
