## List of instrinsics in AlgEtQ/AlgEt.m:

> <pre><b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ</pre>
---
*Given a sequence of number fields returns the étale algebra corresponding to the direct product.*

> <pre><b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ</pre>
---
*Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.*

> <pre><b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ</pre>
---
*Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.*


## List of instrinsics in AlgEtQ/AlgEtAttributes.m:

> <pre><b>Print</b>(A::AlgEtQ)</pre>
---
*Prints the defining polynomial or the components defining A.*

> <pre><b>DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt</pre>
---
*Returns the defining polynomial of A, if the corresponding number fields are distinct.*

> <pre><b>Components</b>(A::AlgEtQ) -> SeqEnum</pre>
---
*Returns the number fields of which A is a product of,together with embeddings and projections.*

> <pre><b>Dimension</b>(A::AlgEtQ)->RngInt</pre>
---
*Dimension of A.*

> <pre><b>AbsoluteDimension</b>(A::AlgEtQ)->RngInt</pre>
---
*Dimension of A over the prime field.*

> <pre><b>HasBaseField</b>(A::AlgEtQ) -> BoolElt,FldNum</pre>
---
*Returns whether A has common base field. If this is the case it returns it.*

> <pre><b>BaseField</b>(A::AlgEtQ) -> FldNum</pre>
---
*Returns the common base field of the Algebra, if it exists.*

> <pre><b>PrimeField</b>(A::AlgEtQ) -> FldNum</pre>
---
*Returns the prime field of the Algebra.*

> <pre><b>'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt</pre>
---
*A1 eq A2.*


## List of instrinsics in AlgEtQ/Homs.m:

> <pre><b>HomsToC</b>(A::AlgEtQ : Precision:=30)->SeqEnum[Map]</pre>
---
*returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30*


## List of instrinsics in AlgEtQ/Elt.m:

> <pre><b>Print</b>(x::AlgEtQElt)</pre>
---
*Print the AlgEtQElt.*

> <pre><b>Parent</b>(x::AlgEtQElt) -> AlgEtQ</pre>
---
*Returns the algebra to which the elemenet belongs to.*

> <pre><b>Algebra</b>(x::AlgEtQElt) -> AlgEtQ</pre>
---
*Returns the algebra to which the elemenet belongs to.*

> <pre><b>Components</b>(x::AlgEtQElt) -> SeqEnum</pre>
---
*Given an element x returns its components, which are elements of number fields.*

> <pre><b>AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum</pre>
---
*Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.*

> <pre><b>IsCoercible</b>(A::AlgEtQ, x::.) -> BoolElt, .</pre>
---
*Return whether x is coercible into A and the result of the coercion if so.*

> <pre><b>'!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt</pre>
---
*Coerce x into A.*

> <pre><b>One</b>(A::AlgEtQ) -> AlgEtQElt</pre>
---
*The multiplicative neutral element of A.*

> <pre><b>Zero</b>(A::AlgEtQ) -> AlgEtQElt</pre>
---
*The additive neutral element of A.*

> <pre><b>IsUnit</b>(x::AlgEtQElt) -> BoolElt</pre>
---
*Returns wheter x is a unit in A.*

> <pre><b>IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt</pre>
---
*Returns wheter x is a not unit in A.*

> <pre><b>IsZeroDivisor2</b>(x::AlgEtQElt) -> BoolElt</pre>
---
*Returns wheter x is a not unit in A.*

> <pre><b>Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre>
---
*Random element of A. The Coefficients are bounded by the positive integer bd.*

> <pre><b>Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt</pre>
---
*Random element of A. The Coefficients are bounded by VarArg bd (default 3).*

> <pre><b>RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre>
---
*Random unit of A. The Coefficients are bounded by the positive integer bd.*

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt</pre>
---
*Is x1=x2 ?*

> <pre><b>'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt</pre>
---
*Is x1=x2 ?*

> <pre><b>'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt</pre>
---
*Is x1=x2 ?*

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt</pre>
---
*Is x1=x2 ?*

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt</pre>
---
*Is x1=x2 ?*

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
---
*x1+x2.*

> <pre><b>'-'</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
---
*-x.*

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
---
*x1-x2.*

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
---
*x1\*x2.*

> <pre><b>Inverse</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
---
*1/x.*

> <pre><b>'^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt</pre>
---
*x^n.*

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
---
*x1/x2.*

> <pre><b>'&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
---
*Given a sequence of AlgEtQElt returns the sum of the entries.*

> <pre><b>'&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
---
*Given a sequence of AlgEtQElt returns the product of the entries.*

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
---
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

> <pre><b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
---
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

> <pre><b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
---
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt</pre>
---
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt</pre>
---
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

> <pre><b>MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
---
*Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.*

> <pre><b>MinimalPolynomial</b>(x::AlgEtQElt, F::Rng) -> RngUPolElt</pre>
---
*Returns the minimal polynommial over the ring F of the element x.*

> <pre><b>AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
---
*Returns the minimal polynommial over the prime field of the element x.*

> <pre><b>IsIntegral</b>(x::AlgEtQElt) -> BoolElt</pre>
---
*Returns whether the element x is integral (over the integers).*

> <pre><b>Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt</pre>
---
*Evaluate the polynomial f at the element a.*

> <pre><b>PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt</pre>
---
*Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.*

> <pre><b>PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]</pre>
---
*Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.*

> <pre><b>Basis</b>(A::AlgEtQ) -> SeqEnum</pre>
---
*Returns a basis of the algebra over the common base field.*

> <pre><b>AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum</pre>
---
*Returns a basis of the algebra over the prime field.*

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum</pre>
---
*Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.*

> <pre><b>OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
---
*Returns the orthogonal ideampotent element of the étale algebra A.*

> <pre><b>Idempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
---
*Returns the ideampotent element of the étale algebra A.*


## List of instrinsics in AlgEtQ/TraceNorm.m:

> <pre><b>Trace</b>(x::AlgEtQElt) -> Any</pre>
---
*Returns the trace of the element x of an étale algebra.*

> <pre><b>Norm</b>(x::AlgEtQElt) -> Any</pre>
---
*Returns the norm of the element x of an étale algebra.*

> <pre><b>AbsoluteTrace</b>(x::AlgEtQElt) -> Any</pre>
---
*Returns the absolute trace of the element x of an étale algebra.*

> <pre><b>AbsoluteNorm</b>(x::AlgEtQElt) -> Any</pre>
---
*Returns the absolute norm of the element x of an étale algebra.*

> <pre><b>TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Returns the trace dual ideal of an ideal in an order in an etale algebra.*

> <pre><b>TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Returns the trace dual ideal of an order in an etale algebra.*


## List of instrinsics in AlgEtQ/Ord.m:

> <pre><b>Print</b>(A::AlgEtQOrd)</pre>
---
*Print the order.*

> <pre><b>IsCoercible</b>(S::AlgEtQOrd, x::.) -> BoolElt, Any</pre>
---
*Return whether x is coercible into S and the result if so.*

> <pre><b>Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd</pre>
---
*Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped. The vararg CheckIsKnownOrder determines if we check if the order is already known, i.e. in the attribute Algebra`KnownOrders. The default value is true.*

> <pre><b>Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd</pre>
---
*Given a sequence of order in the number fiedls defining the etale algebra A, generates the product order.*

> <pre><b>Algebra</b>(S::AlgEtQOrd) -> AlgEtQ</pre>
---
*Returns the algebra of the order.*

> <pre><b>myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]</pre>
---
*Hash function for AlgEtQOrd.*

> <pre><b>ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
---
*Return a Z-basis of the order.*

> <pre><b>Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
---
*Return a set of generators of the order.*

> <pre><b>'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
---
*Checks equality of orders in an etale Algebra.*

> <pre><b>'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt</pre>
---
*Inclusion of elements.*

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum</pre>
---
*AbsoluteCoordinates with respect to the ZBasis.*

> <pre><b>'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt</pre>
---
*Inclusion of elements.*

> <pre><b>'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt</pre>
---
*Inclusion of elements.*

> <pre><b>One</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
---
*Unit element of S.*

> <pre><b>Zero</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
---
*Zero element of S.*

> <pre><b>Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
---
*Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

> <pre><b>Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
---
*Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

> <pre><b>IsKnownOrder</b>(~R::AlgEtQOrd)</pre>
---
*This procedure checks wheter the order R is already in the list of known orders of the algebra A of definition of R. If so then it replaces R with the copy stored in the attribute KnownOrders. If not it adds it to KnownOrders. This is done to avoid creating multiple copies of the same order.*

> <pre><b>EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd</pre>
---
*Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial.*

> <pre><b>ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd</pre>
---
*Given a product of number field A, returns the order consisting of the product of the equation orders of the number fields.*

> <pre><b>MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd</pre>
---
*Returns the maximal order of the étale algebra A.*

> <pre><b>IsMaximal</b>(S::AlgEtQOrd) -> BoolElt</pre>
---
*Returns wheter the given order is the maximal order of the étale algebra.*

> <pre><b>IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup</pre>
---
*Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.*

> <pre><b>Index</b>(T::AlgEtQOrd) -> FldRatElt</pre>
---
*Given an order T computes its index with respect to the basis of the algebra of T as a free Z-module.*

> <pre><b>Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any</pre>
---
*Given two orders T \subset S, returns [S:T] = #S/T.*

> <pre><b>'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt</pre>
---
*Checks if the first argument is inside the second.*

> <pre><b>'*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
---
*Returns the order generated by the orders O1 and O2.*

> <pre><b>'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
---
*Intersection of orders.*

> <pre><b>MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd</pre>
---
*Returns the multiplicator ring of an order R, that is R itself.*


## List of instrinsics in AlgEtQ/Quotients.m:

> <pre><b>Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
---
*Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order.*

> <pre><b>Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map</pre>
---
*Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.*

> <pre><b>Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
---
*Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.*

> <pre><b>ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map</pre>
---
*Given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S). We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.*

> <pre><b>ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map</pre>
---
*Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.*

> <pre><b>PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt</pre>
---
*Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^\*.*

> <pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
---
*Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

> <pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
---
*Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

> <pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
---
*Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

> <pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
---
*Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*


## List of instrinsics in AlgEtQ/OverOrders.m:

> <pre><b>IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt</pre>
---
*Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.*

> <pre><b>MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]</pre>
---
*Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

> <pre><b>MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]</pre>
---
*Computes the minimal overorders of R.*

> <pre><b>OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]</pre>
---
*Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

> <pre><b>OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]</pre>
---
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*

> <pre><b>FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]</pre>
---
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*


## List of instrinsics in AlgEtQ/GraphOverOrders.m:

> <pre><b>GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir</pre>
---
*Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].*


## List of instrinsics in AlgEtQ/Idl.m:

> <pre><b>Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gens.*

> <pre><b>Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl</pre>
---
*Given an order S which is a product of orders S_i in the number fiedls generting the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal direct product of the I_i.*

> <pre><b>Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>'*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Creates an ideal of S, generated by gen.*

> <pre><b>Print</b>(I::AlgEtQIdl)</pre>
---
*Prints the ideal.*

> <pre><b>'!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is in S, then IT=I.*

> <pre><b>Algebra</b>(I::AlgEtQIdl) -> AlgEtQ</pre>
---
*Returns the étale algebra in which the ideal lives.*

> <pre><b>Order</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
---
*Returns the order of definition of the ideal.*

> <pre><b>ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]</pre>
---
*Returns a Z-basis of the ideal.*

> <pre><b>Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre>
---
*Returns the generators of the ideal.*

> <pre><b>myHash</b>(I::AlgEtQIdl)->RngInt</pre>
---
*Hash function.*

> <pre><b>'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre>
---
*Equality testing.*

> <pre><b>'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre>
---
*Equality testing.*

> <pre><b>'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt</pre>
---
*Return if I eq S. I needs to be an ideal of S.*

> <pre><b>'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
---
*Return if I eq S. I needs to be an ideal of S.*

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum</pre>
---
*AbsoluteCoordiantes with respect to the ZBasis.*

> <pre><b>'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt</pre>
---
*Returns if x is in I.*

> <pre><b>'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt</pre>
---
*Returns if x is in I.*

> <pre><b>'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt</pre>
---
*Returns if x is in I.*

> <pre><b>'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
---
*Given an ideal I of S, return if S subseteq I.*

> <pre><b>'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt</pre>
---
*Given an ideal I of S, return if I subseteq S.*

> <pre><b>'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt</pre>
---
*Checks if the first argument is inside the second. The ideals need to be fractional.*

> <pre><b>Index</b>(T::AlgEtQIdl) -> FldRatElt</pre>
---
*Given an ideal T computes its index with respect to the basis of the algebra of T as a free Z-module.*

> <pre><b>Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any</pre>
---
*Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I].*

> <pre><b>Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any</pre>
---
*Given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I].*

> <pre><b>OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Given an order S returns the ideal 1\*S which will be cached.*

> <pre><b>Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl</pre>
---
*Computes the conductor of an order, defined as he colon ideal (O:OK), where OK is the maximal order of the algebra.*

> <pre><b>'+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
---
*Returns the sum of two ideals.*

> <pre><b>'*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
---
*Product of two ideals.*

> <pre><b>'*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl</pre>
---
*Returns x\*I.*

> <pre><b>'^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl</pre>
---
*nth power of an ideal.*

> <pre><b>'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl</pre>
---
*Given an ideal I of S, return S cap I.*

> <pre><b>'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Given an ideal I of S, return S cap I.*

> <pre><b>'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*Given ideals I and J, return J cap I.*

> <pre><b>'&+'</b>(seq::SeqEnum[AlgEtQIdl])->AlgEtQIdl</pre>
---
*Returns the sum of the fractional ideals in the sequence.*

> <pre><b>ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl</pre>
---
*Computes the colon ideal (I:J) (as an O-ideal) of two O-idealsi.*

> <pre><b>ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl</pre>
---
*Computes the colon ideal (1\*O:J) (as an O-ideal).*

> <pre><b>ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl</pre>
---
*Computes the colon ideal (I:1\*O) (as an O-ideal).*

> <pre><b>IsInvertible</b>(I::AlgEtQIdl) ->BoolElt</pre>
---
*Checks if the ideal I is invertible in its order of definition O.*

> <pre><b>Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl</pre>
---
*Computes the inverse of an ideal of a maximal order.*

> <pre><b>MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
---
*Given a fractional R-ideal I computes its multiplicator ring (I:I). If the overorders of R are known the corresponding overorder is returned, in order to preserve the known attributes.*

> <pre><b>IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup</pre>
---
*Return if the argument is a product of ideals in number fields, and if so return also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).*

> <pre><b>Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
---
*Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

> <pre><b>Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
---
*Returns a random (small coefficient) element of I. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

> <pre><b>IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt</pre>
---
*Given two integral ideals I and J of an order S, returns whether I+J=R.*

> <pre><b>IsIntegral</b>(I::AlgEtQIdl) -> BoolElt</pre>
---
*Returns wheter the ideal I of S is integral, that is I \subseteq S.*

> <pre><b>MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt</pre>
---
*Given a fractional S ideal I, returns the ideal d\*I,d when d is the smallest integer such that d\*I is integral in S.*

> <pre><b>MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt</pre>
---
*Returns the smallest integer contained in the ideal I.*

> <pre><b>CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl</pre>
---
*Returns an element x such that x\*I is an integral ideal coprime with J, togheter with the product x\*I. The first ideal must be invertible and the second should be integral.*


## List of instrinsics in AlgEtQ/ZBasisLLL.m:

> <pre><b>ZBasisLLL</b>(S::AlgEtQOrd)</pre>
---
*A procedure that replaces the ZBasis with an LLL-reduced one.*

> <pre><b>ZBasisLLL</b>(S::AlgEtQIdl)</pre>
---
*A procedure that replaces the ZBasis with an LLL-reduced one.*


## List of instrinsics in AlgEtQ/Completion.m:

> <pre><b>Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map</pre>
---
*Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage (it acts a bit weird).*


## List of instrinsics in AlgEtQ/ComplexConj.m:

> <pre><b>HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt</pre>
---
*Returns if the algebra is the product of CM fields.*

> <pre><b>ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
---
*If A is a product of CM fields, it returns the complex conjugate of the argument.*

> <pre><b>IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd</pre>
---
*Returns wheter O is conjugate stable and the complex conjugate.*

> <pre><b>ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd</pre>
---
*It returns the complex conjugate of the argument.*

> <pre><b>IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl</pre>
---
*Returns wheter O is conjugate stable and the complex conjugate.*

> <pre><b>ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
---
*If A is a product of CM fields, it returns the complex conjugate of the argument.*


## List of instrinsics in AlgEtQ/ComplexMult.m:

> <pre><b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
---
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

> <pre><b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
---
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

> <pre><b>CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
---
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

> <pre><b>CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
---
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

> <pre><b>Print</b>( PHI :: AlgEtQCMType)</pre>
---
*Print the AlgEtQCMType.*

> <pre><b>CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
---
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

> <pre><b>CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
---
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

> <pre><b>Homs</b>( PHI::AlgEtQCMType : prec:=30 )->SeqEnum[Map]</pre>
---
*Given a AlgEtQCMType PHI returns the sequence of maps to CC defining it.*

> <pre><b>'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=30)->BoolElt</pre>
---
*Returns whether two cm types are equal. This happens if and only if the ration of (any) two CMPositiveElements is totally real and totally positive.*

> <pre><b>Precision</b>(PHI :: AlgEtQCMType)->RngIntElt</pre>
---
*Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

> <pre><b>ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType</pre>
---
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

> <pre><b>ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )</pre>
---
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

> <pre><b>AllCMTypes</b>(A::AlgEtQ : Precision := 30 ) -> SeqEnum[AlgEtQCMType]</pre>
---
*Returns all the AlgEtQCMTypes of A.*


## List of instrinsics in AlgEtQ/IntermediateIdeals.m:

> <pre><b>MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

> <pre><b>IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.*

> <pre><b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. They are produced recursively using from the minimal ones.*

> <pre><b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

> <pre><b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

> <pre><b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

> <pre><b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

> <pre><b>IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]</pre>
---
*Given ideals J subset I over the same order, and a positive integer N, it returns all the ideals K such that J subset K subset I and [I:K]=N. These are produced by recursively searching for maximal submodules.*

> <pre><b>MinimalIntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl : primes:=[])->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I.*

> <pre><b>IntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
---
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones*


## List of instrinsics in AlgEtQ/IdealsOfIndex.m:

> <pre><b>IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
---
*Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.*

> <pre><b>IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
---
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

> <pre><b>IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]</pre>
---
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

> <pre><b>IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
---
*Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".*

> <pre><b>IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
---
*Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".*


## List of instrinsics in AlgEtQ/ShortEltSmallRep.m:

> <pre><b>ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt</pre>
---
*Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).*

> <pre><b>SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt</pre>
---
*Given a fractional R-ideal I, it returns an isomorphic ideal a\*I, and the element a, such that a\*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.*


## List of instrinsics in AlgEtQ/MinimalGenerators.m:

> <pre><b>TwoGeneratingSet</b>(I::AlgEtQIdl)</pre>
---
*A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.*


## List of instrinsics in AlgEtQ/CRT.m:

> <pre><b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt</pre>
---
*Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.*

> <pre><b>ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt</pre>
---
*Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.*


## List of instrinsics in AlgEtQ/PicardGroup.m:

> <pre><b>ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map</pre>
---
*Returns the group (S/I)^\* and a map (S/I)^\* -> S. It is required S to be maximal.*

> <pre><b>IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgAssElt</pre>
---
*Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*

> <pre><b>PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
---
*Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".*

> <pre><b>UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
---
*Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".*

> <pre><b>IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgAssElt</pre>
---
*Checks if I=x\*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*


## List of instrinsics in AlgEtQ/FactPrimes.m:

> <pre><b>Factorization</b>(I::AlgEtQIdl) -> Tup</pre>
---
*Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.*

> <pre><b>PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgAssEtOrdIdl]</pre>
---
*Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.*

> <pre><b>SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgAssEtOrdIdl]</pre>
---
*Returns the non-invertible primes of the order.*

> <pre><b>NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx</pre>
---
*Returns the non-invertible primes of the order.*

> <pre><b>IsPrime</b>(I::AlgEtQIdl) -> BoolElt</pre>
---
*Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal.*

> <pre><b>IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
---
*Check if the order is Bass at the prime ideal P.*

> <pre><b>IsBass</b>(S::AlgEtQOrd) -> BoolElt</pre>
---
*Check if the order is Bass.*

> <pre><b>IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
---
*Check if the order is Gorenstein at the prime ideal P.*


## List of instrinsics in AlgEtQ/TotRealTotPos.m:

> <pre><b>IsTotallyReal</b>(a::AlgEtQElt) -> BoolElt</pre>
---
*Returns whther a is totally real.*

> <pre><b>IsTotallyRealPositive</b>(a::AlgEtQElt) -> BoolElt</pre>
---
*Returns whether a is totally positive, that is, totally real and with positive image in C.*

> <pre><b>TotallyRealSubAlgebra</b>(K::AlgEtQ) -> AlgEtQ,Map</pre>
---
*Given a CM algebra K returns the unique totally real subalgebra, and an embedding.*

> <pre><b>TotallyRealUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre>
---
*Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^\*.*

> <pre><b>TotallyRealPositiveUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre>
---
*Given an order S in a CM étale algebra A. Returns the groups of totally positive units of S, as a subgroup of S^\*.*


## List of instrinsics in AlgEtQ/PrintSave.m:

> <pre><b>PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt</pre>
---
*Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.*

> <pre><b>PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt</pre>
---
*Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.*

> <pre><b>LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd</pre>
---
*Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics.*


## List of instrinsics in AlgEtQ/LowCohenMacaulayType.m:

> <pre><b>NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum</pre>
---
*Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.*

> <pre><b>CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt</pre>
---
*Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P\*S^t where S^t is the trace dual of S.*

> <pre><b>CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt</pre>
---
*Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P\*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.*


## List of instrinsics in AlgEtQ/WkClasses.m:

> <pre><b>WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum</pre>
---
*Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*

> <pre><b>WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum</pre>
---
*Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*


## List of instrinsics in AlgEtQ/WkTesting.m:

> <pre><b>IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt</pre>
---
*Checks if I and J are weakly equivalent 1 \in (I:J)\*(J:I). This function does not require that the ideals are defined over the same order.*

> <pre><b>IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
---
*Check if the two orders are weakly equivalent, that is equal.*

> <pre><b>IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt</pre>
---
*Checks if the second argument is weakly equivalent to the first argument.*

> <pre><b>IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt</pre>
---
*Checks if the second argument is weakly equivalent to the first argument.*

> <pre><b>IsGorenstein</b>(O::AlgEtQOrd)->BoolElt</pre>
---
*Checks if the order O is Gorenstein.*


## List of instrinsics in AlgEtQ/IdealClassMonoid.m:

> <pre><b>ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
---
*returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S).*

> <pre><b>ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
---
*returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals.*


