## List of instrinsics in AlgEtQ/AlgEt.m:

<p>
<b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ
</p>
*Given a sequence of number fields returns the étale algebra corresponding to the direct product.*

<p>
<b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ
</p>
*Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.*

<p>
<b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ
</p>
*Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.*


## List of instrinsics in AlgEtQ/AlgEtAttributes.m:

<p>
<b>Print</b>(A::AlgEtQ)
</p>
*Prints the defining polynomial or the components defining A.*

<p>
<b>DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt
</p>
*Returns the defining polynomial of A, if the corresponding number fields are distinct.*

<p>
<b>Components</b>(A::AlgEtQ) -> SeqEnum
</p>
*Returns the number fields of which A is a product of,together with embeddings and projections.*

<p>
<b>Dimension</b>(A::AlgEtQ)->RngInt
</p>
*Dimension of A.*

<p>
<b>AbsoluteDimension</b>(A::AlgEtQ)->RngInt
</p>
*Dimension of A over the prime field.*

<p>
<b>HasBaseField</b>(A::AlgEtQ) -> BoolElt,FldNum
</p>
*Returns whether A has common base field. If this is the case it returns it.*

<p>
<b>BaseField</b>(A::AlgEtQ) -> FldNum
</p>
*Returns the common base field of the Algebra, if it exists.*

<p>
<b>PrimeField</b>(A::AlgEtQ) -> FldNum
</p>
*Returns the prime field of the Algebra.*

<p>
<b>'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt
</p>
*A1 eq A2.*


## List of instrinsics in AlgEtQ/Homs.m:

<p>
<b>HomsToC</b>(A::AlgEtQ : Precision:=30)->SeqEnum[Map]
</p>
*returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30*


## List of instrinsics in AlgEtQ/Elt.m:

<p>
<b>Print</b>(x::AlgEtQElt)
</p>
*Print the AlgEtQElt.*

<p>
<b>Parent</b>(x::AlgEtQElt) -> AlgEtQ
</p>
*Returns the algebra to which the elemenet belongs to.*

<p>
<b>Algebra</b>(x::AlgEtQElt) -> AlgEtQ
</p>
*Returns the algebra to which the elemenet belongs to.*

<p>
<b>Components</b>(x::AlgEtQElt) -> SeqEnum
</p>
*Given an element x returns its components, which are elements of number fields.*

<p>
<b>AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum
</p>
*Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.*

<p>
<b>IsCoercible</b>(A::AlgEtQ, x::.) -> BoolElt, .
</p>
*Return whether x is coercible into A and the result of the coercion if so.*

<p>
<b>'!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt
</p>
*Coerce x into A.*

<p>
<b>One</b>(A::AlgEtQ) -> AlgEtQElt
</p>
*The multiplicative neutral element of A.*

<p>
<b>Zero</b>(A::AlgEtQ) -> AlgEtQElt
</p>
*The additive neutral element of A.*

<p>
<b>IsUnit</b>(x::AlgEtQElt) -> BoolElt
</p>
*Returns wheter x is a unit in A.*

<p>
<b>IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt
</p>
*Returns wheter x is a not unit in A.*

<p>
<b>IsZeroDivisor2</b>(x::AlgEtQElt) -> BoolElt
</p>
*Returns wheter x is a not unit in A.*

<p>
<b>Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt
</p>
*Random element of A. The Coefficients are bounded by the positive integer bd.*

<p>
<b>Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt
</p>
*Random element of A. The Coefficients are bounded by VarArg bd (default 3).*

<p>
<b>RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt
</p>
*Random unit of A. The Coefficients are bounded by the positive integer bd.*

<p>
<b>'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt
</p>
*Is x1=x2 ?*

<p>
<b>'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt
</p>
*Is x1=x2 ?*

<p>
<b>'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt
</p>
*Is x1=x2 ?*

<p>
<b>'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt
</p>
*Is x1=x2 ?*

<p>
<b>'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt
</p>
*Is x1=x2 ?*

<p>
<b>'+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
</p>
*x1+x2.*

<p>
<b>'-'</b>(x::AlgEtQElt) -> AlgEtQElt
</p>
*-x.*

<p>
<b>'-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
</p>
*x1-x2.*

<p>
<b>'*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>'*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
</p>
*x1*x2.*

<p>
<b>Inverse</b>(x::AlgEtQElt) -> AlgEtQElt
</p>
*1/x.*

<p>
<b>'^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt
</p>
*x^n.*

<p>
<b>'/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
</p>
*x1/x2.*

<p>
<b>'&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt
</p>
*Given a sequence of AlgEtQElt returns the sum of the entries.*

<p>
<b>'&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt
</p>
*Given a sequence of AlgEtQElt returns the product of the entries.*

<p>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
</p>
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<p>
<b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
</p>
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<p>
<b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
</p>
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<p>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt
</p>
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<p>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt
</p>
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<p>
<b>MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt
</p>
*Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.*

<p>
<b>MinimalPolynomial</b>(x::AlgEtQElt, F::Rng) -> RngUPolElt
</p>
*Returns the minimal polynommial over the ring F of the element x.*

<p>
<b>AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt
</p>
*Returns the minimal polynommial over the prime field of the element x.*

<p>
<b>IsIntegral</b>(x::AlgEtQElt) -> BoolElt
</p>
*Returns whether the element x is integral (over the integers).*

<p>
<b>Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt
</p>
*Evaluate the polynomial f at the element a.*

<p>
<b>PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt
</p>
*Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.*

<p>
<b>PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]
</p>
*Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.*

<p>
<b>Basis</b>(A::AlgEtQ) -> SeqEnum
</p>
*Returns a basis of the algebra over the common base field.*

<p>
<b>AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum
</p>
*Returns a basis of the algebra over the prime field.*

<p>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum
</p>
*Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.*

<p>
<b>OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum
</p>
*Returns the orthogonal ideampotent element of the étale algebra A.*

<p>
<b>Idempotents</b>(A::AlgEtQ) -> SeqEnum
</p>
*Returns the ideampotent element of the étale algebra A.*


## List of instrinsics in AlgEtQ/TraceNorm.m:

<p>
<b>Trace</b>(x::AlgEtQElt) -> Any
</p>
*Returns the trace of the element x of an étale algebra.*

<p>
<b>Norm</b>(x::AlgEtQElt) -> Any
</p>
*Returns the norm of the element x of an étale algebra.*

<p>
<b>AbsoluteTrace</b>(x::AlgEtQElt) -> Any
</p>
*Returns the absolute trace of the element x of an étale algebra.*

<p>
<b>AbsoluteNorm</b>(x::AlgEtQElt) -> Any
</p>
*Returns the absolute norm of the element x of an étale algebra.*

<p>
<b>TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Returns the trace dual ideal of an ideal in an order in an etale algebra.*

<p>
<b>TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl
</p>
*Returns the trace dual ideal of an order in an etale algebra.*


## List of instrinsics in AlgEtQ/Ord.m:

<p>
<b>Print</b>(A::AlgEtQOrd)
</p>
*Print the order.*

<p>
<b>IsCoercible</b>(S::AlgEtQOrd, x::.) -> BoolElt, Any
</p>
*Return whether x is coercible into S and the result if so.*

<p>
<b>Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd
</p>
*Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped. The vararg CheckIsKnownOrder determines if we check if the order is already known, i.e. in the attribute Algebra`KnownOrders. The default value is true.*

<p>
<b>Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd
</p>
*Given a sequence of order in the number fiedls defining the etale algebra A, generates the product order.*

<p>
<b>Algebra</b>(S::AlgEtQOrd) -> AlgEtQ
</p>
*Returns the algebra of the order.*

<p>
<b>myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]
</p>
*Hash function for AlgEtQOrd.*

<p>
<b>ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]
</p>
*Return a Z-basis of the order.*

<p>
<b>Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]
</p>
*Return a set of generators of the order.*

<p>
<b>'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
</p>
*Checks equality of orders in an etale Algebra.*

<p>
<b>'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt
</p>
*Inclusion of elements.*

<p>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum
</p>
*AbsoluteCoordinates with respect to the ZBasis.*

<p>
<b>'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt
</p>
*Inclusion of elements.*

<p>
<b>'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt
</p>
*Inclusion of elements.*

<p>
<b>One</b>(S::AlgEtQOrd)->AlgEtQElt
</p>
*Unit element of S.*

<p>
<b>Zero</b>(S::AlgEtQOrd)->AlgEtQElt
</p>
*Zero element of S.*

<p>
<b>Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt
</p>
*Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

<p>
<b>Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt
</p>
*Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

<p>
<b>IsKnownOrder</b>(~R::AlgEtQOrd)
</p>
*This procedure checks wheter the order R is already in the list of known orders of the algebra A of definition of R. If so then it replaces R with the copy stored in the attribute KnownOrders. If not it adds it to KnownOrders. This is done to avoid creating multiple copies of the same order.*

<p>
<b>EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd
</p>
*Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial.*

<p>
<b>ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd
</p>
*Given a product of number field A, returns the order consisting of the product of the equation orders of the number fields.*

<p>
<b>MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd
</p>
*Returns the maximal order of the étale algebra A.*

<p>
<b>IsMaximal</b>(S::AlgEtQOrd) -> BoolElt
</p>
*Returns wheter the given order is the maximal order of the étale algebra.*

<p>
<b>IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup
</p>
*Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.*

<p>
<b>Index</b>(T::AlgEtQOrd) -> FldRatElt
</p>
*Given an order T computes its index with respect to the basis of the algebra of T as a free Z-module.*

<p>
<b>Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any
</p>
*Given two orders T \subset S, returns [S:T] = #S/T.*

<p>
<b>'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt
</p>
*Checks if the first argument is inside the second.*

<p>
<b>'*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd
</p>
*Returns the order generated by the orders O1 and O2.*

<p>
<b>'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd
</p>
*Intersection of orders.*

<p>
<b>MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd
</p>
*Returns the multiplicator ring of an order R, that is R itself.*


## List of instrinsics in AlgEtQ/Quotients.m:

<p>
<b>Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map
</p>
*Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order.*

<p>
<b>Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map
</p>
*Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.*

<p>
<b>Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map
</p>
*Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.*

<p>
<b>ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map
</p>
*Given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S). We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.*

<p>
<b>ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map
</p>
*Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.*

<p>
<b>PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt
</p>
*Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^*.*

<p>
<b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map
</p>
*Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

<p>
<b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map
</p>
*Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

<p>
<b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map
</p>
*Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

<p>
<b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map
</p>
*Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*


## List of instrinsics in AlgEtQ/OverOrders.m:

<p>
<b>IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt
</p>
*Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.*

<p>
<b>MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]
</p>
*Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

<p>
<b>MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]
</p>
*Computes the minimal overorders of R.*

<p>
<b>OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
</p>
*Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

<p>
<b>OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
</p>
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*

<p>
<b>FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]
</p>
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*


## List of instrinsics in AlgEtQ/GraphOverOrders.m:

<p>
<b>GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir
</p>
*Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].*


## List of instrinsics in AlgEtQ/Idl.m:

<p>
<b>Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gens.*

<p>
<b>Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl
</p>
*Given an order S which is a product of orders S_i in the number fiedls generting the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal direct product of the I_i.*

<p>
<b>Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>'*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl
</p>
*Creates an ideal of S, generated by gen.*

<p>
<b>Print</b>(I::AlgEtQIdl)
</p>
*Prints the ideal.*

<p>
<b>'!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is in S, then IT=I.*

<p>
<b>Algebra</b>(I::AlgEtQIdl) -> AlgEtQ
</p>
*Returns the étale algebra in which the ideal lives.*

<p>
<b>Order</b>(I::AlgEtQIdl) -> AlgEtQOrd
</p>
*Returns the order of definition of the ideal.*

<p>
<b>ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]
</p>
*Returns a Z-basis of the ideal.*

<p>
<b>Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]
</p>
*Returns the generators of the ideal.*

<p>
<b>myHash</b>(I::AlgEtQIdl)->RngInt
</p>
*Hash function.*

<p>
<b>'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt
</p>
*Equality testing.*

<p>
<b>'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt
</p>
*Equality testing.*

<p>
<b>'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt
</p>
*Return if I eq S. I needs to be an ideal of S.*

<p>
<b>'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt
</p>
*Return if I eq S. I needs to be an ideal of S.*

<p>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum
</p>
*AbsoluteCoordiantes with respect to the ZBasis.*

<p>
<b>'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt
</p>
*Returns if x is in I.*

<p>
<b>'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt
</p>
*Returns if x is in I.*

<p>
<b>'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt
</p>
*Returns if x is in I.*

<p>
<b>'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt
</p>
*Given an ideal I of S, return if S subseteq I.*

<p>
<b>'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt
</p>
*Given an ideal I of S, return if I subseteq S.*

<p>
<b>'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt
</p>
*Checks if the first argument is inside the second. The ideals need to be fractional.*

<p>
<b>Index</b>(T::AlgEtQIdl) -> FldRatElt
</p>
*Given an ideal T computes its index with respect to the basis of the algebra of T as a free Z-module.*

<p>
<b>Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any
</p>
*Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I].*

<p>
<b>Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any
</p>
*Given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I].*

<p>
<b>OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl
</p>
*Given an order S returns the ideal 1*S which will be cached.*

<p>
<b>Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl
</p>
*Computes the conductor of an order, defined as he colon ideal (O:OK), where OK is the maximal order of the algebra.*

<p>
<b>'+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl
</p>
*Returns the sum of two ideals.*

<p>
<b>'*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl
</p>
*Product of two ideals.*

<p>
<b>'*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl
</p>
*Returns x*I.*

<p>
<b>'^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl
</p>
*nth power of an ideal.*

<p>
<b>'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl
</p>
*Given an ideal I of S, return S cap I.*

<p>
<b>'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl
</p>
*Given an ideal I of S, return S cap I.*

<p>
<b>'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl
</p>
*Given ideals I and J, return J cap I.*

<p>
<b>'&+'</b>(seq::SeqEnum[AlgEtQIdl])->AlgEtQIdl
</p>
*Returns the sum of the fractional ideals in the sequence.*

<p>
<b>ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl
</p>
*Computes the colon ideal (I:J) (as an O-ideal) of two O-idealsi.*

<p>
<b>ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl
</p>
*Computes the colon ideal (1*O:J) (as an O-ideal).*

<p>
<b>ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl
</p>
*Computes the colon ideal (I:1*O) (as an O-ideal).*

<p>
<b>IsInvertible</b>(I::AlgEtQIdl) ->BoolElt
</p>
*Checks if the ideal I is invertible in its order of definition O.*

<p>
<b>Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl
</p>
*Computes the inverse of an ideal of a maximal order.*

<p>
<b>MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd
</p>
*Given a fractional R-ideal I computes its multiplicator ring (I:I). If the overorders of R are known the corresponding overorder is returned, in order to preserve the known attributes.*

<p>
<b>IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup
</p>
*Return if the argument is a product of ideals in number fields, and if so return also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).*

<p>
<b>Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt
</p>
*Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

<p>
<b>Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt
</p>
*Returns a random (small coefficient) element of I. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

<p>
<b>IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt
</p>
*Given two integral ideals I and J of an order S, returns whether I+J=R.*

<p>
<b>IsIntegral</b>(I::AlgEtQIdl) -> BoolElt
</p>
*Returns wheter the ideal I of S is integral, that is I \subseteq S.*

<p>
<b>MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt
</p>
*Given a fractional S ideal I, returns the ideal d*I,d when d is the smallest integer such that d*I is integral in S.*

<p>
<b>MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt
</p>
*Returns the smallest integer contained in the ideal I.*

<p>
<b>CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl
</p>
*Returns an element x such that x*I is an integral ideal coprime with J, togheter with the product x*I. The first ideal must be invertible and the second should be integral.*


## List of instrinsics in AlgEtQ/ZBasisLLL.m:

<p>
<b>ZBasisLLL</b>(S::AlgEtQOrd)
</p>
*A procedure that replaces the ZBasis with an LLL-reduced one.*

<p>
<b>ZBasisLLL</b>(S::AlgEtQIdl)
</p>
*A procedure that replaces the ZBasis with an LLL-reduced one.*


## List of instrinsics in AlgEtQ/Completion.m:

<p>
<b>Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map
</p>
*Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage (it acts a bit weird).*


## List of instrinsics in AlgEtQ/ComplexConj.m:

<p>
<b>HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt
</p>
*Returns if the algebra is the product of CM fields.*

<p>
<b>ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt
</p>
*If A is a product of CM fields, it returns the complex conjugate of the argument.*

<p>
<b>IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd
</p>
*Returns wheter O is conjugate stable and the complex conjugate.*

<p>
<b>ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd
</p>
*It returns the complex conjugate of the argument.*

<p>
<b>IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl
</p>
*Returns wheter O is conjugate stable and the complex conjugate.*

<p>
<b>ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl
</p>
*If A is a product of CM fields, it returns the complex conjugate of the argument.*


## List of instrinsics in AlgEtQ/ComplexMult.m:

<p>
<b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType
</p>
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

<p>
<b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType
</p>
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

<p>
<b>CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType
</p>
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

<p>
<b>CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType
</p>
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

<p>
<b>Print</b>( PHI :: AlgEtQCMType)
</p>
*Print the AlgEtQCMType.*

<p>
<b>CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt
</p>
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

<p>
<b>CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt
</p>
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

<p>
<b>Homs</b>( PHI::AlgEtQCMType : prec:=30 )->SeqEnum[Map]
</p>
*Given a AlgEtQCMType PHI returns the sequence of maps to CC defining it.*

<p>
<b>'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=30)->BoolElt
</p>
*Returns whether two cm types are equal. This happens if and only if the ration of (any) two CMPositiveElements is totally real and totally positive.*

<p>
<b>Precision</b>(PHI :: AlgEtQCMType)->RngIntElt
</p>
*Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

<p>
<b>ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType
</p>
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

<p>
<b>ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )
</p>
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

<p>
<b>AllCMTypes</b>(A::AlgEtQ : Precision := 30 ) -> SeqEnum[AlgEtQCMType]
</p>
*Returns all the AlgEtQCMTypes of A.*


## List of instrinsics in AlgEtQ/IntermediateIdeals.m:

<p>
<b>MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<p>
<b>IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.*

<p>
<b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. They are produced recursively using from the minimal ones.*

<p>
<b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<p>
<b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<p>
<b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

<p>
<b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

<p>
<b>IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]
</p>
*Given ideals J subset I over the same order, and a positive integer N, it returns all the ideals K such that J subset K subset I and [I:K]=N. These are produced by recursively searching for maximal submodules.*

<p>
<b>MinimalIntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl : primes:=[])->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I.*

<p>
<b>IntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
</p>
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones*


## List of instrinsics in AlgEtQ/IdealsOfIndex.m:

<p>
<b>IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]
</p>
*Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.*

<p>
<b>IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]
</p>
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

<p>
<b>IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]
</p>
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

<p>
<b>IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]
</p>
*Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".*

<p>
<b>IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]
</p>
*Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".*


## List of instrinsics in AlgEtQ/ShortEltSmallRep.m:

<p>
<b>ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt
</p>
*Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).*

<p>
<b>SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt
</p>
*Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.*


## List of instrinsics in AlgEtQ/MinimalGenerators.m:

<p>
<b>TwoGeneratingSet</b>(I::AlgEtQIdl)
</p>
*A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.*


## List of instrinsics in AlgEtQ/CRT.m:

<p>
<b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt
</p>
*Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.*

<p>
<b>ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt
</p>
*Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.*


## List of instrinsics in AlgEtQ/PicardGroup.m:

<p>
<b>ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map
</p>
*Returns the group (S/I)^* and a map (S/I)^* -> S. It is required S to be maximal.*

<p>
<b>IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgAssElt
</p>
*Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*

<p>
<b>PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
</p>
*Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".*

<p>
<b>UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
</p>
*Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".*

<p>
<b>IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgAssElt
</p>
*Checks if I=x*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*


## List of instrinsics in AlgEtQ/FactPrimes.m:

<p>
<b>Factorization</b>(I::AlgEtQIdl) -> Tup
</p>
*Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.*

<p>
<b>PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgAssEtOrdIdl]
</p>
*Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.*

<p>
<b>SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgAssEtOrdIdl]
</p>
*Returns the non-invertible primes of the order.*

<p>
<b>NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx
</p>
*Returns the non-invertible primes of the order.*

<p>
<b>IsPrime</b>(I::AlgEtQIdl) -> BoolElt
</p>
*Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal.*

<p>
<b>IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
</p>
*Check if the order is Bass at the prime ideal P.*

<p>
<b>IsBass</b>(S::AlgEtQOrd) -> BoolElt
</p>
*Check if the order is Bass.*

<p>
<b>IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
</p>
*Check if the order is Gorenstein at the prime ideal P.*


## List of instrinsics in AlgEtQ/TotRealTotPos.m:

<p>
<b>IsTotallyReal</b>(a::AlgEtQElt) -> BoolElt
</p>
*Returns whther a is totally real.*

<p>
<b>IsTotallyRealPositive</b>(a::AlgEtQElt) -> BoolElt
</p>
*Returns whether a is totally positive, that is, totally real and with positive image in C.*

<p>
<b>TotallyRealSubAlgebra</b>(K::AlgEtQ) -> AlgEtQ,Map
</p>
*Given a CM algebra K returns the unique totally real subalgebra, and an embedding.*

<p>
<b>TotallyRealUnitGroup</b>(S::AlgEtQOrd) -> Grp
</p>
*Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^*.*

<p>
<b>TotallyRealPositiveUnitGroup</b>(S::AlgEtQOrd) -> Grp
</p>
*Given an order S in a CM étale algebra A. Returns the groups of totally positive units of S, as a subgroup of S^*.*


## List of instrinsics in AlgEtQ/PrintSave.m:

<p>
<b>PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt
</p>
*Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.*

<p>
<b>PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt
</p>
*Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.*

<p>
<b>LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd
</p>
*Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics.*


## List of instrinsics in AlgEtQ/LowCohenMacaulayType.m:

<p>
<b>NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum
</p>
*Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.*

<p>
<b>CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt
</p>
*Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P*S^t where S^t is the trace dual of S.*

<p>
<b>CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt
</p>
*Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.*


## List of instrinsics in AlgEtQ/WkClasses.m:

<p>
<b>WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum
</p>
*Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*

<p>
<b>WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum
</p>
*Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*


## List of instrinsics in AlgEtQ/WkTesting.m:

<p>
<b>IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt
</p>
*Checks if I and J are weakly equivalent 1 \in (I:J)*(J:I). This function does not require that the ideals are defined over the same order.*

<p>
<b>IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
</p>
*Check if the two orders are weakly equivalent, that is equal.*

<p>
<b>IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt
</p>
*Checks if the second argument is weakly equivalent to the first argument.*

<p>
<b>IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt
</p>
*Checks if the second argument is weakly equivalent to the first argument.*

<p>
<b>IsGorenstein</b>(O::AlgEtQOrd)->BoolElt
</p>
*Checks if the order O is Gorenstein.*


## List of instrinsics in AlgEtQ/IdealClassMonoid.m:

<p>
<b>ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
</p>
*returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S).*

<p>
<b>ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
</p>
*returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals.*


