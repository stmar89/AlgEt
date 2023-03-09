## List of instrinsics in AlgEtQ/AlgEt.m:

> <dt><pre><b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ</pre></dt>
<dd>Given a sequence of number fields returns the étale algebra corresponding to the direct product.</dd>

> <dt><pre><b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ</pre></dt>
<dd>Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.</dd>

> <dt><pre><b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ</pre></dt>
<dd>Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.</dd>


## List of instrinsics in AlgEtQ/AlgEtAttributes.m:

> <dt><pre><b>Print</b>(A::AlgEtQ)</pre></dt>
<dd>Prints the defining polynomial or the components defining A.</dd>

> <dt><pre><b>DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt</pre></dt>
<dd>Returns the defining polynomial of A, if the corresponding number fields are distinct.</dd>

> <dt><pre><b>Components</b>(A::AlgEtQ) -> SeqEnum</pre></dt>
<dd>Returns the number fields of which A is a product of,together with embeddings and projections.</dd>

> <dt><pre><b>Dimension</b>(A::AlgEtQ)->RngInt</pre></dt>
<dd>Dimension of A.</dd>

> <dt><pre><b>AbsoluteDimension</b>(A::AlgEtQ)->RngInt</pre></dt>
<dd>Dimension of A over the prime field.</dd>

> <dt><pre><b>HasBaseField</b>(A::AlgEtQ) -> BoolElt,FldNum</pre></dt>
<dd>Returns whether A has common base field. If this is the case it returns it.</dd>

> <dt><pre><b>BaseField</b>(A::AlgEtQ) -> FldNum</pre></dt>
<dd>Returns the common base field of the Algebra, if it exists.</dd>

> <dt><pre><b>PrimeField</b>(A::AlgEtQ) -> FldNum</pre></dt>
<dd>Returns the prime field of the Algebra.</dd>

> <dt><pre><b>'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt</pre></dt>
<dd>A1 eq A2.</dd>


## List of instrinsics in AlgEtQ/Homs.m:

> <dt><pre><b>HomsToC</b>(A::AlgEtQ : Precision:=30)->SeqEnum[Map]</pre></dt>
<dd>returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30</dd>


## List of instrinsics in AlgEtQ/Elt.m:

> <dt><pre><b>Print</b>(x::AlgEtQElt)</pre></dt>
<dd>Print the AlgEtQElt.</dd>

> <dt><pre><b>Parent</b>(x::AlgEtQElt) -> AlgEtQ</pre></dt>
<dd>Returns the algebra to which the elemenet belongs to.</dd>

> <dt><pre><b>Algebra</b>(x::AlgEtQElt) -> AlgEtQ</pre></dt>
<dd>Returns the algebra to which the elemenet belongs to.</dd>

> <dt><pre><b>Components</b>(x::AlgEtQElt) -> SeqEnum</pre></dt>
<dd>Given an element x returns its components, which are elements of number fields.</dd>

> <dt><pre><b>AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum</pre></dt>
<dd>Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.</dd>

> <dt><pre><b>IsCoercible</b>(A::AlgEtQ, x::.) -> BoolElt, .</pre></dt>
<dd>Return whether x is coercible into A and the result of the coercion if so.</dd>

> <dt><pre><b>'!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt</pre></dt>
<dd>Coerce x into A.</dd>

> <dt><pre><b>One</b>(A::AlgEtQ) -> AlgEtQElt</pre></dt>
<dd>The multiplicative neutral element of A.</dd>

> <dt><pre><b>Zero</b>(A::AlgEtQ) -> AlgEtQElt</pre></dt>
<dd>The additive neutral element of A.</dd>

> <dt><pre><b>IsUnit</b>(x::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns wheter x is a unit in A.</dd>

> <dt><pre><b>IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns wheter x is a not unit in A.</dd>

> <dt><pre><b>IsZeroDivisor2</b>(x::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns wheter x is a not unit in A.</dd>

> <dt><pre><b>Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>Random element of A. The Coefficients are bounded by the positive integer bd.</dd>

> <dt><pre><b>Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt</pre></dt>
<dd>Random element of A. The Coefficients are bounded by VarArg bd (default 3).</dd>

> <dt><pre><b>RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>Random unit of A. The Coefficients are bounded by the positive integer bd.</dd>

> <dt><pre><b>'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Is x1=x2 ?</dd>

> <dt><pre><b>'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Is x1=x2 ?</dd>

> <dt><pre><b>'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Is x1=x2 ?</dd>

> <dt><pre><b>'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt</pre></dt>
<dd>Is x1=x2 ?</dd>

> <dt><pre><b>'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt</pre></dt>
<dd>Is x1=x2 ?</dd>

> <dt><pre><b>'+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre></dt>
<dd>x1+x2.</dd>

> <dt><pre><b>'-'</b>(x::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>-x.</dd>

> <dt><pre><b>'-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre></dt>
<dd>x1-x2.</dd>

> <dt><pre><b>'*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>'*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre></dt>
<dd>x1\*x2.</dd>

> <dt><pre><b>Inverse</b>(x::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>1/x.</dd>

> <dt><pre><b>'^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>x^n.</dd>

> <dt><pre><b>'/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre></dt>
<dd>x1/x2.</dd>

> <dt><pre><b>'&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre></dt>
<dd>Given a sequence of AlgEtQElt returns the sum of the entries.</dd>

> <dt><pre><b>'&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre></dt>
<dd>Given a sequence of AlgEtQElt returns the product of the entries.</dd>

> <dt><pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre></dt>
<dd>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</dd>

> <dt><pre><b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre></dt>
<dd>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</dd>

> <dt><pre><b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre></dt>
<dd>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</dd>

> <dt><pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt</pre></dt>
<dd>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</dd>

> <dt><pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt</pre></dt>
<dd>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</dd>

> <dt><pre><b>MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre></dt>
<dd>Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.</dd>

> <dt><pre><b>MinimalPolynomial</b>(x::AlgEtQElt, F::Rng) -> RngUPolElt</pre></dt>
<dd>Returns the minimal polynommial over the ring F of the element x.</dd>

> <dt><pre><b>AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre></dt>
<dd>Returns the minimal polynommial over the prime field of the element x.</dd>

> <dt><pre><b>IsIntegral</b>(x::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns whether the element x is integral (over the integers).</dd>

> <dt><pre><b>Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>Evaluate the polynomial f at the element a.</dd>

> <dt><pre><b>PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt</pre></dt>
<dd>Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.</dd>

> <dt><pre><b>PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]</pre></dt>
<dd>Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.</dd>

> <dt><pre><b>Basis</b>(A::AlgEtQ) -> SeqEnum</pre></dt>
<dd>Returns a basis of the algebra over the common base field.</dd>

> <dt><pre><b>AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum</pre></dt>
<dd>Returns a basis of the algebra over the prime field.</dd>

> <dt><pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum</pre></dt>
<dd>Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.</dd>

> <dt><pre><b>OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum</pre></dt>
<dd>Returns the orthogonal ideampotent element of the étale algebra A.</dd>

> <dt><pre><b>Idempotents</b>(A::AlgEtQ) -> SeqEnum</pre></dt>
<dd>Returns the ideampotent element of the étale algebra A.</dd>


## List of instrinsics in AlgEtQ/TraceNorm.m:

> <dt><pre><b>Trace</b>(x::AlgEtQElt) -> Any</pre></dt>
<dd>Returns the trace of the element x of an étale algebra.</dd>

> <dt><pre><b>Norm</b>(x::AlgEtQElt) -> Any</pre></dt>
<dd>Returns the norm of the element x of an étale algebra.</dd>

> <dt><pre><b>AbsoluteTrace</b>(x::AlgEtQElt) -> Any</pre></dt>
<dd>Returns the absolute trace of the element x of an étale algebra.</dd>

> <dt><pre><b>AbsoluteNorm</b>(x::AlgEtQElt) -> Any</pre></dt>
<dd>Returns the absolute norm of the element x of an étale algebra.</dd>

> <dt><pre><b>TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Returns the trace dual ideal of an ideal in an order in an etale algebra.</dd>

> <dt><pre><b>TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Returns the trace dual ideal of an order in an etale algebra.</dd>


## List of instrinsics in AlgEtQ/Ord.m:

> <dt><pre><b>Print</b>(A::AlgEtQOrd)</pre></dt>
<dd>Print the order.</dd>

> <dt><pre><b>IsCoercible</b>(S::AlgEtQOrd, x::.) -> BoolElt, Any</pre></dt>
<dd>Return whether x is coercible into S and the result if so.</dd>

> <dt><pre><b>Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd</pre></dt>
<dd>Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped. The vararg CheckIsKnownOrder determines if we check if the order is already known, i.e. in the attribute Algebra`KnownOrders. The default value is true.</dd>

> <dt><pre><b>Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd</pre></dt>
<dd>Given a sequence of order in the number fiedls defining the etale algebra A, generates the product order.</dd>

> <dt><pre><b>Algebra</b>(S::AlgEtQOrd) -> AlgEtQ</pre></dt>
<dd>Returns the algebra of the order.</dd>

> <dt><pre><b>myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]</pre></dt>
<dd>Hash function for AlgEtQOrd.</dd>

> <dt><pre><b>ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre></dt>
<dd>Return a Z-basis of the order.</dd>

> <dt><pre><b>Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre></dt>
<dd>Return a set of generators of the order.</dd>

> <dt><pre><b>'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre></dt>
<dd>Checks equality of orders in an etale Algebra.</dd>

> <dt><pre><b>'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Inclusion of elements.</dd>

> <dt><pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum</pre></dt>
<dd>AbsoluteCoordinates with respect to the ZBasis.</dd>

> <dt><pre><b>'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Inclusion of elements.</dd>

> <dt><pre><b>'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Inclusion of elements.</dd>

> <dt><pre><b>One</b>(S::AlgEtQOrd)->AlgEtQElt</pre></dt>
<dd>Unit element of S.</dd>

> <dt><pre><b>Zero</b>(S::AlgEtQOrd)->AlgEtQElt</pre></dt>
<dd>Zero element of S.</dd>

> <dt><pre><b>Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre></dt>
<dd>Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</dd>

> <dt><pre><b>Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre></dt>
<dd>Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</dd>

> <dt><pre><b>IsKnownOrder</b>(~R::AlgEtQOrd)</pre></dt>
<dd>This procedure checks wheter the order R is already in the list of known orders of the algebra A of definition of R. If so then it replaces R with the copy stored in the attribute KnownOrders. If not it adds it to KnownOrders. This is done to avoid creating multiple copies of the same order.</dd>

> <dt><pre><b>EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd</pre></dt>
<dd>Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial.</dd>

> <dt><pre><b>ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd</pre></dt>
<dd>Given a product of number field A, returns the order consisting of the product of the equation orders of the number fields.</dd>

> <dt><pre><b>MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd</pre></dt>
<dd>Returns the maximal order of the étale algebra A.</dd>

> <dt><pre><b>IsMaximal</b>(S::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Returns wheter the given order is the maximal order of the étale algebra.</dd>

> <dt><pre><b>IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup</pre></dt>
<dd>Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.</dd>

> <dt><pre><b>Index</b>(T::AlgEtQOrd) -> FldRatElt</pre></dt>
<dd>Given an order T computes its index with respect to the basis of the algebra of T as a free Z-module.</dd>

> <dt><pre><b>Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any</pre></dt>
<dd>Given two orders T \subset S, returns [S:T] = #S/T.</dd>

> <dt><pre><b>'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Checks if the first argument is inside the second.</dd>

> <dt><pre><b>'*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre></dt>
<dd>Returns the order generated by the orders O1 and O2.</dd>

> <dt><pre><b>'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre></dt>
<dd>Intersection of orders.</dd>

> <dt><pre><b>MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd</pre></dt>
<dd>Returns the multiplicator ring of an order R, that is R itself.</dd>


## List of instrinsics in AlgEtQ/Quotients.m:

> <dt><pre><b>Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre></dt>
<dd>Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order.</dd>

> <dt><pre><b>Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map</pre></dt>
<dd>Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.</dd>

> <dt><pre><b>Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre></dt>
<dd>Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.</dd>

> <dt><pre><b>ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map</pre></dt>
<dd>Given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S). We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.</dd>

> <dt><pre><b>ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map</pre></dt>
<dd>Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.</dd>

> <dt><pre><b>PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt</pre></dt>
<dd>Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^\*.</dd>

> <dt><pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre></dt>
<dd>Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).</dd>

> <dt><pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre></dt>
<dd>Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).</dd>

> <dt><pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre></dt>
<dd>Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).</dd>

> <dt><pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre></dt>
<dd>Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image).</dd>


## List of instrinsics in AlgEtQ/OverOrders.m:

> <dt><pre><b>IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.</dd>

> <dt><pre><b>MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]</pre></dt>
<dd>Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.</dd>

> <dt><pre><b>MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]</pre></dt>
<dd>Computes the minimal overorders of R.</dd>

> <dt><pre><b>OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]</pre></dt>
<dd>Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.</dd>

> <dt><pre><b>OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]</pre></dt>
<dd>We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.</dd>

> <dt><pre><b>FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]</pre></dt>
<dd>We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.</dd>


## List of instrinsics in AlgEtQ/GraphOverOrders.m:

> <dt><pre><b>GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir</pre></dt>
<dd>Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].</dd>


## List of instrinsics in AlgEtQ/Idl.m:

> <dt><pre><b>Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gens.</dd>

> <dt><pre><b>Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl</pre></dt>
<dd>Given an order S which is a product of orders S_i in the number fiedls generting the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal direct product of the I_i.</dd>

> <dt><pre><b>Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>'*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Creates an ideal of S, generated by gen.</dd>

> <dt><pre><b>Print</b>(I::AlgEtQIdl)</pre></dt>
<dd>Prints the ideal.</dd>

> <dt><pre><b>'!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is in S, then IT=I.</dd>

> <dt><pre><b>Algebra</b>(I::AlgEtQIdl) -> AlgEtQ</pre></dt>
<dd>Returns the étale algebra in which the ideal lives.</dd>

> <dt><pre><b>Order</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre></dt>
<dd>Returns the order of definition of the ideal.</dd>

> <dt><pre><b>ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]</pre></dt>
<dd>Returns a Z-basis of the ideal.</dd>

> <dt><pre><b>Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre></dt>
<dd>Returns the generators of the ideal.</dd>

> <dt><pre><b>myHash</b>(I::AlgEtQIdl)->RngInt</pre></dt>
<dd>Hash function.</dd>

> <dt><pre><b>'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre></dt>
<dd>Equality testing.</dd>

> <dt><pre><b>'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre></dt>
<dd>Equality testing.</dd>

> <dt><pre><b>'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Return if I eq S. I needs to be an ideal of S.</dd>

> <dt><pre><b>'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Return if I eq S. I needs to be an ideal of S.</dd>

> <dt><pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum</pre></dt>
<dd>AbsoluteCoordiantes with respect to the ZBasis.</dd>

> <dt><pre><b>'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt</pre></dt>
<dd>Returns if x is in I.</dd>

> <dt><pre><b>'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt</pre></dt>
<dd>Returns if x is in I.</dd>

> <dt><pre><b>'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt</pre></dt>
<dd>Returns if x is in I.</dd>

> <dt><pre><b>'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Given an ideal I of S, return if S subseteq I.</dd>

> <dt><pre><b>'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Given an ideal I of S, return if I subseteq S.</dd>

> <dt><pre><b>'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Checks if the first argument is inside the second. The ideals need to be fractional.</dd>

> <dt><pre><b>Index</b>(T::AlgEtQIdl) -> FldRatElt</pre></dt>
<dd>Given an ideal T computes its index with respect to the basis of the algebra of T as a free Z-module.</dd>

> <dt><pre><b>Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any</pre></dt>
<dd>Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I].</dd>

> <dt><pre><b>Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any</pre></dt>
<dd>Given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I].</dd>

> <dt><pre><b>OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Given an order S returns the ideal 1\*S which will be cached.</dd>

> <dt><pre><b>Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl</pre></dt>
<dd>Computes the conductor of an order, defined as he colon ideal (O:OK), where OK is the maximal order of the algebra.</dd>

> <dt><pre><b>'+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre></dt>
<dd>Returns the sum of two ideals.</dd>

> <dt><pre><b>'*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre></dt>
<dd>Product of two ideals.</dd>

> <dt><pre><b>'*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl</pre></dt>
<dd>Returns x\*I.</dd>

> <dt><pre><b>'^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl</pre></dt>
<dd>nth power of an ideal.</dd>

> <dt><pre><b>'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl</pre></dt>
<dd>Given an ideal I of S, return S cap I.</dd>

> <dt><pre><b>'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Given an ideal I of S, return S cap I.</dd>

> <dt><pre><b>'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>Given ideals I and J, return J cap I.</dd>

> <dt><pre><b>'&+'</b>(seq::SeqEnum[AlgEtQIdl])->AlgEtQIdl</pre></dt>
<dd>Returns the sum of the fractional ideals in the sequence.</dd>

> <dt><pre><b>ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl</pre></dt>
<dd>Computes the colon ideal (I:J) (as an O-ideal) of two O-idealsi.</dd>

> <dt><pre><b>ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl</pre></dt>
<dd>Computes the colon ideal (1\*O:J) (as an O-ideal).</dd>

> <dt><pre><b>ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl</pre></dt>
<dd>Computes the colon ideal (I:1\*O) (as an O-ideal).</dd>

> <dt><pre><b>IsInvertible</b>(I::AlgEtQIdl) ->BoolElt</pre></dt>
<dd>Checks if the ideal I is invertible in its order of definition O.</dd>

> <dt><pre><b>Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl</pre></dt>
<dd>Computes the inverse of an ideal of a maximal order.</dd>

> <dt><pre><b>MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre></dt>
<dd>Given a fractional R-ideal I computes its multiplicator ring (I:I). If the overorders of R are known the corresponding overorder is returned, in order to preserve the known attributes.</dd>

> <dt><pre><b>IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup</pre></dt>
<dd>Return if the argument is a product of ideals in number fields, and if so return also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).</dd>

> <dt><pre><b>Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre></dt>
<dd>Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</dd>

> <dt><pre><b>Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre></dt>
<dd>Returns a random (small coefficient) element of I. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</dd>

> <dt><pre><b>IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Given two integral ideals I and J of an order S, returns whether I+J=R.</dd>

> <dt><pre><b>IsIntegral</b>(I::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Returns wheter the ideal I of S is integral, that is I \subseteq S.</dd>

> <dt><pre><b>MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt</pre></dt>
<dd>Given a fractional S ideal I, returns the ideal d\*I,d when d is the smallest integer such that d\*I is integral in S.</dd>

> <dt><pre><b>MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt</pre></dt>
<dd>Returns the smallest integer contained in the ideal I.</dd>

> <dt><pre><b>CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl</pre></dt>
<dd>Returns an element x such that x\*I is an integral ideal coprime with J, togheter with the product x\*I. The first ideal must be invertible and the second should be integral.</dd>


## List of instrinsics in AlgEtQ/ZBasisLLL.m:

> <dt><pre><b>ZBasisLLL</b>(S::AlgEtQOrd)</pre></dt>
<dd>A procedure that replaces the ZBasis with an LLL-reduced one.</dd>

> <dt><pre><b>ZBasisLLL</b>(S::AlgEtQIdl)</pre></dt>
<dd>A procedure that replaces the ZBasis with an LLL-reduced one.</dd>


## List of instrinsics in AlgEtQ/Completion.m:

> <dt><pre><b>Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map</pre></dt>
<dd>Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage (it acts a bit weird).</dd>


## List of instrinsics in AlgEtQ/ComplexConj.m:

> <dt><pre><b>HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt</pre></dt>
<dd>Returns if the algebra is the product of CM fields.</dd>

> <dt><pre><b>ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt</pre></dt>
<dd>If A is a product of CM fields, it returns the complex conjugate of the argument.</dd>

> <dt><pre><b>IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd</pre></dt>
<dd>Returns wheter O is conjugate stable and the complex conjugate.</dd>

> <dt><pre><b>ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd</pre></dt>
<dd>It returns the complex conjugate of the argument.</dd>

> <dt><pre><b>IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl</pre></dt>
<dd>Returns wheter O is conjugate stable and the complex conjugate.</dd>

> <dt><pre><b>ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre></dt>
<dd>If A is a product of CM fields, it returns the complex conjugate of the argument.</dd>


## List of instrinsics in AlgEtQ/ComplexMult.m:

> <dt><pre><b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre></dt>
<dd>Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.</dd>

> <dt><pre><b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre></dt>
<dd>Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.</dd>

> <dt><pre><b>CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre></dt>
<dd>Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.</dd>

> <dt><pre><b>CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre></dt>
<dd>Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.</dd>

> <dt><pre><b>Print</b>( PHI :: AlgEtQCMType)</pre></dt>
<dd>Print the AlgEtQCMType.</dd>

> <dt><pre><b>CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre></dt>
<dd>Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).</dd>

> <dt><pre><b>CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre></dt>
<dd>Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).</dd>

> <dt><pre><b>Homs</b>( PHI::AlgEtQCMType : prec:=30 )->SeqEnum[Map]</pre></dt>
<dd>Given a AlgEtQCMType PHI returns the sequence of maps to CC defining it.</dd>

> <dt><pre><b>'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=30)->BoolElt</pre></dt>
<dd>Returns whether two cm types are equal. This happens if and only if the ration of (any) two CMPositiveElements is totally real and totally positive.</dd>

> <dt><pre><b>Precision</b>(PHI :: AlgEtQCMType)->RngIntElt</pre></dt>
<dd>Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</dd>

> <dt><pre><b>ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType</pre></dt>
<dd>Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</dd>

> <dt><pre><b>ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )</pre></dt>
<dd>Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</dd>

> <dt><pre><b>AllCMTypes</b>(A::AlgEtQ : Precision := 30 ) -> SeqEnum[AlgEtQCMType]</pre></dt>
<dd>Returns all the AlgEtQCMTypes of A.</dd>


## List of instrinsics in AlgEtQ/IntermediateIdeals.m:

> <dt><pre><b>MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.</dd>

> <dt><pre><b>IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.</dd>

> <dt><pre><b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. They are produced recursively using from the minimal ones.</dd>

> <dt><pre><b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.</dd>

> <dt><pre><b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.</dd>

> <dt><pre><b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.</dd>

> <dt><pre><b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.</dd>

> <dt><pre><b>IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given ideals J subset I over the same order, and a positive integer N, it returns all the ideals K such that J subset K subset I and [I:K]=N. These are produced by recursively searching for maximal submodules.</dd>

> <dt><pre><b>MinimalIntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl : primes:=[])->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I.</dd>

> <dt><pre><b>IntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre></dt>
<dd>Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones</dd>


## List of instrinsics in AlgEtQ/IdealsOfIndex.m:

> <dt><pre><b>IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre></dt>
<dd>Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.</dd>

> <dt><pre><b>IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre></dt>
<dd>Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.</dd>

> <dt><pre><b>IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]</pre></dt>
<dd>Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.</dd>

> <dt><pre><b>IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre></dt>
<dd>Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".</dd>

> <dt><pre><b>IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre></dt>
<dd>Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".</dd>


## List of instrinsics in AlgEtQ/ShortEltSmallRep.m:

> <dt><pre><b>ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt</pre></dt>
<dd>Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).</dd>

> <dt><pre><b>SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt</pre></dt>
<dd>Given a fractional R-ideal I, it returns an isomorphic ideal a\*I, and the element a, such that a\*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.</dd>


## List of instrinsics in AlgEtQ/MinimalGenerators.m:

> <dt><pre><b>TwoGeneratingSet</b>(I::AlgEtQIdl)</pre></dt>
<dd>A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.</dd>


## List of instrinsics in AlgEtQ/CRT.m:

> <dt><pre><b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt</pre></dt>
<dd>Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.</dd>

> <dt><pre><b>ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt</pre></dt>
<dd>Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.</dd>


## List of instrinsics in AlgEtQ/PicardGroup.m:

> <dt><pre><b>ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map</pre></dt>
<dd>Returns the group (S/I)^\* and a map (S/I)^\* -> S. It is required S to be maximal.</dd>

> <dt><pre><b>IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgAssElt</pre></dt>
<dd>Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".</dd>

> <dt><pre><b>PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre></dt>
<dd>Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".</dd>

> <dt><pre><b>UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre></dt>
<dd>Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".</dd>

> <dt><pre><b>IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgAssElt</pre></dt>
<dd>Checks if I=x\*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".</dd>


## List of instrinsics in AlgEtQ/FactPrimes.m:

> <dt><pre><b>Factorization</b>(I::AlgEtQIdl) -> Tup</pre></dt>
<dd>Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.</dd>

> <dt><pre><b>PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgAssEtOrdIdl]</pre></dt>
<dd>Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.</dd>

> <dt><pre><b>SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgAssEtOrdIdl]</pre></dt>
<dd>Returns the non-invertible primes of the order.</dd>

> <dt><pre><b>NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx</pre></dt>
<dd>Returns the non-invertible primes of the order.</dd>

> <dt><pre><b>IsPrime</b>(I::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal.</dd>

> <dt><pre><b>IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Check if the order is Bass at the prime ideal P.</dd>

> <dt><pre><b>IsBass</b>(S::AlgEtQOrd) -> BoolElt</pre></dt>
<dd>Check if the order is Bass.</dd>

> <dt><pre><b>IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre></dt>
<dd>Check if the order is Gorenstein at the prime ideal P.</dd>


## List of instrinsics in AlgEtQ/TotRealTotPos.m:

> <dt><pre><b>IsTotallyReal</b>(a::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns whther a is totally real.</dd>

> <dt><pre><b>IsTotallyRealPositive</b>(a::AlgEtQElt) -> BoolElt</pre></dt>
<dd>Returns whether a is totally positive, that is, totally real and with positive image in C.</dd>

> <dt><pre><b>TotallyRealSubAlgebra</b>(K::AlgEtQ) -> AlgEtQ,Map</pre></dt>
<dd>Given a CM algebra K returns the unique totally real subalgebra, and an embedding.</dd>

> <dt><pre><b>TotallyRealUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre></dt>
<dd>Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^\*.</dd>

> <dt><pre><b>TotallyRealPositiveUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre></dt>
<dd>Given an order S in a CM étale algebra A. Returns the groups of totally positive units of S, as a subgroup of S^\*.</dd>


## List of instrinsics in AlgEtQ/PrintSave.m:

> <dt><pre><b>PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt</pre></dt>
<dd>Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.</dd>

> <dt><pre><b>PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt</pre></dt>
<dd>Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.</dd>

> <dt><pre><b>LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd</pre></dt>
<dd>Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics.</dd>


## List of instrinsics in AlgEtQ/LowCohenMacaulayType.m:

> <dt><pre><b>NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum</pre></dt>
<dd>Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.</dd>

> <dt><pre><b>CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt</pre></dt>
<dd>Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P\*S^t where S^t is the trace dual of S.</dd>

> <dt><pre><b>CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt</pre></dt>
<dd>Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P\*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.</dd>


## List of instrinsics in AlgEtQ/WkClasses.m:

> <dt><pre><b>WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum</pre></dt>
<dd>Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.</dd>

> <dt><pre><b>WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum</pre></dt>
<dd>Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.</dd>


## List of instrinsics in AlgEtQ/WkTesting.m:

> <dt><pre><b>IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt</pre></dt>
<dd>Checks if I and J are weakly equivalent 1 \in (I:J)\*(J:I). This function does not require that the ideals are defined over the same order.</dd>

> <dt><pre><b>IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre></dt>
<dd>Check if the two orders are weakly equivalent, that is equal.</dd>

> <dt><pre><b>IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt</pre></dt>
<dd>Checks if the second argument is weakly equivalent to the first argument.</dd>

> <dt><pre><b>IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt</pre></dt>
<dd>Checks if the second argument is weakly equivalent to the first argument.</dd>

> <dt><pre><b>IsGorenstein</b>(O::AlgEtQOrd)->BoolElt</pre></dt>
<dd>Checks if the order O is Gorenstein.</dd>


## List of instrinsics in AlgEtQ/IdealClassMonoid.m:

> <dt><pre><b>ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre></dt>
<dd>returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S).</dd>

> <dt><pre><b>ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre></dt>
<dd>returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals.</dd>


