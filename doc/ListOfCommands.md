# List of instrinsics in AlgEt.m:
--

<pre>
<b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEt
</pre>

*Given a sequence of number fields returns the étale algebra corresponding to the direct product.*

<pre>
<b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEt
</pre>

*Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.*

<pre>
<b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEt
</pre>

*Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.*

<pre>
<b>EtaleAlgebra</b>(f::RngUPolElt[FldNum]) -> AlgEt
</pre>

*Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.*


# List of instrinsics in AlgEtAttributes.m:
--

<pre>
<b>Print</b>(A::AlgEt)
</pre>

*Prints the defining polynomial or the components defining A.*

<pre>
<b>DefiningPolynomial</b>(A::AlgEt) -> RngUPolElt
</pre>

*Returns the defining polynomial of A, if the corresponding number fields are distinct.*

<pre>
<b>Components</b>(A::AlgEt) -> SeqEnum
</pre>

*Returns the number fields of which A is a product of,together with embeddings and projections*

<pre>
<b>Dimension</b>(A::AlgEt)->RngInt
</pre>

*Dimension of A*

<pre>
<b>AbsoluteDimension</b>(A::AlgEt)->RngInt
</pre>

*Dimension of A over the prime field*

<pre>
<b>HasBaseField</b>(A::AlgEt) -> BoolElt,FldNum
</pre>

*Returns whether A has common base field. If this is the case it returns it.*

<pre>
<b>BaseField</b>(A::AlgEt) -> FldNum
</pre>

*Returns the common base field of the Algebra, if it exists.*

<pre>
<b>PrimeField</b>(A::AlgEt) -> FldNum
</pre>

*Returns the prime field of the Algebra*

<pre>
<b>'eq'</b>(A1::AlgEt,A2::AlgEt) -> BoolElt
</pre>

*A1 eq A2*


# List of instrinsics in Homs.m:
--

<pre>
<b>HomsToC</b>(A::AlgEt : Precision:=30)->SeqEnum[Map]
</pre>

*returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30*


# List of instrinsics in Elt.m:
--

<pre>
<b>Print</b>(x::AlgEtElt)
</pre>

*Print the AlgEtElt.*

<pre>
<b>Parent</b>(x::AlgEtElt) -> AlgEt
</pre>

*Returns the algebra to which the elemenet belongs to.*

<pre>
<b>Algebra</b>(x::AlgEtElt) -> AlgEt
</pre>

*Returns the algebra to which the elemenet belongs to.*

<pre>
<b>Components</b>(x::AlgEtElt) -> SeqEnum
</pre>

*Given an element x returns its components, which are elements of number fields.*

<pre>
<b>AbsoluteCoordinates</b>(x::AlgEtElt) -> SeqEnum
</pre>

*Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.*

<pre>
<b>IsCoercible</b>(A::AlgEt, x::.) -> BoolElt, .
</pre>

*Return whether x is coercible into A and the result of the coercion if so.*

<pre>
<b>'!'</b>(A::AlgEt, x::.) -> AlgEtElt
</pre>

*Coerce x into A.*

<pre>
<b>One</b>(A::AlgEt) -> AlgEtElt
</pre>

*The multiplicative neutral element of A.*

<pre>
<b>Zero</b>(A::AlgEt) -> AlgEtElt
</pre>

*The additive neutral element of A.*

<pre>
<b>IsUnit</b>(x::AlgEtElt) -> BoolElt
</pre>

*Returns wheter x is a unit in A.*

<pre>
<b>IsZeroDivisor</b>(x::AlgEtElt) -> BoolElt
</pre>

*Returns wheter x is a not unit in A.*

<pre>
<b>IsZeroDivisor2</b>(x::AlgEtElt) -> BoolElt
</pre>

*Returns wheter x is a not unit in A.*

<pre>
<b>Random</b>(A::AlgEt , bd::RngIntElt) -> AlgEtElt
</pre>

*Random element of A. The Coefficients are bounded by the positive integer bd.*

<pre>
<b>Random</b>(A::AlgEt : bd:=3) -> AlgEtElt
</pre>

*Random element of A. The Coefficients are bounded by VarArg bd (default 3).*

<pre>
<b>RandomUnit</b>(A::AlgEt , bd::RngIntElt) -> AlgEtElt
</pre>

*Random unit of A. The Coefficients are bounded by the positive integer bd.*

<pre>
<b>'eq'</b>(x1::AlgEtElt,x2::AlgEtElt) -> BoolElt
</pre>

*Is x1=x2 ?*

<pre>
<b>'eq'</b>(x1::RngIntElt,x2::AlgEtElt) -> BoolElt
</pre>

*Is x1=x2 ?*

<pre>
<b>'eq'</b>(x1::FldRatElt,x2::AlgEtElt) -> BoolElt
</pre>

*Is x1=x2 ?*

<pre>
<b>'eq'</b>(x1::AlgEtElt,x2::RngIntElt) -> BoolElt
</pre>

*Is x1=x2 ?*

<pre>
<b>'eq'</b>(x1::AlgEtElt,x2::FldRatElt) -> BoolElt
</pre>

*Is x1=x2 ?*

<pre>
<b>'+'</b>(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::.,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::AlgEtElt,x2::.) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::RngIntElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::FldRatElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::AlgEtElt,x2::RngIntElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'+'</b>(x1::AlgEtElt,x2::FldRatElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'-'</b>(x::AlgEtElt) -> AlgEtElt
</pre>

*-x*

<pre>
<b>'-'</b>(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1-x2*

<pre>
<b>'-'</b>(x1::.,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1-x2*

<pre>
<b>'-'</b>(x1::AlgEtElt,x2::.) -> AlgEtElt
</pre>

*x1-x2*

<pre>
<b>'-'</b>(x1::RngIntElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'-'</b>(x1::FldRatElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'-'</b>(x1::AlgEtElt,x2::RngIntElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'-'</b>(x1::AlgEtElt,x2::FldRatElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'*'</b>(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1*x2*

<pre>
<b>'*'</b>(x1::.,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1*x2*

<pre>
<b>'*'</b>(x1::AlgEtElt,x2::.) -> AlgEtElt
</pre>

*x1*x2*

<pre>
<b>'*'</b>(x1::RngIntElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'*'</b>(x1::FldRatElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'*'</b>(x1::AlgEtElt,x2::RngIntElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'*'</b>(x1::AlgEtElt,x2::FldRatElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>Inverse</b>(x::AlgEtElt) -> AlgEtElt
</pre>

*1/x*

<pre>
<b>'^'</b>(x::AlgEtElt,n::RngIntElt) -> AlgEtElt
</pre>

*x^n*

<pre>
<b>'/'</b>(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1/x2*

<pre>
<b>'/'</b>(x1::.,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1/x2*

<pre>
<b>'/'</b>(x1::AlgEtElt,x2::.) -> AlgEtElt
</pre>

*x1/x2*

<pre>
<b>'/'</b>(x1::RngIntElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'/'</b>(x1::FldRatElt,x2::AlgEtElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'/'</b>(x1::AlgEtElt,x2::RngIntElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'/'</b>(x1::AlgEtElt,x2::FldRatElt) -> AlgEtElt
</pre>

*x1+x2*

<pre>
<b>'&+'</b>(seq::SeqEnum[AlgEtElt]) -> AlgEtElt
</pre>

*Given a sequence of AlgEtElt returns the sum of the entries.*

<pre>
<b>'&*'</b>(seq::SeqEnum[AlgEtElt]) -> AlgEtElt
</pre>

*Given a sequence of AlgEtElt returns the product of the entries.*

<pre>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtElt],bs::SeqEnum[AlgEtElt]) -> AlgEtElt
</pre>

*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<pre>
<b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtElt]) -> AlgEtElt
</pre>

*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<pre>
<b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtElt]) -> AlgEtElt
</pre>

*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<pre>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtElt],bs::SeqEnum[RngIntElt]) -> AlgEtElt
</pre>

*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<pre>
<b>SumOfProducts</b>(as::SeqEnum[AlgEtElt],bs::SeqEnum[FldRatElt]) -> AlgEtElt
</pre>

*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

<pre>
<b>MinimalPolynomial</b>(x::AlgEtElt) -> RngUPolElt
</pre>

*Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.*

<pre>
<b>MinimalPolynomial</b>(x::AlgEtElt, F::Rng) -> RngUPolElt
</pre>

*Returns the minimal polynommial over the ring F of the element x.*

<pre>
<b>AbsoluteMinimalPolynomial</b>(x::AlgEtElt) -> RngUPolElt
</pre>

*Returns the minimal polynommial over the prime field of the element x.*

<pre>
<b>IsIntegral</b>(x::AlgEtElt) -> BoolElt
</pre>

*Returns whether the element x is integral (over the integers).*

<pre>
<b>Evaluate</b>(f::RngUPolElt,a::AlgEtElt) -> AlgEtElt
</pre>

*Evaluate the polynomial f at the element a.*

<pre>
<b>PrimitiveElement</b>(A::AlgEt) -> AlgEtElt
</pre>

*Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.*

<pre>
<b>PowerBasis</b>(A::AlgEt) -> SeqEnum[AlgEtElt]
</pre>

*Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A*

<pre>
<b>Basis</b>(A::AlgEt) -> SeqEnum
</pre>

*Returns a basis of the algebra over the common base field.*

<pre>
<b>AbsoluteBasis</b>(A::AlgEt) -> SeqEnum
</pre>

*Returns a basis of the algebra over the prime field.*

<pre>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtElt] , basis::SeqEnum[AlgEtElt]) -> SeqEnum
</pre>

*Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.*

<pre>
<b>OrthogonalIdempotents</b>(A::AlgEt) -> SeqEnum
</pre>

*Returns the orthogonal ideampotent element of the étale algebra A*

<pre>
<b>Idempotents</b>(A::AlgEt) -> SeqEnum
</pre>

*Returns the ideampotent element of the étale algebra A*

<pre>
<b>CoordinatesOverBaseRing</b>(seq::SeqEnum[AlgEtElt] , basis::SeqEnum[AlgEtElt]) -> SeqEnum
</pre>

*//TODO I am not sure what this function should do*

<pre>
<b>'.'</b>(A::AlgEt,i::RngIntElt) -> AlgEtElt
</pre>

*A.i returns the ith generator of the Algebra A*


# List of instrinsics in TraceNorm.m:
--

<pre>
<b>Trace</b>(x::AlgEtElt) -> Any
</pre>

*Returns the trace of the element x of an étale algebra.*

<pre>
<b>Norm</b>(x::AlgEtElt) -> Any
</pre>

*Returns the norm of the element x of an étale algebra.*

<pre>
<b>AbsoluteTrace</b>(x::AlgEtElt) -> Any
</pre>

*Returns the absolute trace of the element x of an étale algebra.*

<pre>
<b>AbsoluteNorm</b>(x::AlgEtElt) -> Any
</pre>

*Returns the absolute norm of the element x of an étale algebra.*

<pre>
<b>TraceDualIdeal</b>(I::AlgEtIdl) -> AlgEtIdl
</pre>

*Returns the trace dual ideal of an ideal in an order in an etale algebra.*

<pre>
<b>TraceDualIdeal</b>(O::AlgEtOrd) -> AlgEtIdl
</pre>

*Returns the trace dual ideal of an order in an etale algebra.*


# List of instrinsics in Ord.m:
--

<pre>
<b>Print</b>(A::AlgEtOrd)
</pre>

*Print the order.*

<pre>
<b>IsCoercible</b>(S::AlgEtOrd, x::.) -> BoolElt, Any
</pre>

*Return whether x is coercible into S and the result if so.*

<pre>
<b>Order</b>( gens::SeqEnum[AlgEtElt] : Check:=100 ) -> AlgEtOrd
</pre>

*Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped.*

<pre>
<b>Order</b>(A::AlgEt , orders::Tup) -> AlgEtOrd
</pre>

*Given a sequence of order in the number fiedls defining the etale algebra A, generates the product order.*

<pre>
<b>OrderOver</b>( gens::SeqEnum[AlgEtElt] , S::RngOrd : Check:=100 ) -> AlgEtOrd
</pre>

*Given elements gens from an etale algebra A over a base field F and an order S in F, returns the order in A generated by gens over S. The parameter Check (default <true,100>) determines whether the programs checks if gens create a multiplicatively closed lattice, and if not adds elements until it is so, for Check[2] times.*

<pre>
<b>Algebra</b>(S::AlgEtOrd) -> AlgEt
</pre>

*Returns the algebra of the order.*

<pre>
<b>myHash</b>(S::AlgEtOrd)->SeqEnum[RngInt]
</pre>

*hash function for AlgEtOrd*

<pre>
<b>ZBasis</b>(S::AlgEtOrd)->SeqEnum[AlgEtElt]
</pre>

*Return a Z-basis of the order.*

<pre>
<b>Generators</b>(S::AlgEtOrd)->SeqEnum[AlgEtElt]
</pre>

*Return a set of generators of the order.*

<pre>
<b>'eq'</b>(O1::AlgEtOrd,O2::AlgEtOrd)->BoolElt
</pre>

*checks equality of orders in an etale Algebra*

<pre>
<b>'in'</b>(x::AlgEtElt,O::AlgEtOrd) -> BoolElt
</pre>

*inclusion of elements*

<pre>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtElt],O::AlgEtOrd) -> SeqEnum
</pre>

*AbsoluteCoordinates with respect to the ZBasis*

<pre>
<b>'in'</b>(x::RngIntElt,O::AlgEtOrd) -> BoolElt
</pre>

*inclusion of elements*

<pre>
<b>'in'</b>(x::FldRatElt,O::AlgEtOrd) -> BoolElt
</pre>

*inclusion of elements*

<pre>
<b>One</b>(S::AlgEtOrd)->AlgEtElt
</pre>

*unit element of S*

<pre>
<b>Zero</b>(S::AlgEtOrd)->AlgEtElt
</pre>

*zero element of S*

<pre>
<b>Random</b>(O::AlgEtOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtElt
</pre>

*Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false*

<pre>
<b>Random</b>(O::AlgEtOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtElt
</pre>

*Returns a random (small coefficient) element of O. 
  The range of the random coefficients can be increased by giving the optional argument CoeffRange.
  One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false*

<pre>
<b>EquationOrder</b>(A::AlgEt) -> AlgEtOrd
</pre>

*Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial*

<pre>
<b>ProductOfEquationOrders</b>(A::AlgEt)->AlgEtOrd
</pre>

*Given a product of number field A, returns the order consisting of the product of the equation orders of the number fields.*

<pre>
<b>MaximalOrder</b>(A::AlgEt)->AlgEtOrd
</pre>

*Returns the maximal order of the étale algebra A.*

<pre>
<b>IsMaximal</b>(S::AlgEtOrd) -> BoolElt
</pre>

*Returns wheter the given order is the maximal order of the étale algebra.*

<pre>
<b>IsProductOfOrders</b>(O::AlgEtOrd)->BoolElt, Tup
</pre>

*Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.*

<pre>
<b>Index</b>(T::AlgEtOrd) -> FldRatElt
</pre>

*given an order T computes its index with respect to the basis of the algebra of T as a free Z-module*

<pre>
<b>Index</b>(S::AlgEtOrd, T::AlgEtOrd) -> Any
</pre>

*given two orders T \subset S, returns [S:T] = #S/T*

<pre>
<b>'subset'</b>(O1 :: AlgEtOrd, O2 :: AlgEtOrd) -> BoolElt
</pre>

*Checks if the first argument is inside the second.*

<pre>
<b>'*'</b>(O1::AlgEtOrd,O2::AlgEtOrd)->BoolElt
</pre>

*checks equality of orders in an etale Algebra*

<pre>
<b>'meet'</b>(O1::AlgEtOrd,O2::AlgEtOrd)->BoolElt
</pre>

*checks equality of orders in an etale Algebra*

<pre>
<b>MultiplicatorRing</b>(R::AlgEtOrd) -> AlgEtOrd
</pre>

*Returns the MultiplicatorRing of an order R, that is R itself.*

<pre>
<b>ListToSequence</b>(L::List)->SeqEnum
</pre>

*given a list of elements returns the same sequence*

<pre>
<b>Discriminant</b>(R::AlgEtOrd) -> RngInt
</pre>

*returns the discriminant of the order*


# List of instrinsics in Quotients.m:
--

<pre>
<b>Quotient</b>(I::AlgEtIdl, zbJ::SeqEnum[AlgEtElt]) -> GrpAb, Map
</pre>

*Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order.*

<pre>
<b>Quotient</b>(I::AlgEtIdl, J::AlgEtIdl) -> GrpAb, Map
</pre>

*given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J*

<pre>
<b>Quotient</b>(S::AlgEtOrd, zbJ::SeqEnum[AlgEtElt]) -> GrpAb, Map
</pre>

*Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.*

<pre>
<b>ResidueRing</b>(S::AlgEtOrd,I::AlgEtIdl) -> GrpAb , Map
</pre>

*given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S)*

<pre>
<b>ResidueField</b>(P::AlgEtIdl) -> FldFin, Map
</pre>

*given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.*

<pre>
<b>PrimitiveElementResidueField</b>(P::AlgEtIdl)->AlgEtElt
</pre>

*Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^*.*

<pre>
<b>QuotientVS</b>(I::AlgEtIdl, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
</pre>

*Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)*

<pre>
<b>QuotientVS</b>(I::AlgEtOrd, J::AlgEtOrd, P::AlgEtIdl) -> ModRng, Map
</pre>

*Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)*

<pre>
<b>QuotientVS</b>(I::AlgEtOrd, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
</pre>

*Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)*

<pre>
<b>QuotientVS</b>(I::AlgEtIdl, J::AlgEtOrd, P::AlgEtIdl) -> ModRng, Map
</pre>

*Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)*

<pre>
<b>QuotientVS</b>(I::AlgEtIdl, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
</pre>

*let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image)*


# List of instrinsics in OverOrders.m:
--

<pre>
<b>MinimalOverOrders</b>(R::AlgEtOrd : singular_primes := [], orders := {@ @}) -> SetIndx[AlgEtOrd]
</pre>

*@ @*

<pre>
<b>FindOverOrders_Minimal</b>(R::AlgEtOrd) -> SetIndx[AlgEtOrd]
</pre>

*Given an order R returns all the over orders by a recursive search of the minimal overordes.
  Based on "On the computations of overorders" by TommyHofmann and Carlo Sircana*

<pre>
<b>FindOverOrders</b>(E::AlgEtOrd:  populateoo_in_oo := false) -> SetIndx[AlgEtOrd]
</pre>

*returns all the overorders of E, and populates*

<pre>
<b>pMaximalOrder</b>(O::AlgEtOrd, p::RngIntElt) -> AlgEtOrd
</pre>

*given O, retuns the maximal p over order*

<pre>
<b>FindOverOrders</b>(E::AlgEtOrd, O::AlgEtOrd) -> SetIndx[AlgEtOrd]
</pre>

*given E subset O, returns the sequence of orders between E and O*

<pre>
<b>FindOverOrders_Naive</b>(E::AlgEtOrd) -> SetIndx[AlgEtOrd]
</pre>

*returns all the overorders of E*

<pre>
<b>FindOverOrders_Naive</b>(E::AlgEtOrd, O::AlgEtOrd) -> SetIndx[AlgEtOrd]
</pre>

*given E subset O, returns the sequence of orders between E and O*


# List of instrinsics in Idl.m:
--

<pre>
<b>Ideal</b>(S::AlgEtOrd, gens::SeqEnum) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gens.*

<pre>
<b>Ideal</b>(S::AlgEtOrd, idls::Tup) -> AlgEtIdl
</pre>

*Given an order S which is a product of orders S_i in the number fiedls generting the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal direct product of the I_i.*

<pre>
<b>Ideal</b>(S::AlgEtOrd, gen::Any) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(S::AlgEtOrd, gen::AlgEtElt) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(gen::AlgEtElt, S::AlgEtOrd) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(S::AlgEtOrd, gen::RngIntElt) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(gen::RngIntElt, S::AlgEtOrd) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(S::AlgEtOrd, gen::FldRatElt) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>'*'</b>(gen::FldRatElt, S::AlgEtOrd) -> AlgEtIdl
</pre>

*Creates an ideal of S, generated by gen.*

<pre>
<b>Print</b>(I::AlgEtIdl)
</pre>

*Prints the ideal.*

<pre>
<b>'!!'</b>(T::AlgEtOrd,I::AlgEtIdl) -> AlgEtIdl
</pre>

*Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is in S, then IT=I*

<pre>
<b>Algebra</b>(I::AlgEtIdl) -> AlgEt
</pre>

*Returns the étale algebra in which the ideal lives.*

<pre>
<b>Order</b>(I::AlgEtIdl) -> AlgEtOrd
</pre>

*Returns the order of definition of the ideal.*

<pre>
<b>ZBasis</b>(I::AlgEtIdl)->SeqEnum[AlgEtElt]
</pre>

*Returns a Z-basis of the ideal.*

<pre>
<b>Generators</b>(I::AlgEtIdl) -> SeqEnum[AlgEtElt]
</pre>

*Returns the generators of the ideal.*

<pre>
<b>myHash</b>(I::AlgEtIdl)->RngInt
</pre>

*Hash function*

<pre>
<b>'eq'</b>(I::AlgEtIdl , J::AlgEtIdl ) -> BoolElt
</pre>

*Equality testing.*

<pre>
<b>'ne'</b>(I::AlgEtIdl , J::AlgEtIdl ) -> BoolElt
</pre>

*Equality testing*

<pre>
<b>'eq'</b>(I::AlgEtIdl, S::AlgEtOrd) -> BoolElt
</pre>

*Return if I eq S. I needs to be an ideal of S.*

<pre>
<b>'eq'</b>(S::AlgEtOrd,I::AlgEtIdl) -> BoolElt
</pre>

*Return if I eq S. I needs to be an ideal of S.*

<pre>
<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtElt],I::AlgEtIdl) -> SeqEnum
</pre>

*AbsoluteCoordiantes with respect to the ZBasis*

<pre>
<b>'in'</b>(x::AlgEtElt , I::AlgEtIdl ) -> BoolElt
</pre>

*Returns if x is in I.*

<pre>
<b>'in'</b>(x::RngIntElt , I::AlgEtIdl ) -> BoolElt
</pre>

*Returns if x is in I.*

<pre>
<b>'in'</b>(x::FldRatElt , I::AlgEtIdl ) -> BoolElt
</pre>

*Returns if x is in I.*

<pre>
<b>'subset'</b>(S::AlgEtOrd,I::AlgEtIdl) -> BoolElt
</pre>

*given an ideal I of S, return if S subseteq I*

<pre>
<b>'subset'</b>(I::AlgEtIdl,S::AlgEtOrd) -> BoolElt
</pre>

*given an ideal I of S, return if I subseteq S*

<pre>
<b>'subset'</b>(I1 :: AlgEtIdl, I2 :: AlgEtIdl) -> BoolElt
</pre>

*Checks if the first argument is inside the second. The ideals need to be fractional*

<pre>
<b>Index</b>(T::AlgEtIdl) -> FldRatElt
</pre>

*Given an ideal T computes its index with respect to the basis of the algebra of T as a free Z-module.*

<pre>
<b>Index</b>(J::AlgEtIdl, I::AlgEtIdl) -> Any
</pre>

*Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I]*

<pre>
<b>Index</b>(S::AlgEtOrd, I::AlgEtIdl) -> Any
</pre>

*given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I]*

<pre>
<b>OneIdeal</b>(S::AlgEtOrd) -> AlgEtIdl
</pre>

*Given an order S returns the ideal 1*S which will be cached*

<pre>
<b>Conductor</b>(O::AlgEtOrd) ->AlgEtOrdIdl
</pre>

*computes the conductor of an order*

<pre>
<b>'+'</b>(I::AlgEtIdl , J::AlgEtIdl ) -> AlgEtIdl
</pre>

*returns the sum of two ideals*

<pre>
<b>'*'</b>(I::AlgEtIdl , J::AlgEtIdl ) -> AlgEtIdl
</pre>

*Product of two ideals.*

<pre>
<b>'*'</b>(I::AlgEtIdl , x::AlgEtElt ) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'*'</b>(x::AlgEtElt, I::AlgEtIdl) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'*'</b>(x::RngIntElt, I::AlgEtIdl) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'*'</b>(I::AlgEtIdl, x::RngIntElt) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'*'</b>(x::FldRatElt, I::AlgEtIdl) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'*'</b>(I::AlgEtIdl, x::FldRatElt) -> AlgEtIdl
</pre>

*Returns x*I.*

<pre>
<b>'^'</b>(I::AlgEtIdl, n::RngIntElt) -> AlgEtIdl
</pre>

*nth power of an ideal.*

<pre>
<b>'meet'</b>(I::AlgEtIdl, S::AlgEtOrd) -> AlgEtIdl
</pre>

*given an ideal I of S, return S cap I*

<pre>
<b>'meet'</b>(S::AlgEtOrd,I::AlgEtIdl) -> AlgEtIdl
</pre>

*given an ideal I of S, return S cap I*

<pre>
<b>'meet'</b>(I::AlgEtIdl, J::AlgEtIdl) -> AlgEtIdl
</pre>

*given ideals I and J, return J cap I*

<pre>
<b>'&+'</b>(seq::SeqEnum[AlgEtIdl])->AlgEtIdl
</pre>

*Returns the sum of the fractional ideals in the sequence.*

<pre>
<b>ColonIdeal</b>(I::AlgEtIdl,J::AlgEtIdl)->AlgEtIdl
</pre>

*Computes the colon ideal (I:J) (as an O-ideal) of two O-idealsi.*

<pre>
<b>ColonIdeal</b>(O::AlgEtOrd,J::AlgEtIdl)->AlgEtIdl
</pre>

*computes the colon ideal (1*O:J) (as an O-ideal)*

<pre>
<b>ColonIdeal</b>(I::AlgEtIdl,O::AlgEtOrd)->AlgEtIdl
</pre>

*computes the colon ideal (I:1*O) (as an O-ideal)*

<pre>
<b>IsInvertible</b>(I::AlgEtIdl) ->BoolElt
</pre>

*checks if the ideal I is invertible in its order of definition O*

<pre>
<b>Inverse</b>(I::AlgEtIdl) ->AlgEtIdl
</pre>

*computes the inverse of an ideal of a maximal order*

<pre>
<b>MultiplicatorRing</b>(I::AlgEtIdl) -> AlgEtOrd
</pre>

*Given a fractional R-ideal I computes its multiplicator ring (I:I). If the overorders of R are known the corresponding overorder is returned, in order to preserve the known attributes.*

<pre>
<b>IsProductOfIdeals</b>(I::AlgEtIdl) -> BoolElt, Tup
</pre>

*Return if the argument is a product of ideals in number fields, and if so return also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).*

<pre>
<b>Random</b>(I::AlgEtIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtElt
</pre>

*Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false*

<pre>
<b>Random</b>(I::AlgEtIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtElt
</pre>

*Returns a random (small coefficient) element of I. 
  The range of the random coefficients can be increased by giving the optional argument CoeffRange.
  One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false*

<pre>
<b>IsCoprime</b>(I::AlgEtIdl,J::AlgEtIdl) -> BoolElt
</pre>

*given two integral ideals I and J of an order S, returns whether I+J=R*

<pre>
<b>IsIntegral</b>(I::AlgEtIdl) -> BoolElt
</pre>

*returns wheter the ideal I of S is integral, that is I \subseteq S*

<pre>
<b>MakeIntegral</b>(I::AlgEtIdl) -> AlgEtIdl,RngIntElt
</pre>

*given a fractional S ideal I, returns the ideal d*I,d when d is the smallest integer such that d*I is integral in S*

<pre>
<b>MinimalInteger</b>(I::AlgEtIdl) -> RngIntElt
</pre>

*returns the smallest integer contained in the ideal I*

<pre>
<b>CoprimeRepresentative</b>(I::AlgEtIdl,J::AlgEtIdl) -> AlgEtElt,AlgEtIdl
</pre>

*return an element x such that x*I is an integral ideal coprime with J, togheter with the product x*I.. The first ideal must be invertible and the second should be integral*


# List of instrinsics in Completion.m:
--

<pre>
<b>Completion</b>(P::AlgEtIdl : MinPrecision:=20) -> FldPad,Map
</pre>

*Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage (it acts a bit weird).*


# List of instrinsics in ComplexConj.m:
--

<pre>
<b>HasComplexConjugate</b>(A::AlgEt) -> BoolElt
</pre>

*Returns if the algebra is the product of CM fields*

<pre>
<b>ComplexConjugate</b>(x::AlgEtElt) -> AlgEtElt
</pre>

*If A is a product of CM fields, it returns the complex conjugate of the argument*

<pre>
<b>IsConjugateStable</b>(O::AlgEtOrd) -> BoolElt,AlgEtOrd
</pre>

*Returns wheter O is conjugate stable and the complex conjugate.*

<pre>
<b>ComplexConjugate</b>(O::AlgEtOrd) -> AlgEtOrd
</pre>

*It returns the complex conjugate of the argument.*

<pre>
<b>IsConjugateStable</b>(I::AlgEtIdl) -> BoolElt,AlgEtIdl
</pre>

*Returns wheter O is conjugate stable and the complex conjugate.*

<pre>
<b>ComplexConjugate</b>(I::AlgEtIdl) -> AlgEtIdl
</pre>

*if A is a product of CM fields, it returns the complex conjugate of the argument*


# List of instrinsics in ComplexMult.m:
--

<pre>
<b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtCMType
</pre>

*given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType*

<pre>
<b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtCMType
</pre>

*given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType*

<pre>
<b>CMType</b>( b::AlgEtElt  ) -> AlgEtCMType
</pre>

*given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive*

<pre>
<b>CreateCMType</b>( b::AlgEtElt  ) -> AlgEtCMType
</pre>

*given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive*

<pre>
<b>Print</b>( PHI :: AlgEtCMType)
</pre>

*print the AlgEtCMType*

<pre>
<b>CMPositiveElement</b>( PHI::AlgEtCMType )->AlgEtElt
</pre>

*given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI)*

<pre>
<b>CMPosElt</b>( PHI::AlgEtCMType )->AlgEtElt
</pre>

*given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI)*

<pre>
<b>Homs</b>( PHI::AlgEtCMType : prec:=30 )->SeqEnum[Map]
</pre>

*given a AlgEtCMType PHI returns the sequence of maps to CC defining it*

<pre>
<b>'eq'</b>(PHI1 :: AlgEtCMType, PHI2::AlgEtCMType : prec:=30)->BoolElt
</pre>

*returns whether two cm types are equal*

<pre>
<b>Precision</b>(PHI :: AlgEtCMType)->RngIntElt
</pre>

*the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision)*

<pre>
<b>ChangePrecision</b>(PHI0 :: AlgEtCMType, prec::RngIntElt )->AlgEtCMType
</pre>

*changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision)*

<pre>
<b>ChangePrecision</b>(~PHI :: AlgEtCMType, prec::RngIntElt )
</pre>

*changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision)*

<pre>
<b>AllCMTypes</b>(A::AlgEt : Precision := 30 ) -> SeqEnum[AlgEtCMType]
</pre>

*Returns all the AlgEtCMTypes of A*


# List of instrinsics in IntermediateIdeals.m:
--

<pre>
<b>MinimalIntermediateIdeals</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<pre>
<b>IntermediateIdeals</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones*

<pre>
<b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. 
  They are produced recursively using from the minimal ones*

<pre>
<b>MaximalIntermediateIdeals</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<pre>
<b>MaximalIntermediateIdeals</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

<pre>
<b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtIdl,J::AlgEtIdl, O::AlgEtOrd)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. 
  Note that we need O subset (J:J).
  They are produced recursively using from the maximal ones*

<pre>
<b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtIdl,J::AlgEtIdl, O::AlgEtOrd)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S.
  Note that we need O subset (J:J).
  They are produced recursively using from the maximal ones.*

<pre>
<b>MinimalIntermediateIdealsVS</b>(I::AlgEtIdl,J::AlgEtIdl : primes:=[])->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I.*

<pre>
<b>IntermediateIdealsVS</b>(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
</pre>

*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones*


# List of instrinsics in ShortEltSmallRep.m:
--

<pre>
<b>ShortestElement</b>(I::AlgEtIdl) ->AlgEtElt
</pre>

*Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). Note that if I decomposes into a direct sums I=I1+I2, then the element returned is the sum of the shortest elements of I1 and I2.*

<pre>
<b>SmallRepresentative</b>(I::AlgEtIdl) ->AlgEtIdl,AlgEtElt
</pre>

*Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortestElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.*


# List of instrinsics in MinimalGenerators.m:
--

<pre>
<b>TwoGeneratingSet</b>(I::AlgEtIdl)
</pre>

*A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.*


# List of instrinsics in CRT.m:
--

<pre>
<b>ChineseRemainderTheorem</b>(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
</pre>

*Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.*

<pre>
<b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtIdl],as::SeqEnum[AlgEtElt])-> AlgEtElt
</pre>

*Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.*

<pre>
<b>ChineseRemainderTheorem</b>(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
</pre>

*Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.*


# List of instrinsics in PicardGroup.m:
--

<pre>
<b>ResidueRingUnits</b>(S::AlgEtOrd,I::AlgEtIdl) -> GrpAb,Map
</pre>

*returns the group (S/I)^* and a map (S/I)^* -> S. It is required S to be maximal*

<pre>
<b>IsPrincipal</b>(I1::AlgEtIdl : GRH:=false )->BoolElt, AlgAssElt
</pre>

*Return if the argument is a principal ideal; if so the function returns also the generator.
    The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*

<pre>
<b>PicardGroup</b>( S::AlgEtOrd : GRH:=false ) -> GrpAb, Map
</pre>

*Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes
    The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".*

<pre>
<b>UnitGroup</b>(S::AlgEtOrd : GRH:=false ) -> GrpAb, Map
</pre>

*Return the unit group of a order in a etale algebra
    The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".*

<pre>
<b>IsIsomorphic</b>(I::AlgEtIdl, J::AlgEtIdl : GRH:=false ) -> BoolElt, AlgAssElt
</pre>

*Checks if I=x*J, for some x. If so, also x is returned
    The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false"*


# List of instrinsics in FactPrimes.m:
--

<pre>
<b>Factorization</b>(I::AlgEtIdl) -> Tup
</pre>

*Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.*

<pre>
<b>PrimesAbove</b>(I::AlgEtIdl) -> SeqEnum[AlgAssEtOrdIdl]
</pre>

*Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.*

<pre>
<b>IsPrime</b>(I::AlgEtIdl) -> BoolElt
</pre>

*Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal*

<pre>
<b>IsBassAtPrime</b>(S::AlgEtOrd,P::AlgEtIdl) -> BoolElt
</pre>

*Check if the order is Bass at the prime ideal P.*

<pre>
<b>IsBass</b>(S::AlgEtOrd) -> BoolElt
</pre>

*Check if the order is Bass.*

<pre>
<b>IsGorensteinAtPrime</b>(S::AlgEtOrd,P::AlgEtIdl) -> BoolElt
</pre>

*Check if the order is Gorenstein at the prime ideal P.*


# List of instrinsics in TotRealTotPos.m:
--

<pre>
<b>IsTotallyReal</b>(a::AlgEtElt) -> BoolElt
</pre>

*returns whther a is totally real*

<pre>
<b>IsTotallyRealPositive</b>(a::AlgEtElt) -> BoolElt
</pre>

*returns whether a is totally positive, that is, totally real and with positive image in C*

<pre>
<b>TotallyRealSubAlgebra</b>(K::AlgEt) -> AlgEt,Map
</pre>

*Given a CM algebra K returns the unique totally real subalgebra, and an embedding*

<pre>
<b>TotallyRealUnitGroup</b>(S::AlgEtOrd) -> Grp
</pre>

*Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^**

<pre>
<b>TotallyRealPositiveUnitGroup</b>(S::AlgEtOrd) -> Grp
</pre>

*Given an order S in a CM étale algebra A
      returns the groups of totally positive units of S, as a subgroup of S^**


# List of instrinsics in PrintSave.m:
--

<pre>
<b>PrintSeqAlgEtElt</b>(seq::SeqEnum[AlgEtElt]) -> SeqEnum,MonStgElt
</pre>

*Given a sequence of elements of an AlgEt, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.*

<pre>
<b>PrintWKICM</b>(R::AlgEtOrd) -> MonStgElt
</pre>

*Given an order R in an AlgEt, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM*

<pre>
<b>LoadWKICM</b>(str::MonStgElt) -> AlgEtOrd
</pre>

*Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics.*


# List of instrinsics in LowCohenMacaulayType.m:
--

<pre>
<b>NonGorensteinPrimes</b>(S::AlgEtOrd)->SeqEnum,SeqEnum
</pre>

*Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S*

<pre>
<b>CohenMacaulayType</b>(S::AlgEtOrd)->RngIntElt
</pre>

*Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.*


# List of instrinsics in WkClasses.m:
--

<pre>
<b>WKICM_bar</b>(S::AlgEtOrd : Method:="Auto") -> SeqEnum
</pre>

*Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*

<pre>
<b>WKICM</b>(E::AlgEtOrd : Method:="Auto")->SeqEnum
</pre>

*Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*

<pre>
<b>WKICM_bar_intermediate_idls</b>(S::AlgEtOrd) -> SeqEnum[AlgEtIdl]
</pre>

*???*


# List of instrinsics in WkTesting.m:
--

<pre>
<b>IsWeakEquivalent</b>(I::AlgEtIdl,J::AlgEtIdl)->BoolElt
</pre>

*checks if 1 \in (I:J)*(J:I). This function does not require that the ideals are defined over the same order.*

<pre>
<b>IsWeakEquivalent</b>(O1::AlgEtOrd,O2::AlgEtOrd)->BoolElt
</pre>

*check if the two orders are weakly equivalent, that is equal*

<pre>
<b>IsWeakEquivalent</b>(O::AlgEtOrd,J::AlgEtIdl)->BoolElt
</pre>

*checks if the second argument is weakly equivalent to the first argument*

<pre>
<b>IsWeakEquivalent</b>(J::AlgEtIdl,O::AlgEtOrd)->BoolElt
</pre>

*checks if the second argument is weakly equivalent to the first argument*

<pre>
<b>IsGorenstein</b>(O::AlgEtOrd)->BoolElt
</pre>

*checks if the order O is Gorenstein*


# List of instrinsics in IdealClassMonoid.m:
--

<pre>
<b>ICM_bar</b>(S::AlgEtOrd : GRH:=false ) -> SeqEnum
</pre>

*returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S)*

<pre>
<b>ICM</b>(S::AlgEtOrd : GRH:=false ) -> SeqEnum
</pre>

*returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals*


