## List of instrinsics in AlgEtQ/AlgEt.m:

`<b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ`
*Given a sequence of number fields returns the étale algebra corresponding to the direct product.*

`<b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ`
*Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.*

`<b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ`
*Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.*


## List of instrinsics in AlgEtQ/AlgEtAttributes.m:

`<b>Print</b>(A::AlgEtQ)`
*Prints the defining polynomial or the components defining A.*

`<b>DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt`
*Returns the defining polynomial of A, if the corresponding number fields are distinct.*

`<b>Components</b>(A::AlgEtQ) -> SeqEnum`
*Returns the number fields of which A is a product of,together with embeddings and projections.*

`<b>Dimension</b>(A::AlgEtQ)->RngInt`
*Dimension of A.*

`<b>AbsoluteDimension</b>(A::AlgEtQ)->RngInt`
*Dimension of A over the prime field.*

`<b>HasBaseField</b>(A::AlgEtQ) -> BoolElt,FldNum`
*Returns whether A has common base field. If this is the case it returns it.*

`<b>BaseField</b>(A::AlgEtQ) -> FldNum`
*Returns the common base field of the Algebra, if it exists.*

`<b>PrimeField</b>(A::AlgEtQ) -> FldNum`
*Returns the prime field of the Algebra.*

`<b>'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt`
*A1 eq A2.*


## List of instrinsics in AlgEtQ/Homs.m:

`<b>HomsToC</b>(A::AlgEtQ : Precision:=30)->SeqEnum[Map]`
*returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30*


## List of instrinsics in AlgEtQ/Elt.m:

`<b>Print</b>(x::AlgEtQElt)`
*Print the AlgEtQElt.*

`<b>Parent</b>(x::AlgEtQElt) -> AlgEtQ`
*Returns the algebra to which the elemenet belongs to.*

`<b>Algebra</b>(x::AlgEtQElt) -> AlgEtQ`
*Returns the algebra to which the elemenet belongs to.*

`<b>Components</b>(x::AlgEtQElt) -> SeqEnum`
*Given an element x returns its components, which are elements of number fields.*

`<b>AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum`
*Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.*

`<b>IsCoercible</b>(A::AlgEtQ, x::.) -> BoolElt, .`
*Return whether x is coercible into A and the result of the coercion if so.*

`<b>'!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt`
*Coerce x into A.*

`<b>One</b>(A::AlgEtQ) -> AlgEtQElt`
*The multiplicative neutral element of A.*

`<b>Zero</b>(A::AlgEtQ) -> AlgEtQElt`
*The additive neutral element of A.*

`<b>IsUnit</b>(x::AlgEtQElt) -> BoolElt`
*Returns wheter x is a unit in A.*

`<b>IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt`
*Returns wheter x is a not unit in A.*

`<b>IsZeroDivisor2</b>(x::AlgEtQElt) -> BoolElt`
*Returns wheter x is a not unit in A.*

`<b>Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt`
*Random element of A. The Coefficients are bounded by the positive integer bd.*

`<b>Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt`
*Random element of A. The Coefficients are bounded by VarArg bd (default 3).*

`<b>RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt`
*Random unit of A. The Coefficients are bounded by the positive integer bd.*

`<b>'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt`
*Is x1=x2 ?*

`<b>'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt`
*Is x1=x2 ?*

`<b>'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt`
*Is x1=x2 ?*

`<b>'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt`
*Is x1=x2 ?*

`<b>'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt`
*Is x1=x2 ?*

`<b>'+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt`
*x1+x2.*

`<b>'+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt`
*x1+x2.*

`<b>'-'</b>(x::AlgEtQElt) -> AlgEtQElt`
*-x.*

`<b>'-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt`
*x1-x2.*

`<b>'-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt`
*x1-x2.*

`<b>'*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt`
*x1*x2.*

`<b>'*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt`
*x1*x2.*

`<b>Inverse</b>(x::AlgEtQElt) -> AlgEtQElt`
*1/x.*

`<b>'^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt`
*x^n.*

`<b>'/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt`
*x1/x2.*

`<b>'/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt`
*x1/x2.*

`<b>'&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt`
*Given a sequence of AlgEtQElt returns the sum of the entries.*

`<b>'&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt`
*Given a sequence of AlgEtQElt returns the product of the entries.*

`<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt`
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

`<b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt`
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

`<b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt`
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

`<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt`
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

`<b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt`
*Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.*

`<b>MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt`
*Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.*

`<b>MinimalPolynomial</b>(x::AlgEtQElt, F::Rng) -> RngUPolElt`
*Returns the minimal polynommial over the ring F of the element x.*

`<b>AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt`
*Returns the minimal polynommial over the prime field of the element x.*

`<b>IsIntegral</b>(x::AlgEtQElt) -> BoolElt`
*Returns whether the element x is integral (over the integers).*

`<b>Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt`
*Evaluate the polynomial f at the element a.*

`<b>PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt`
*Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.*

`<b>PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]`
*Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.*

`<b>Basis</b>(A::AlgEtQ) -> SeqEnum`
*Returns a basis of the algebra over the common base field.*

`<b>AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum`
*Returns a basis of the algebra over the prime field.*

`<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum`
*Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.*

`<b>OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum`
*Returns the orthogonal ideampotent element of the étale algebra A.*

`<b>Idempotents</b>(A::AlgEtQ) -> SeqEnum`
*Returns the ideampotent element of the étale algebra A.*


## List of instrinsics in AlgEtQ/TraceNorm.m:

`<b>Trace</b>(x::AlgEtQElt) -> Any`
*Returns the trace of the element x of an étale algebra.*

`<b>Norm</b>(x::AlgEtQElt) -> Any`
*Returns the norm of the element x of an étale algebra.*

`<b>AbsoluteTrace</b>(x::AlgEtQElt) -> Any`
*Returns the absolute trace of the element x of an étale algebra.*

`<b>AbsoluteNorm</b>(x::AlgEtQElt) -> Any`
*Returns the absolute norm of the element x of an étale algebra.*

`<b>TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl`
*Returns the trace dual ideal of an ideal in an order in an etale algebra.*

`<b>TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl`
*Returns the trace dual ideal of an order in an etale algebra.*


## List of instrinsics in AlgEtQ/Ord.m:

`<b>Print</b>(A::AlgEtQOrd)`
*Print the order.*

`<b>IsCoercible</b>(S::AlgEtQOrd, x::.) -> BoolElt, Any`
*Return whether x is coercible into S and the result if so.*

`<b>Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd`
*Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped. The vararg CheckIsKnownOrder determines if we check if the order is already known, i.e. in the attribute Algebra`KnownOrders. The default value is true.*

`<b>Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd`
*Given a sequence of order in the number fiedls defining the etale algebra A, generates the product order.*

`<b>Algebra</b>(S::AlgEtQOrd) -> AlgEtQ`
*Returns the algebra of the order.*

`<b>myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]`
*Hash function for AlgEtQOrd.*

`<b>ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]`
*Return a Z-basis of the order.*

`<b>Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]`
*Return a set of generators of the order.*

`<b>'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt`
*Checks equality of orders in an etale Algebra.*

`<b>'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt`
*Inclusion of elements.*

`<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum`
*AbsoluteCoordinates with respect to the ZBasis.*

`<b>'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt`
*Inclusion of elements.*

`<b>'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt`
*Inclusion of elements.*

`<b>One</b>(S::AlgEtQOrd)->AlgEtQElt`
*Unit element of S.*

`<b>Zero</b>(S::AlgEtQOrd)->AlgEtQElt`
*Zero element of S.*

`<b>Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt`
*Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

`<b>Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt`
*Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

`<b>IsKnownOrder</b>(~R::AlgEtQOrd)`
*This procedure checks wheter the order R is already in the list of known orders of the algebra A of definition of R. If so then it replaces R with the copy stored in the attribute KnownOrders. If not it adds it to KnownOrders. This is done to avoid creating multiple copies of the same order.*

`<b>EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd`
*Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial.*

`<b>ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd`
*Given a product of number field A, returns the order consisting of the product of the equation orders of the number fields.*

`<b>MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd`
*Returns the maximal order of the étale algebra A.*

`<b>IsMaximal</b>(S::AlgEtQOrd) -> BoolElt`
*Returns wheter the given order is the maximal order of the étale algebra.*

`<b>IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup`
*Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.*

`<b>Index</b>(T::AlgEtQOrd) -> FldRatElt`
*Given an order T computes its index with respect to the basis of the algebra of T as a free Z-module.*

`<b>Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any`
*Given two orders T \subset S, returns [S:T] = #S/T.*

`<b>'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt`
*Checks if the first argument is inside the second.*

`<b>'*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd`
*Returns the order generated by the orders O1 and O2.*

`<b>'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd`
*Intersection of orders.*

`<b>MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd`
*Returns the multiplicator ring of an order R, that is R itself.*


## List of instrinsics in AlgEtQ/Quotients.m:

`<b>Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map`
*Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order.*

`<b>Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map`
*Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.*

`<b>Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map`
*Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.*

`<b>ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map`
*Given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S). We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.*

`<b>ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map`
*Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.*

`<b>PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt`
*Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^*.*

`<b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map`
*Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

`<b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map`
*Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

`<b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map`
*Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*

`<b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map`
*Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image).*


## List of instrinsics in AlgEtQ/OverOrders.m:

`<b>IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt`
*Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.*

`<b>MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]`
*Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

`<b>MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]`
*Computes the minimal overorders of R.*

`<b>OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]`
*Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.*

`<b>OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]`
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*

`<b>FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]`
*We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.*


## List of instrinsics in AlgEtQ/GraphOverOrders.m:

`<b>GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir`
*Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].*


## List of instrinsics in AlgEtQ/Idl.m:

`<b>Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl`
*Creates an ideal of S, generated by gens.*

`<b>Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl`
*Given an order S which is a product of orders S_i in the number fiedls generting the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal direct product of the I_i.*

`<b>Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>'*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl`
*Creates an ideal of S, generated by gen.*

`<b>Print</b>(I::AlgEtQIdl)`
*Prints the ideal.*

`<b>'!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl`
*Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is in S, then IT=I.*

`<b>Algebra</b>(I::AlgEtQIdl) -> AlgEtQ`
*Returns the étale algebra in which the ideal lives.*

`<b>Order</b>(I::AlgEtQIdl) -> AlgEtQOrd`
*Returns the order of definition of the ideal.*

`<b>ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]`
*Returns a Z-basis of the ideal.*

`<b>Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]`
*Returns the generators of the ideal.*

`<b>myHash</b>(I::AlgEtQIdl)->RngInt`
*Hash function.*

`<b>'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt`
*Equality testing.*

`<b>'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt`
*Equality testing.*

`<b>'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt`
*Return if I eq S. I needs to be an ideal of S.*

`<b>'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt`
*Return if I eq S. I needs to be an ideal of S.*

`<b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum`
*AbsoluteCoordiantes with respect to the ZBasis.*

`<b>'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt`
*Returns if x is in I.*

`<b>'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt`
*Returns if x is in I.*

`<b>'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt`
*Returns if x is in I.*

`<b>'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt`
*Given an ideal I of S, return if S subseteq I.*

`<b>'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt`
*Given an ideal I of S, return if I subseteq S.*

`<b>'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt`
*Checks if the first argument is inside the second. The ideals need to be fractional.*

`<b>Index</b>(T::AlgEtQIdl) -> FldRatElt`
*Given an ideal T computes its index with respect to the basis of the algebra of T as a free Z-module.*

`<b>Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any`
*Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I].*

`<b>Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any`
*Given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I].*

`<b>OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl`
*Given an order S returns the ideal 1*S which will be cached.*

`<b>Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl`
*Computes the conductor of an order, defined as he colon ideal (O:OK), where OK is the maximal order of the algebra.*

`<b>'+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl`
*Returns the sum of two ideals.*

`<b>'*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl`
*Product of two ideals.*

`<b>'*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl`
*Returns x*I.*

`<b>'*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl`
*Returns x*I.*

`<b>'*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl`
*Returns x*I.*

`<b>'*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl`
*Returns x*I.*

`<b>'*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl`
*Returns x*I.*

`<b>'*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl`
*Returns x*I.*

`<b>'^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl`
*nth power of an ideal.*

`<b>'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl`
*Given an ideal I of S, return S cap I.*

`<b>'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl`
*Given an ideal I of S, return S cap I.*

`<b>'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl`
*Given ideals I and J, return J cap I.*

`<b>'&+'</b>(seq::SeqEnum[AlgEtQIdl])->AlgEtQIdl`
*Returns the sum of the fractional ideals in the sequence.*

`<b>ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl`
*Computes the colon ideal (I:J) (as an O-ideal) of two O-idealsi.*

`<b>ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl`
*Computes the colon ideal (1*O:J) (as an O-ideal).*

`<b>ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl`
*Computes the colon ideal (I:1*O) (as an O-ideal).*

`<b>IsInvertible</b>(I::AlgEtQIdl) ->BoolElt`
*Checks if the ideal I is invertible in its order of definition O.*

`<b>Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl`
*Computes the inverse of an ideal of a maximal order.*

`<b>MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd`
*Given a fractional R-ideal I computes its multiplicator ring (I:I). If the overorders of R are known the corresponding overorder is returned, in order to preserve the known attributes.*

`<b>IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup`
*Return if the argument is a product of ideals in number fields, and if so return also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).*

`<b>Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt`
*Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

`<b>Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt`
*Returns a random (small coefficient) element of I. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.*

`<b>IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt`
*Given two integral ideals I and J of an order S, returns whether I+J=R.*

`<b>IsIntegral</b>(I::AlgEtQIdl) -> BoolElt`
*Returns wheter the ideal I of S is integral, that is I \subseteq S.*

`<b>MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt`
*Given a fractional S ideal I, returns the ideal d*I,d when d is the smallest integer such that d*I is integral in S.*

`<b>MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt`
*Returns the smallest integer contained in the ideal I.*

`<b>CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl`
*Returns an element x such that x*I is an integral ideal coprime with J, togheter with the product x*I. The first ideal must be invertible and the second should be integral.*


## List of instrinsics in AlgEtQ/ZBasisLLL.m:

`<b>ZBasisLLL</b>(S::AlgEtQOrd)`
*A procedure that replaces the ZBasis with an LLL-reduced one.*

`<b>ZBasisLLL</b>(S::AlgEtQIdl)`
*A procedure that replaces the ZBasis with an LLL-reduced one.*


## List of instrinsics in AlgEtQ/Completion.m:

`<b>Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map`
*Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage (it acts a bit weird).*


## List of instrinsics in AlgEtQ/ComplexConj.m:

`<b>HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt`
*Returns if the algebra is the product of CM fields.*

`<b>ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt`
*If A is a product of CM fields, it returns the complex conjugate of the argument.*

`<b>IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd`
*Returns wheter O is conjugate stable and the complex conjugate.*

`<b>ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd`
*It returns the complex conjugate of the argument.*

`<b>IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl`
*Returns wheter O is conjugate stable and the complex conjugate.*

`<b>ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl`
*If A is a product of CM fields, it returns the complex conjugate of the argument.*


## List of instrinsics in AlgEtQ/ComplexMult.m:

`<b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType`
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

`<b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType`
*Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.*

`<b>CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType`
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

`<b>CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType`
*Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.*

`<b>Print</b>( PHI :: AlgEtQCMType)`
*Print the AlgEtQCMType.*

`<b>CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt`
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

`<b>CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt`
*Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).*

`<b>Homs</b>( PHI::AlgEtQCMType : prec:=30 )->SeqEnum[Map]`
*Given a AlgEtQCMType PHI returns the sequence of maps to CC defining it.*

`<b>'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=30)->BoolElt`
*Returns whether two cm types are equal. This happens if and only if the ration of (any) two CMPositiveElements is totally real and totally positive.*

`<b>Precision</b>(PHI :: AlgEtQCMType)->RngIntElt`
*Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

`<b>ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType`
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

`<b>ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )`
*Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).*

`<b>AllCMTypes</b>(A::AlgEtQ : Precision := 30 ) -> SeqEnum[AlgEtQCMType]`
*Returns all the AlgEtQCMTypes of A.*


## List of instrinsics in AlgEtQ/IntermediateIdeals.m:

`<b>MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

`<b>IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.*

`<b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. They are produced recursively using from the minimal ones.*

`<b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

`<b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.*

`<b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

`<b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.*

`<b>IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]`
*Given ideals J subset I over the same order, and a positive integer N, it returns all the ideals K such that J subset K subset I and [I:K]=N. These are produced by recursively searching for maximal submodules.*

`<b>MinimalIntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl : primes:=[])->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I.*

`<b>IntermediateIdealsVS</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]`
*Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones*


## List of instrinsics in AlgEtQ/IdealsOfIndex.m:

`<b>IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]`
*Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.*

`<b>IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]`
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

`<b>IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]`
*Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.*

`<b>IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]`
*Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".*

`<b>IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]`
*Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".*


## List of instrinsics in AlgEtQ/ShortEltSmallRep.m:

`<b>ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt`
*Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).*

`<b>SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt`
*Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.*


## List of instrinsics in AlgEtQ/MinimalGenerators.m:

`<b>TwoGeneratingSet</b>(I::AlgEtQIdl)`
*A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.*


## List of instrinsics in AlgEtQ/CRT.m:

`<b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt`
*Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.*

`<b>ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt`
*Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.*


## List of instrinsics in AlgEtQ/PicardGroup.m:

`<b>ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map`
*Returns the group (S/I)^* and a map (S/I)^* -> S. It is required S to be maximal.*

`<b>IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgAssElt`
*Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*

`<b>PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map`
*Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".*

`<b>UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map`
*Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".*

`<b>IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgAssElt`
*Checks if I=x*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".*


## List of instrinsics in AlgEtQ/FactPrimes.m:

`<b>Factorization</b>(I::AlgEtQIdl) -> Tup`
*Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.*

`<b>PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgAssEtOrdIdl]`
*Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.*

`<b>SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgAssEtOrdIdl]`
*Returns the non-invertible primes of the order.*

`<b>NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx`
*Returns the non-invertible primes of the order.*

`<b>IsPrime</b>(I::AlgEtQIdl) -> BoolElt`
*Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal.*

`<b>IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt`
*Check if the order is Bass at the prime ideal P.*

`<b>IsBass</b>(S::AlgEtQOrd) -> BoolElt`
*Check if the order is Bass.*

`<b>IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt`
*Check if the order is Gorenstein at the prime ideal P.*


## List of instrinsics in AlgEtQ/TotRealTotPos.m:

`<b>IsTotallyReal</b>(a::AlgEtQElt) -> BoolElt`
*Returns whther a is totally real.*

`<b>IsTotallyRealPositive</b>(a::AlgEtQElt) -> BoolElt`
*Returns whether a is totally positive, that is, totally real and with positive image in C.*

`<b>TotallyRealSubAlgebra</b>(K::AlgEtQ) -> AlgEtQ,Map`
*Given a CM algebra K returns the unique totally real subalgebra, and an embedding.*

`<b>TotallyRealUnitGroup</b>(S::AlgEtQOrd) -> Grp`
*Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^*.*

`<b>TotallyRealPositiveUnitGroup</b>(S::AlgEtQOrd) -> Grp`
*Given an order S in a CM étale algebra A. Returns the groups of totally positive units of S, as a subgroup of S^*.*


## List of instrinsics in AlgEtQ/PrintSave.m:

`<b>PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt`
*Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.*

`<b>PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt`
*Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.*

`<b>LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd`
*Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics.*


## List of instrinsics in AlgEtQ/LowCohenMacaulayType.m:

`<b>NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum`
*Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.*

`<b>CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt`
*Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P*S^t where S^t is the trace dual of S.*

`<b>CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt`
*Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.*


## List of instrinsics in AlgEtQ/WkClasses.m:

`<b>WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum`
*Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*

`<b>WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum`
*Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.*


## List of instrinsics in AlgEtQ/WkTesting.m:

`<b>IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt`
*Checks if I and J are weakly equivalent 1 \in (I:J)*(J:I). This function does not require that the ideals are defined over the same order.*

`<b>IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt`
*Check if the two orders are weakly equivalent, that is equal.*

`<b>IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt`
*Checks if the second argument is weakly equivalent to the first argument.*

`<b>IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt`
*Checks if the second argument is weakly equivalent to the first argument.*

`<b>IsGorenstein</b>(O::AlgEtQOrd)->BoolElt`
*Checks if the order O is Gorenstein.*


## List of instrinsics in AlgEtQ/IdealClassMonoid.m:

`<b>ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum`
*returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S).*

`<b>ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum`
*returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals.*


