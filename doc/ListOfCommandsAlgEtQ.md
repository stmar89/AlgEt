## List of instrinsics in AlgEtQ/AlgEt.m:

> <pre><b>EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ</pre>
<em>Given a sequence of number fields returns the étale algebra corresponding to the direct product. Note: the number fields with DefiningPolynomial of degree one should be created with the vararg DoLinearExtention set to true.</em>

> <pre><b>EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ</pre>
<em>Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.</em>

> <pre><b>EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ</pre>
<em>Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.</em>


## List of instrinsics in AlgEtQ/AlgEtAttributes.m:

> <pre><b>Print</b>(A::AlgEtQ)</pre>
<em>Prints the defining polynomial or the components defining A.</em>

> <pre><b>DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt</pre>
<em>Returns the defining polynomial of A, if the corresponding number fields are distinct.</em>

> <pre><b>Components</b>(A::AlgEtQ) -> SeqEnum,SeqEnum,SeqEnum</pre>
<em>Returns the number fields of which A is a product of, together with embeddings and projections.</em>

> <pre><b>Dimension</b>(A::AlgEtQ)->RngInt</pre>
<em>Dimension of A.</em>

> <pre><b>AbsoluteDimension</b>(A::AlgEtQ)->RngInt</pre>
<em>Dimension of A over the prime field.</em>

> <pre><b>HasBaseField</b>(A::AlgEtQ) -> BoolElt,FldNum</pre>
<em>Returns whether A has common base field. If this is the case it returns it.</em>

> <pre><b>BaseField</b>(A::AlgEtQ) -> FldNum</pre>
<em>Returns the common base field of the Algebra, if it exists.</em>

> <pre><b>PrimeField</b>(A::AlgEtQ) -> FldNum</pre>
<em>Returns the prime field of the Algebra.</em>

> <pre><b>'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt</pre>
<em>Equality testing.</em>


## List of instrinsics in AlgEtQ/NF.m:

> <pre><b>IsNumberField</b>(A::AlgEtQ) -> BoolElt,FldNum,Map</pre>
<em>Given an étale algebra over Q returns wheter it is a number field, and if so the number field and an isomorphism from the étale algebra to the number field.</em>


## List of instrinsics in AlgEtQ/Elt.m:

> <pre><b>Print</b>(x::AlgEtQElt)</pre>
<em>Print the element.</em>

> <pre><b>Parent</b>(x::AlgEtQElt) -> AlgEtQ</pre>
<em>Returns the algebra to which the element belongs to.</em>

> <pre><b>Algebra</b>(x::AlgEtQElt) -> AlgEtQ</pre>
<em>Returns the algebra to which the element belongs to.</em>

> <pre><b>Components</b>(x::AlgEtQElt) -> SeqEnum</pre>
<em>Given an element, returns its components, which are elements of number fields.</em>

> <pre><b>AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum</pre>
<em>Given an element, returns the coordinates relative to the absolute basis, which are elements of the prime rational field.</em>

> <pre><b>IsCoercible</b>(A::AlgEtQ, x::.) -> BoolElt, .</pre>
<em>Return whether the element is coercible into A and the result of the coercion if so.</em>

> <pre><b>'!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt</pre>
<em>Coerce x into A.</em>

> <pre><b>One</b>(A::AlgEtQ) -> AlgEtQElt</pre>
<em>The multiplicative neutral element of A.</em>

> <pre><b>Zero</b>(A::AlgEtQ) -> AlgEtQElt</pre>
<em>The additive neutral element of A.</em>

> <pre><b>IsUnit</b>(x::AlgEtQElt) -> BoolElt</pre>
<em>Returns wheter x is a unit in A.</em>

> <pre><b>IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt</pre>
<em>Returns wheter x is a not unit in A.</em>

> <pre><b>Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre>
<em>Random element of A. The Coefficients are bounded by the positive integer bd.</em>

> <pre><b>Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt</pre>
<em>Random element of A. The Coefficients are bounded by the VarArg bd (default 3).</em>

> <pre><b>RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre>
<em>Random unit of A. The Coefficients are bounded by the positive integer bd.</em>

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
<em>x1+x2.</em>

> <pre><b>'-'</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
<em>-x.</em>

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
<em>x1-x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
<em>x1\*x2.</em>

> <pre><b>Inverse</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
<em>1/x.</em>

> <pre><b>'^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt</pre>
<em>x^n.</em>

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
<em>x1/x2.</em>

> <pre><b>'&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
<em>Given a sequence of AlgEtQElt returns the sum of the entries.</em>

> <pre><b>'&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
<em>Given a sequence of AlgEtQElt returns the product of the entries.</em>

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
<em>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</em>

> <pre><b>SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
<em>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</em>

> <pre><b>SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
<em>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</em>

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt</pre>
<em>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</em>

> <pre><b>SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt</pre>
<em>Given sequences as and bs, such that #as eq #bs, returns &+[as[i]\*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.</em>

> <pre><b>MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
<em>Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.</em>

> <pre><b>MinimalPolynomial</b>(x::AlgEtQElt, F::Rng) -> RngUPolElt</pre>
<em>Returns the minimal polynommial over the ring F of the element x.</em>

> <pre><b>AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
<em>Returns the minimal polynommial over the prime field of the element x.</em>

> <pre><b>IsIntegral</b>(x::AlgEtQElt) -> BoolElt</pre>
<em>Returns whether the element x is integral (over the integers).</em>

> <pre><b>Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt</pre>
<em>Evaluate the polynomial f at the element a.</em>

> <pre><b>PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt</pre>
<em>Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.</em>

> <pre><b>PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]</pre>
<em>Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.</em>

> <pre><b>Basis</b>(A::AlgEtQ) -> SeqEnum</pre>
<em>Returns a basis of the algebra over the common base field.</em>

> <pre><b>AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum</pre>
<em>Returns a basis of the algebra over the prime field.</em>

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum</pre>
<em>Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.</em>

> <pre><b>OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
<em>Returns the orthogonal ideampotent element of the étale algebra A.</em>

> <pre><b>Idempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
<em>Returns the ideampotent element of the étale algebra A.</em>


## List of instrinsics in AlgEtQ/Homs.m:

> <pre><b>HomsToC</b>(A::AlgEtQ : Prec:=Precision(GetDefaultRealField()))->SeqEnum[Map]</pre>
<em>Returns the sequence of homomorphisms from A to the complex field CC. The precision of CC is given by the optional parameter "Prec". Default value is 30</em>

> <pre><b>Hom</b>(A::AlgEtQ , B::AlgEtQ , img::SeqEnum[AlgEtQElt] : CheckMultiplicative:=false, CheckUnital:=false, ComputeInverse:=true)->Map</pre>
<em>Given two étale algebras A and B and a sequence img of elements of B, returns the Q-algebra homomorphism defined by sending the AbsoluteBasis of A to img. The VarArg CheckMultiplicative determines if the multiplicativity of the defined map is checked, while the VarArg CheckUnital determines wheter One(A) is sent to One(B). If the VarArg ComputeInverse is true, it checkes whether the map is invertible and, if so, it defines also the inverse (by assigning preimages).</em>

> <pre><b>NaturalAction</b>(K::AlgEtQ, V::AlgEtQ)->Map</pre>
<em>Let K=K1x...Kn be a product of distinct number fields, and s1,...,sn be strinctly positive integers. Put V=K1^s1x...xKn^sn. It returns the natural action of K on V, that is, the componentwise diagonal.</em>

> <pre><b>DiagonalEmbedding</b>(K::AlgEtQ, V::AlgEtQ)->Map</pre>
<em>Let K=K1x...Kn be a product of distinct number fields, and s1,...,sn be strinctly positive integers. Put V=K1^s1x...xKn^sn. It returns the natural action of K on V, that is, the componentwise diagonal.</em>


## List of instrinsics in AlgEtQ/DirectProduct.m:

> <pre><b>DirectProduct</b>(seq::SeqEnum[AlgEtQ]) -> AlgEtQ,SeqEnum[Map],SeqEnum[Map]</pre>
<em>Given a sequence of étale algebras, it returns the direct product,togheter with canonical inclusions and projections.</em>


## List of instrinsics in AlgEtQ/Ord.m:

> <pre><b>Print</b>(A::AlgEtQOrd)</pre>
<em>Print the order.</em>

> <pre><b>IsCoercible</b>(S::AlgEtQOrd, x::.) -> BoolElt, Any</pre>
<em>Return whether x is coercible into S and the result if so.</em>

> <pre><b>Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd</pre>
<em>Construnct the order generated by gens over the rationals. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the generators. If Check 0 then this step is skipped. The vararg CheckIsKnownOrder determines whether the program checks if the order is already known, i.e. in the attribute KnownOrders of the algebra. This is to avoid the creation of multiple copies of the same order. The default value is true.</em>

> <pre><b>Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd</pre>
<em>Given a sequence of order in the number fields defining the étale algebra A, generates the direct sum order.</em>

> <pre><b>Algebra</b>(S::AlgEtQOrd) -> AlgEtQ</pre>
<em>Returns the algebra of the order.</em>

> <pre><b>myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]</pre>
<em>Hash function for AlgEtQOrd.</em>

> <pre><b>ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
<em>Return a Z-basis of the order.</em>

> <pre><b>Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
<em>Return a set of generators (as a Z-algebra) of the order.</em>

> <pre><b>'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
<em>Checks equality of orders in an étale Algebra.</em>

> <pre><b>'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt</pre>
<em>Inclusion of element in the order.</em>

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum</pre>
<em>AbsoluteCoordinates with respect to the ZBasis.</em>

> <pre><b>'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt</pre>
<em>Inclusion of elements.</em>

> <pre><b>'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt</pre>
<em>Inclusion of elements.</em>

> <pre><b>One</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
<em>Unit element of S.</em>

> <pre><b>Zero</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
<em>Zero element of S.</em>

> <pre><b>Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
<em>Random element of O. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which is set to false by default.</em>

> <pre><b>Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
<em>Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which is set to false by default.</em>

> <pre><b>IsKnownOrder</b>(~R::AlgEtQOrd)</pre>
<em>This procedure checks wheter the order R is already in the list of known orders of the algebra A containing R. If so then it replaces R with the copy stored in the attribute KnownOrders. If not it adds it to KnownOrders. This is done to avoid creating multiple copies of the same order.</em>

> <pre><b>EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd</pre>
<em>Given an étale algebra defined by a polynomial, returns the monogenic order defined by the same polynomial.</em>

> <pre><b>ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd</pre>
<em>Given an étale algebra A, returns the order consisting of the product of the equation orders of the number fields.</em>

> <pre><b>MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd</pre>
<em>Returns the maximal order of the étale algebra A. It is the direct sum of the ring of integers of the number fields composing the algebra.</em>

> <pre><b>IsMaximal</b>(S::AlgEtQOrd) -> BoolElt</pre>
<em>Returns wheter the given order is the maximal order of the étale algebra.</em>

> <pre><b>IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup</pre>
<em>Return if the argument is a product of orders in number fields, and if so return also the sequence of these orders.</em>

> <pre><b>Index</b>(T::AlgEtQOrd) -> FldRatElt</pre>
<em>Given an order T computes its index with respect to the basis of the algebra of T as a free Z-module.</em>

> <pre><b>Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any</pre>
<em>Given two orders T subset S, returns [S:T] = #S/T.</em>

> <pre><b>'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt</pre>
<em>Checks if the first argument is inside the second.</em>

> <pre><b>'*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
<em>Returns the order generated by the orders O1 and O2.</em>

> <pre><b>'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
<em>Intersection of orders.</em>

> <pre><b>MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd</pre>
<em>Returns the multiplicator ring of an order R, that is R itself.</em>


## List of instrinsics in AlgEtQ/Idl.m:

> <pre><b>Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gens.</em>

> <pre><b>Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl</pre>
<em>Given an order S which is a product of orders S_i in the number fields generating the Algebra(S), and a Tup of ideals I_i of S_i, returns the S-ideal corresponding to the direct sum of the I_i.</em>

> <pre><b>Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>'*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Creates an ideal of S, generated by gen.</em>

> <pre><b>Print</b>(I::AlgEtQIdl)</pre>
<em>Prints the ideal.</em>

> <pre><b>'!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Given an S-ideal I and an order T, returns the extension IT as a T-ideal. Note that if T is a subset of S, then IT=I.</em>

> <pre><b>Algebra</b>(I::AlgEtQIdl) -> AlgEtQ</pre>
<em>Returns the étale algebra in which the ideal lives.</em>

> <pre><b>Order</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
<em>Returns the order of definition of the ideal.</em>

> <pre><b>ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]</pre>
<em>Returns a Z-basis of the ideal.</em>

> <pre><b>Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre>
<em>Returns the generators of the ideal.</em>

> <pre><b>myHash</b>(I::AlgEtQIdl)->RngInt</pre>
<em>Hash function.</em>

> <pre><b>'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre>
<em>Equality testing.</em>

> <pre><b>'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre>
<em>Equality testing.</em>

> <pre><b>'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt</pre>
<em>Return if I eq S. I needs to be an ideal of S.</em>

> <pre><b>'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
<em>Return if I eq S. I needs to be an ideal of S.</em>

> <pre><b>AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum</pre>
<em>AbsoluteCoordiantes with respect to the ZBasis.</em>

> <pre><b>'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt</pre>
<em>Returns if x is in I.</em>

> <pre><b>'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt</pre>
<em>Returns if x is in I.</em>

> <pre><b>'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt</pre>
<em>Returns if x is in I.</em>

> <pre><b>'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
<em>Given an ideal I of S, return if S subseteq I.</em>

> <pre><b>'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt</pre>
<em>Given an ideal I of S, return if I subseteq S.</em>

> <pre><b>'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt</pre>
<em>Checks if the first argument is inside the second. The ideals need to be fractional.</em>

> <pre><b>Index</b>(T::AlgEtQIdl) -> FldRatElt</pre>
<em>Given an ideal T computes its index with respect to the basis of the algebra of T as a free Q-module.</em>

> <pre><b>Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any</pre>
<em>Given fractional ideals J and I defined over the same order returns [J:I] = [J:J cap I]/[I : J cap I].</em>

> <pre><b>Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any</pre>
<em>Given and ideal I of an order S returns [S:I] = [S:S cap I]/[I : S cap I].</em>

> <pre><b>OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Given an order S returns the ideal 1\*S which will be cached.</em>

> <pre><b>Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl</pre>
<em>Computes the conductor of an order, defined as he colon ideal (O:OK), where OK is the maximal order of the algebra.</em>

> <pre><b>'+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
<em>Returns the sum of two ideals.</em>

> <pre><b>'*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
<em>Product of two ideals.</em>

> <pre><b>'*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl</pre>
<em>Returns x\*I.</em>

> <pre><b>'^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl</pre>
<em>nth power of an ideal.</em>

> <pre><b>'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Given an ideal I of S, return S cap I.</em>

> <pre><b>'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Given an ideal I of S, return S cap I.</em>

> <pre><b>'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Given ideals I and J, return J cap I.</em>

> <pre><b>'&+'</b>(seq::SeqEnum[AlgEtQIdl])->AlgEtQIdl</pre>
<em>Returns the sum of the fractional ideals in the sequence.</em>

> <pre><b>ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl</pre>
<em>Computes the colon ideal (I:J) (as an O-ideal) of two O-ideals, which is the set of elements x of the algebra such that x\*J subset I.</em>

> <pre><b>ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl</pre>
<em>Computes the colon ideal (1\*O:J) (as an O-ideal).</em>

> <pre><b>ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl</pre>
<em>Computes the colon ideal (I:1\*O) (as an O-ideal).</em>

> <pre><b>IsInvertible</b>(I::AlgEtQIdl) ->BoolElt</pre>
<em>Checks if the ideal I is invertible in its order of definition O.</em>

> <pre><b>Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl</pre>
<em>Computes the inverse of an invertible ideal.</em>

> <pre><b>MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
<em>Given a fractional ideal I computes its multiplicator ring (I:I).</em>

> <pre><b>IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup</pre>
<em>Return if the argument is a product of ideals in the number fields defining the algebra. If so, it returns also the sequence of these ideals (in the appropriate orders). Note: we require the Order(I) to be the MultiplicatorRing(I).</em>

> <pre><b>Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
<em>Random element of I. The Coefficients are bounded by the positive integer bd. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</em>

> <pre><b>Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
<em>Returns a random (small coefficient) element of I. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which by default is set to false.</em>

> <pre><b>IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt</pre>
<em>Given two integral ideals I and J of an order S, returns whether I+J=R.</em>

> <pre><b>IsIntegral</b>(I::AlgEtQIdl) -> BoolElt</pre>
<em>Returns wheter the ideal I of S is integral, that is I subseteq S.</em>

> <pre><b>MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt</pre>
<em>Given a fractional S ideal I, returns the ideal d\*I,d when d is the smallest integer such that d\*I is integral in S. Compare with SmallRepresentative.</em>

> <pre><b>MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt</pre>
<em>Returns the smallest integer contained in the ideal I.</em>

> <pre><b>CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl</pre>
<em>Returns an element x such that x\*I is an integral ideal coprime with J, togheter with the product x\*I. The first ideal must be invertible and the second should be integral.</em>


## List of instrinsics in AlgEtQ/ZBasisLLL.m:

> <pre><b>ZBasisLLL</b>(S::AlgEtQOrd)</pre>
<em>A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.</em>

> <pre><b>ZBasisLLL</b>(S::AlgEtQIdl)</pre>
<em>A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.</em>


## List of instrinsics in AlgEtQ/Quotients.m:

> <pre><b>Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
<em>Given an ideal I and the ZBasis of an ideal or order J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.</em>

> <pre><b>Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map</pre>
<em>Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.</em>

> <pre><b>Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
<em>Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.</em>

> <pre><b>ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map</pre>
<em>Given an integral ideal I of S, returns the abelian group S/I and the quotient map q:S -> S/I (with preimages). Important: the domain of q is the Algebra of S, since the elements of S are expressed as elements of A. We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.</em>

> <pre><b>ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map</pre>
<em>Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.</em>

> <pre><b>PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt</pre>
<em>Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^\*.</em>

> <pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
<em>Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).</em>

> <pre><b>QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
<em>Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).</em>

> <pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
<em>Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).</em>

> <pre><b>QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
<em>Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with preimages).</em>


## List of instrinsics in AlgEtQ/OverOrders.m:

> <pre><b>IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt</pre>
<em>Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.</em>

> <pre><b>MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]</pre>
<em>Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.</em>

> <pre><b>MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]</pre>
<em>Computes the minimal overorders of R.</em>

> <pre><b>OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]</pre>
<em>Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.</em>

> <pre><b>OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]</pre>
<em>We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.</em>

> <pre><b>FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]</pre>
<em>We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.</em>


## List of instrinsics in AlgEtQ/GraphOverOrders.m:

> <pre><b>GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir</pre>
<em>Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].</em>


## List of instrinsics in AlgEtQ/TraceNorm.m:

> <pre><b>Trace</b>(x::AlgEtQElt) -> Any</pre>
<em>Returns the trace of the element x of an étale algebra.</em>

> <pre><b>Norm</b>(x::AlgEtQElt) -> Any</pre>
<em>Returns the norm of the element x of an étale algebra.</em>

> <pre><b>AbsoluteTrace</b>(x::AlgEtQElt) -> Any</pre>
<em>Returns the absolute trace of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Trace.</em>

> <pre><b>AbsoluteNorm</b>(x::AlgEtQElt) -> Any</pre>
<em>Returns the absolute norm of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Norm.</em>

> <pre><b>TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>Returns the trace dual ideal of the ideal I, that is, the set of elements x of the algebra such that Trace(x\*I) is integer-valued.</em>

> <pre><b>TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl</pre>
<em>Returns the trace dual ideal of an order in an etale algebra, that is, the set of elements x of the algebra such that Trace(x\*O) is integer-valued.</em>


## List of instrinsics in AlgEtQ/Completion.m:

> <pre><b>Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map</pre>
<em>Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage.</em>


## List of instrinsics in AlgEtQ/ComplexConj.m:

> <pre><b>HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt</pre>
<em>Returns if the algebra is the product of CM fields.</em>

> <pre><b>ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
<em>If A is a product of CM fields, it returns the complex conjugate of the argument.</em>

> <pre><b>IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd</pre>
<em>Given an order O in a CM-étale algebra, it returns wheter O is conjugate stable and the complex conjugate.</em>

> <pre><b>ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd</pre>
<em>Given an order O in a CM-étale algebra, it returns the complex conjugate of the argument.</em>

> <pre><b>IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl</pre>
<em>Given a fractional ideal I in a CM-étale algebra, it returns wheter I is conjugate stable and the complex conjugate. Note: if the order of I is not conjugate stable, then the second output will be defined over the complex conjugate of the order.</em>

> <pre><b>ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
<em>If A is a product of CM fields, it returns the complex conjugate of the fractional ideal I. Note: if the order of I is not conjugate stable, then the output will be defined over the complex conjugate of the order.</em>


## List of instrinsics in AlgEtQ/ComplexMult.m:

> <pre><b>CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
<em>Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType</em>

> <pre><b>CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
<em>Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.</em>

> <pre><b>CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
<em>Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive, that is, Im(phi(b))>0 for every phi in PHI.</em>

> <pre><b>CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
<em>Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.</em>

> <pre><b>Print</b>( PHI :: AlgEtQCMType)</pre>
<em>Print the AlgEtQCMType.</em>

> <pre><b>CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
<em>Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).</em>

> <pre><b>CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
<em>Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).</em>

> <pre><b>Homs</b>( PHI::AlgEtQCMType : prec:=Precision(GetDefaultRealField()) )->SeqEnum[Map]</pre>
<em>Given a AlgEtQCMType PHI returns the sequence of maps to the complex field. The vararg prec (default value 30) determines the precision of the codomains of the maps.</em>

> <pre><b>'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=Precision(GetDefaultRealField()))->BoolElt</pre>
<em>Returns whether two cm types are equal. This happens if and only if the quotient of (any) two CMPositiveElements is totally real and totally positive.</em>

> <pre><b>Precision</b>(PHI :: AlgEtQCMType)->RngIntElt</pre>
<em>Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</em>

> <pre><b>ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType</pre>
<em>Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</em>

> <pre><b>ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )</pre>
<em>Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).</em>

> <pre><b>AllCMTypes</b>(A::AlgEtQ : Prec := Precision(GetDefaultRealField()) ) -> SeqEnum[AlgEtQCMType]</pre>
<em>Returns all the AlgEtQCMTypes of A. The vararg Prec determined the precision of the codomain of the maps defining the CMTypes.</em>


## List of instrinsics in AlgEtQ/IntermediateIdeals.m:

> <pre><b>MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals J subseteq I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subsetneq K subseteq I. Note J is never in the output.</em>

> <pre><b>IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I. They are produced recursively from the minimal ones. The output includes I and J.</em>

> <pre><b>IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I and (K:K)=S. They are produced recursively from the minimal ones. The output includes I, if (I:I)=S, and J, if (J:J)=S.</em>

> <pre><b>MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals J subseteq I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subseteq K subsetneq I. Note I is never in the output, while J is in the output if and only if the S-module I/J is simple, in which case the output consists only of J.</em>

> <pre><b>IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals I and J and an order O such that 
- S subset O,  
- J subseteq I, and 
- O subset (I:I),
returns all the fractional S-ideals K such that 
- J subseteq K subseteq I, and 
- O!!K = I. 
Note that the output always contains I. The output is produced by a recursiv use of MaximalIntermediateIdeals.</em>

> <pre><b>IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
<em>Given fractional S-ideals I and J and an order O such that 
- S subseteq O,  
- J subseteq I, and 
- O subseteq (I:I), 
returns all the fractional S-ideals K satisfying
- J subseteq K subseteq I, 
- O!!K = I, and 
- (K:K) eq S. 
In particular, the output contains I if and only if O = (I:I) = S. The output is produced by a recursive use of MaximalIntermediateIdeals.</em>

> <pre><b>IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]</pre>
<em>Given ideals J subseteq I over the same order, and a positive integer N, it returns all the ideals K such that 
- J subseteq K subseteq I, and 
- [I:K]=N. 
The output is produced by a recursive use of MaximalIntermediateIdeals.</em>


## List of instrinsics in AlgEtQ/IdealsOfIndex.m:

> <pre><b>IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
<em>Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.</em>

> <pre><b>IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
<em>Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.</em>

> <pre><b>IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]</pre>
<em>Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.</em>

> <pre><b>IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
<em>Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".</em>

> <pre><b>IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
<em>Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".</em>


## List of instrinsics in AlgEtQ/ShortEltSmallRep.m:

> <pre><b>ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt</pre>
<em>Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by randomly picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).</em>

> <pre><b>SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt</pre>
<em>Given a fractional R-ideal I, it returns an isomorphic ideal a\*I, and the element a, such that a\*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.</em>


## List of instrinsics in AlgEtQ/MinimalGenerators.m:

> <pre><b>TwoGeneratingSet</b>(I::AlgEtQIdl)</pre>
<em>A procedure that given an invertible ideal I put in the attibute Generators of I two non-zerodivisors in I that generate I. If I is known to be principal nothing is done.</em>


## List of instrinsics in AlgEtQ/CRT.m:

> <pre><b>ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt</pre>
<em>Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.</em>

> <pre><b>ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt</pre>
<em>Given a sequence `Is` of ideals of S, pairwise coprime, and a sequence `as` of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.</em>

> <pre><b>ChineseRemainderTheoremFunctions</b>(Is::SeqEnum[AlgEtQIdl])-> AlgEtQElt</pre>
<em>Given a sequence `Is` of N ideals of S, pairwise coprime, returns a function S->S^N representing the natural isomorphism S/&\*(Is) -> \prod_(I in Is) S/I and a function S^N-S representing the inverse.</em>


## List of instrinsics in AlgEtQ/PicardGroup.m:

> <pre><b>ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map</pre>
<em>Returns the group (S/I)^\* and a map (S/I)^\* -> S. The MultiplicatorRing(I) must be the maximal order.</em>

> <pre><b>ResidueRingUnits</b>(I::AlgEtQIdl) -> GrpAb,Map</pre>
<em>Returns the group (S/I)^\* and a map (S/I)^\* -> S, where S=Order(I) and the multiplicator ring of I is maximal.</em>

> <pre><b>ResidueRingUnitsSubgroupGenerators</b>(F::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre>
<em>Returns generators of (S/F)^\* where F is an ideal of the order S.</em>

> <pre><b>IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgEtQElt</pre>
<em>Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".</em>

> <pre><b>PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
<em>Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".</em>

> <pre><b>UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
<em>Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".</em>

> <pre><b>IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgEtQElt</pre>
<em>Checks if I=x\*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".</em>


## List of instrinsics in AlgEtQ/FactPrimes.m:

> <pre><b>Factorization</b>(I::AlgEtQIdl) -> Tup</pre>
<em>Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.</em>

> <pre><b>PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQOrdIdl]</pre>
<em>Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.</em>

> <pre><b>PlacesAboveRationalPrime</b>(E::AlgEtQ,p::RngIntElt)->SeqEnum[AlgEtQIdl]</pre>
<em>Returns the maximal ideals of maximal order of the algebra E above the rational prime p.</em>

> <pre><b>SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgEtQOrdIdl]</pre>
<em>Returns the non-invertible primes of the order.</em>

> <pre><b>NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx</pre>
<em>Returns the non-invertible primes of the order.</em>

> <pre><b>IsPrime</b>(I::AlgEtQIdl) -> BoolElt</pre>
<em>Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.</em>

> <pre><b>IsMaximal</b>(I::AlgEtQIdl) -> BoolElt</pre>
<em>Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.</em>

> <pre><b>IsMaximalIdeal</b>(I::AlgEtQIdl) -> BoolElt</pre>
<em>Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.</em>

> <pre><b>IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
<em>Check if the order is Bass at the prime ideal P, that is, if every overorder of S is Gorenstein at the primes above P.</em>

> <pre><b>IsBass</b>(S::AlgEtQOrd) -> BoolElt</pre>
<em>Check if the order is Bass, that is, if every overorder of S is Gorenstein.</em>

> <pre><b>IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
<em>Check if the order is Gorenstein at the prime ideal P, that is, if every fractional ideal I with (I:I)=S is locally principal at P.</em>

> <pre><b>IsGorenstein</b>(O::AlgEtQOrd)->BoolElt</pre>
<em>Checks if the order O is Gorenstein, that is if the TraceDualIdeal of O is invertible, or equivalently, if all fractional ideals I with (I:I)=O are invertible.</em>


## List of instrinsics in AlgEtQ/PrimesAttributes.m:

> <pre><b>Valuation</b>(x::AlgEtQElt,P::AlgEtQIdl)->RngIntElt</pre>
<em>Given an element x and a maximal ideal P of the maximal order, it returns the valuation of x at P.</em>

> <pre><b>Valuation</b>(I::AlgEtQIdl,P::AlgEtQIdl)->RngIntElt</pre>
<em>Given a fractional ideal I and a maximal ideal P, both of the maximal order, it returns the valuation of I at P.</em>

> <pre><b>InertiaDegree</b>(P::AlgEtQIdl)->RngIntElt</pre>
<em>Given a maximal ideal P of the maximal order O above the rational prime p, it returns the inertia degree of P, that is, the index of the finite field extension GF(p)->O/P.</em>

> <pre><b>RamificationIndex</b>(P::AlgEtQIdl)->RngIntElt</pre>
<em>Given a maximal ideal P of the maximal order O, it returns the reamification index of P.</em>

> <pre><b>Uniformizers</b>(PPs::SeqEnum[AlgEtQIdl])->SeqEnum</pre>
<em>Given a sequence of maximal ideals P of the maximal order, it returns a sequence of elements t_P such that t_P is a uniformizer of P and a unit at every other prime in the sequence.</em>


## List of instrinsics in AlgEtQ/LowCohenMacaulayType.m:

> <pre><b>NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum</pre>
<em>Given an order S it returns two sequences: the first containis the primes at which S is locally not Gorenstein; the second containis the CohenMacaulay types of S at these primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.</em>

> <pre><b>CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt</pre>
<em>Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P\*S^t where S^t is the trace dual of S.</em>

> <pre><b>CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt</pre>
<em>Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P\*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.</em>


## List of instrinsics in AlgEtQ/WkClasses.m:

> <pre><b>WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]</pre>
<em>Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.</em>

> <pre><b>WeakEquivalenceClassesWithPrescribedMultiplicatorRing</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]</pre>
<em>Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.</em>

> <pre><b>WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum[AlgEtQIdl]</pre>
<em>Computes the weak equivalence class monoid of E. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.</em>


## List of instrinsics in AlgEtQ/WkTesting.m:

> <pre><b>IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt</pre>
<em>Checks if I and J are weakly equivalent, that is, if 1 \in (I:J)\*(J:I), or equivalently, if I and J are locally equivalent at all prime of their common multiplicator ring. This function does not require that the ideals are defined over the same order.</em>

> <pre><b>IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
<em>Check if the two orders are weakly equivalent, that is equal.</em>

> <pre><b>IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt</pre>
<em>Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.</em>

> <pre><b>IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt</pre>
<em>Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.</em>


## List of instrinsics in AlgEtQ/IdealClassMonoid.m:

> <pre><b>ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
<em>Returns the ideal classes of fractional S-ideals having MultiplicatorRing equal to S. This is the same as the orbit of the action of PicardGroup(S) on WKICM_bar(S).</em>

> <pre><b>ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
<em>Returns the ideal class monoid of the order, that is, a set of representatives for the isomorphism classes of the fractional S-ideals.</em>


## List of instrinsics in AlgEtQ/TotRealTotPos.m:

> <pre><b>IsTotallyReal</b>(a::AlgEtQElt) -> BoolElt</pre>
<em>Returns whether a is totally real.</em>

> <pre><b>IsTotallyRealPositive</b>(a::AlgEtQElt) -> BoolElt</pre>
<em>Returns whether a is totally positive, that is, totally real and with positive image in CC.</em>

> <pre><b>TotallyRealSubAlgebra</b>(K::AlgEtQ) -> AlgEtQ,Map</pre>
<em>Given a CM algebra K, returns the unique totally real subalgebra, with an embedding.</em>

> <pre><b>TotallyRealUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre>
<em>Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^\*.</em>

> <pre><b>TotallyRealPositiveUnitGroup</b>(S::AlgEtQOrd) -> Grp</pre>
<em>Given an order S in a CM étale algebra. Returns the groups of totally positive units of S, as a subgroup of S^\*.</em>


## List of instrinsics in AlgEtQ/PrintSave.m:

> <pre><b>PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt</pre>
<em>Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it returns a string that can be printed to file.</em>

> <pre><b>PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt</pre>
<em>Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered from this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.</em>

> <pre><b>LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd</pre>
<em>Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. These can be recovered with the approriate intrinsics.</em>


