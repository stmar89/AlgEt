# Étale algebras
 An `étale algebra` $A$ over a field $K$ is a finite product of finite separable extensions $K_1,\ldots,K_n$ of $K$.
 Typical examples are:
 - $A = K\times K$ where $K$ is a number field.
 - $A = \dfrac{\mathbb{K}[x]}{f(x)}$ where $f(x)$ is a polynomial with in $K[x]$ and no repeated roots over $\overline{K}$.
 We will refer to the field $K$ as the `prime field` of $A$ and to the fields $K_1,\ldots,K_n$  as the `components` of $A$. 
 If $F$ is a finite extension of $K$ such that $K_1,\ldots,K_n$ are all defined as relative extensions of $F$, we call $F$ the `base field` of $A$.
 If the components $K_1,\ldots,K_n$ of $A$ have distinct defining polynomials, say $f_1(x),\ldots,f_n(x) \in F[x]$ then $A$ is isomorphic to the étale algebra $F[x]/(f(x)$ where $f(x) = f_1(x)\cdot \cdots \cdot f_n(x)$. The polynomial $f(x)$ is then referred to as the `defining polynomial` of $A$.
# Étale algebras over $\mathbb{Q}$
 Currently we have implemented in MAGMA the type `AlgEtQ` which corresponds to étale algebras over the rational field $\mathbb{Q}$. 
 An étale algebra over $\mathbb{Q}$ of type `AlgEtQ` can be created using either a sequence of number fields, given as absolute extensions of $\mathbb{Q}$, or a polynomial in $\mathbb{Z}[x]$ or $\mathbb{Q}[x]$ with no repeated complex roots. 
 For such an algebra, the base field coincides with the prime field $\mathbb{Q}$.
<pre><b> EtaleAlgebra</b>(seq::SeqEnum[FldNum]) -> AlgEtQ</pre>
Given a sequence of number fields which are absolute extensions of the rational field, returns the étale algebra corresponding to the direct product. The number fields with DefiningPolynomial of degree one need to created with DoLinearExtension set to true.
<pre><b> EtaleAlgebra</b>(f::RngUPolElt[RngInt]) -> AlgEtQ</pre><pre><b> EtaleAlgebra</b>(f::RngUPolElt[FldRat]) -> AlgEtQ</pre>
 Given a polynomial with integer of rational coefficients, which is squarefree, that is, with no repeated complex roots, returns the étale algebra which is the product of the number fields defined by the irreducible factors of the polynomial.
## Direct products
<pre><b> DirectProduct</b>(seq::SeqEnum[AlgEtQ]) -> AlgEtQ,SeqEnum[Map],SeqEnum[Map]</pre>
 Given a sequence of étale algebras over $\mathbb{Q}$, returns the direct product, together with sequences of inclusions and projections.
## Conversion to a number field
<pre><b> IsNumberField</b>(A::AlgEtQ) -> BoolElt,FldNum,Map </pre>
 Given an étale algebra $A$ over $\mathbb{Q}$, returns wheter $A$ is a number field, that is, has only one component.  If this is the case, then returns also the number field itself together an isomorphism from $A$ to the number field.
## Attributes
### Components, equality testing and defining polynomial
 Two étale algebras are defined to be equal if the ordered sequence of their components are the same.
<pre><b> Components</b>(A::AlgEtQ) -> SeqEnum,SeqEnum,SeqEnum</pre>
 Returns the number fields of which $A$ is a product of, together with embeddings and projections.
<pre><b> 'eq'</b>(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt</pre>
Two étale algebras are defined to be equal if they have the same components.
<pre><b> DefiningPolynomial</b>(A::AlgEtQ) -> RngUPolElt</pre>
 Returns the defining polynomial of $A$, if the components are distinct number fields.
### Base field and prime field
<pre><b> HasBaseField</b>(A::AlgEtQ) -> BoolElt</pre>
 Returns whether the components of $A$ all have the same base field.
<pre><b> BaseField</b>(A::AlgEtQ) -> FldNum,Map</pre>
 Returns whether the common base field of the components of $A$ if it exists.
<pre><b> PrimeField</b>(A::AlgEtQ) -> FldNum</pre>
Returns the prime field of the étale algebra.
### Dimension
<pre><b> Dimension</b>(A::AlgEtQ)->RngInt</pre>
 Returns the dimension of A over the base field, which in this case is $\mathbb{Q}$.
<pre><b> AbsoluteDimension</b>(A::AlgEtQ)->RngInt</pre>
 Returns the dimension of $A$ over the prime field.
# Elements of étale algebras over $\mathbb{Q}$
 An element $x$ of an étale algebra $A$ over $\mathbb{Q}$ with components $K_1,\ldots,K_n$ is stored as a tuple of elements of the number fields. Such a tuple is referred to as the `components` of the element.
 Note that $x$ is a `unit` of $A$, that is, an invertible element, if all its components are non-zero. Otherwise is a `zero-divisor` of $A$.
 Elements of étale algebra of type `AlgEtQ` have type `AlgEtQElt`.
## Creation 
<pre><b> '!'</b>(A::AlgEtQ, x::.) -> AlgEtQElt</pre>
 Let $A$ be an étale algebra over $\mathbb{Q}$ with components $K_1,\ldots,K_n$.
 Elements of $A$ are created using the operator `!` by giving either
 - a sequence, a tuple, or a list containing $n$ elements $x_1,\ldots,x_n$ with each $x_i$ coercible in $K_i$, or
 - an integer or a rational number.
## Attributes
<pre><b> Parent</b>(x::AlgEtQElt) -> AlgEtQ</pre><pre><b> Algebra</b>(x::AlgEtQElt) -> AlgEtQ</pre>
Returns the algebra to which the element belongs to.
<pre><b> Components</b>(x::AlgEtQElt) -> Tup</pre>
Returns the components of the element.
## Basis and coordinates
<pre><b> Basis</b>(A::AlgEtQ) -> SeqEnum</pre>
Returns a basis of the algebra over the base field.
<pre><b> AbsoluteBasis</b>(A::AlgEtQ) -> SeqEnum</pre>
Returns a basis of the algebra over the prime field.
<pre><b> '.'</b>(A::AlgEtQ,i::RngIntElt)->AlgEtQElt</pre>
 Returns the $i$-th element of the absolute basis.
<pre><b> AbsoluteCoordinates</b>(x::AlgEtQElt) -> SeqEnum</pre>
Given an element, returns the coordinates relative to the absolute basis, which are elements of the prime rational field.
<pre><b> AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum</pre>
Given a sequence of elements and a basis over the prime field, returns a sequence whose entries are the coordinates with respect to the given basis.
## Special elements
<pre><b> One</b>(A::AlgEtQ) -> AlgEtQElt</pre>
The multiplicative neutral element.
<pre><b> Zero</b>(A::AlgEtQ) -> AlgEtQElt</pre>
The additive neutral element.
<pre><b> Idempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
Returns the idempotents of the étale algebra.
<pre><b> OrthogonalIdempotents</b>(A::AlgEtQ) -> SeqEnum</pre>
Returns the orthogonal idempotent elements of the étale algebra, which are the unit elements of the components.
## Primitive element
 Given an étale algebra $A$ over $\mathbb{Q}$ there exists an element $a\in A$ such that $A = \mathbb{Q}[a]$, that is, every element can be written as a polynomial with rational coefficients in $a$. Such an element is called a `primitive element` of $A$. It is characterized by having a minimal polynomial whose degree equals the absolute dimension of $A$.
 The intrinsic `SetPrimitiveElement` allows to specify a primitive element, which is then stored in an attribute of the algebra.
 The intrinsic `PrimitiveElement` produces such an element of the étale algebra $A$ using a deterministic procedure which we now describe:
 Let $N$ be the number of components of $A$, each one having primitive element $a_i$. Set $b_1$ = $a_1$. For $i=2,\ldots,N$, set $b_i = a_i+j$ where $j$ is the smallest non-negative integer such that the minimal polynomial of $a_i+j$ is not in the set of minimal polynomials of the elements $b_1,\ldots,b_{i-1}$. The output is the element of $A$ whose components are $b_1,...,b_N$. In particular, if $A$ is a product of number fields with different defining polynomials, then the output is the element of $A$ whose components are the primitive elements of the components.
<pre><b> PrimitiveElement</b>(A::AlgEtQ) -> AlgEtQElt</pre>
 Returns a primitive element of the étale algebra, produced with a deterministic algorithm. If the components of algebra are distinct number fields, then the output is the element whose components are the primitive elements of the components.
<pre><b> SetPrimitiveElement</b>(a::AlgEtQElt)</pre>
Given an element a in an étale algebra A such that the minimal polynomial of a over Q has degree equal to the Q-dimension of A, it sets a to be the primitive element of A.
<pre><b> PowerBasis</b>(A::AlgEtQ) -> SeqEnum[AlgEtQElt]</pre>
 Returns the basis consisting of powers of the element stored as the primitive element.
## Units and zero divisors
<pre><b> IsUnit</b>(x::AlgEtQElt) -> BoolElt</pre>
Returns whether the element is a unit.
<pre><b> IsZeroDivisor</b>(x::AlgEtQElt) -> BoolElt</pre>
Returns whether the element is a zero-divisor.
## Random elements
<pre><b> Random</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre><pre><b> Random</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt</pre>
Returns a random element of the étale algebra with coefficients are bounded by the positive integer bd.
<pre><b> RandomUnit</b>(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt</pre><pre><b> RandomUnit</b>(A::AlgEtQ : bd:=3) -> AlgEtQElt</pre>
Returns a random unit of the étale algebra with coefficients are bounded by the positive integer bd.
## Equality and operations on elements
 Equality testing and operations are performed component-wise, that is, at the level of the components.
<pre><b> 'eq'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt</pre><pre><b> 'eq'</b>(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt</pre><pre><b> 'eq'</b>(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt</pre><pre><b> 'eq'</b>(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt</pre><pre><b> 'eq'</b>(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt</pre>
 Returns whether the two elements are equal.
<pre><b> '+'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre><pre><b> '+'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
 Returns the sum of the two elements.
<pre><b> '-'</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
 Returns the opposite element.
<pre><b> '-'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre><pre><b> '-'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
 Returns the difference of two elements.
<pre><b> '*'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre><pre><b> '*'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
 Returns the product of two elements.
<pre><b> Inverse</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
 Returns the multiplicative inverse of an unit element of the étale algebra.
<pre><b> '^'</b>(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt</pre>
 Given an element $x$ and an integer $n$ returns $x^n$.
<pre><b> '/'</b>(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::.,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::AlgEtQElt,x2::.) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt</pre><pre><b> '/'</b>(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt</pre>
 Returns the division of an element $x_1$ by the unit $x_2$.
## Operations on sequences of elements.
<pre><b> '&+'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
Returns the sum of the elements of the sequence.
<pre><b> '&*'</b>(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre>
Returns the product of the elements of the sequence.
<pre><b> SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre><pre><b> SumOfProducts</b>(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre><pre><b> SumOfProducts</b>(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt</pre><pre><b> SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt</pre><pre><b> SumOfProducts</b>(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt</pre>
 Given two sequences of $n$ elements $[a_1,\ldots,a_n]$ and $[b_1,\ldots,b_n]$ returns $\sum_{i=1}^n a_i\cdot b_i$.
## Minimal polynomials and integrality testing.
 Given an element $a$ of an étale algebra over $\mathbb{Q}$, the `minimal polynomial` is the polynomial $f(x)$ in $\mathbb{Q}[x]$ of minimal degree such that $f(a)=0$. It is the least common multiple of the minimal polynomials of the components of $a$.
 An element $a$ is called `integral` if its minimal polynomial is monic and has integer coefficients.
<pre><b> MinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
Returns the minimal polynomial of the element over the base field.
<pre><b> AbsoluteMinimalPolynomial</b>(x::AlgEtQElt) -> RngUPolElt</pre>
Returns the minimal polynomial of the element over the prime field.
<pre><b> IsIntegral</b>(x::AlgEtQElt) -> BoolElt</pre>
Returns whether the element is integral (over the integers).
<pre><b> Evaluate</b>(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt</pre>
 Evaluate the polynomial $f(x)$ at the element $a$.
# Example 1
 ```
 _<x>:=PolynomialRing(Integers());
 f:=(x^8+16)*(x^8+81);
 A:=EtaleAlgebra(f);
 // We compute the `canonical` primitive element, which is the class of the variable x in A.
 a:=PrimitiveElement(A); a;
 
 // The algebra A has two components:
 comps,embeddings,projections:=Components(A);
 K1,K2:=Explode(comps);
 // The unit elemenet of each component corresponds to an orthogonal idempotent of A:
 [ embeddings[1](K1!1),embeddings[2](K2!1) ] eq OrthogonalIdempotents(A);
 
 // We conclude this example by showing the use of SumOfProducts and its timings advantages:
 N:=10^5;
 elts1:=[ a+i : i in [1..N] ];
 elts2:=[ a-i : i in [1..N] ];
 time s1:=&+[ elts1[i]*elts2[i] : i in [1..N] ];
 time s2:=SumOfProducts(elts1,elts2);
 s1 eq s2;
 ```
# Example 2
 ```
 //We now consider the étale algebra consisting of two copies of the rational field.
 _<x>:=PolynomialRing(Integers());
 QQ:=NumberField(x-1:DoLinearExtension);
 A:=EtaleAlgebra([QQ,QQ]);
 a:=PrimitiveElement(A); a;
 b:=A!<2,10>;
 SetPrimitiveElement(b);
 ```
# Homomorphisms of étale algebras over $\mathbb{Q}$
 Let $A$ be an étale algebra over $\mathbb{Q}$. By an `homomorphism` from $A$ to some $\mathbb{Q}$-algbra, we mean a unital $\mathbb{Q}$-algbra homomorphisms. 
## Homomorphisms to the complex numbers
 The set of homomorphisms from an étale algebra $A$ to the field of complex numbers consists of the homomorphisms acting as an embedding on a single compoenent and zero on every other component. Such a homomorphism is injective if and only if $A$ has a unique component, that is, is a number field.
<pre><b> HomsToC</b>(A::AlgEtQ : Prec:=Precision(GetDefaultRealField()))->SeqEnum[Map]</pre>
Returns the sequence of homomorphisms from the étale algebra to the complex field. The precision of the target can be set by the vararg "Prec".
## Homomorphisms between étale algebras over $\mathbb{Q}$
<pre><b> Hom</b>(A::AlgEtQ , B::AlgEtQ , img::SeqEnum[AlgEtQElt] : CheckMultiplicative:=false, CheckUnital:=false, ComputeInverse:=true)->Map</pre>
 Given two étale algebras $A$ and $B$ and a sequence of elements $img$ of $B$, returns the homomorphism defined by sending the AbsoluteBasis of A to $img$. The parameter CheckMultiplicative (default false) determines if the multiplicativity of the defined map is checked, while the parameter CheckUnital (default false) determines whether it is unital. If the parameter ComputeInverse (default true) is true, it checkes whether the map is invertible and, if so, it defines also the inverse (by assigning preimages).}
<pre><b> DiagonalEmbedding</b>(K::AlgEtQ, V::AlgEtQ)->Map</pre><pre><b> NaturalAction</b>(K::AlgEtQ, V::AlgEtQ)->Map</pre>
 Given an étale algebra $K$ of the form $K_1\times \cdots \times K_n$ and an ´étale algebra $V$ of the form $K_1^{s_1} \times \cdots \times K_n^{s_n}$, returns the natural componentwise diagonal embedding $K\to V$. 
# Orders in étale algebras over $\mathbb{Q}$
 Let $A$ be an étale algebra over $\mathbb{Q}$. An `order` in $A$ is a subring $R$ of $A$ whose underlying additive group is a free abelian group of rank equal to the absolute dimension of $A$.
 Given an order $R$ in $A$, we say that $S$ is an `overorder` of $R$ if $R\subseteq S$.
 In MAGMA, orders have type `AlgEtQOrd`. Elements will always be considered as elements of the algebra.
## Creation
 Whenever we create an order, we populate the attributes `Generators` and `ZBasis`.
 Unless the paramater `Check` is set to $0$, the `ZBasis` is put in in a canonical form (we use the hermite normal form of the numerators and divide by the denominators). The entries of the corresponding upper triangular matrix uniquely determine the order and are used for hashing it. 
<pre><b> Order</b>( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd</pre>
 Returns the order generated by a given sequence of elements. The parameter `Check` (default $100$) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the given generators. If `Check` is $0$ then this step is skipped. The parameter `CheckIsKnownOrder` determines whether the program checks if the order is already known, i.e. in the attribute `KnownOrders` of the algebra. The default value is true.
<pre><b> Order</b>(A::AlgEtQ , orders::Tup) -> AlgEtQOrd</pre>
Given an étale algebra and a tuple of orders in the components of the étale algebra, it generates the direct sum order.
## Equation orders
 If $a$ is the stored primitive element of $A$ then we refer to the order $\mathbb{Z}[a]$ as the `equation order` of $A$.
<pre><b> EquationOrder</b>(A::AlgEtQ) -> AlgEtQOrd</pre>
Returns the monogenic orded of the étale algebra, which depends on the stored primitive element.
<pre><b> ProductOfEquationOrders</b>(A::AlgEtQ)->AlgEtQOrd</pre>
Given an étale algebra A, returns the order consisting of the product of the equation orders of the number fields.
## Maximal Order
 If $A$ has components $K_1,\ldots,K_n$ then product $\mathcal{O}_A$ of the ring of integers $\mathcal{O}_{K_i}$ of $K_i$ is an order, which contains every other order of $A$. For this reason, we refer to $\mathcal{O}_A$ as the `maximal order` of $A$.
<pre><b> MaximalOrder</b>(A::AlgEtQ)->AlgEtQOrd</pre>
Returns the maximal order of the étale algebra A. It is the direct sum of the ring of integers of the number fields composing the algebra.
<pre><b> IsMaximal</b>(S::AlgEtQOrd) -> BoolElt</pre>
Returns wheter the given order is the maximal order of the étale algebra.
## The attribute `KnownOrders`
 An étale algebra $A$ has an attribute `KnownOrders`, in which all previously constructed orders are stored.
 Whenever the creation intrinsc for orders is called, the program checks if the orders has already been costructed. This is done to prevent the creation of the same order multiple times. This check can be turned off with the approriate parameter.
<pre><b> IsKnownOrder</b>(~R::AlgEtQOrd)</pre>
This procedure checks whether the given order is already in the list of known orders of the parent algebra. If so then it replaces the given order with the stored copy. If not, it adds is added to the storage.
## Main attributes of an order
<pre><b> Algebra</b>(S::AlgEtQOrd) -> AlgEtQ</pre>
Returns the algebra of the order.
<pre><b> ZBasis</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
 Returns the stored $\mathbb{Z}$-basis of the order.
<pre><b> Generators</b>(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]</pre>
 Return a set of generators as a \mathbb{Z}}-algebra of the order.
## Equality
 Equality of orders is perfomed using the `Hash` attribute which is constructed as follows.
 Let $S$ be an order. 
 Let $P$ be the upper triangular Hermite normal form of the integer square matrix $d\cdot M$ where $M$ is the matrix whose rows are the coefficients of a $\mathbb{Z}$-basis of $S$ and $d$ is the least common denominator of its entries.
 The `Hash` of $S$ is defined to be the sequence consisting of the least common denominator of $\frac{1}{d}\cdot P$ and the entries of the upper triangular part of $\frac{1}{d}\cdot P$. 
 This hashing method has no collisions and it is independent of the choice of $\mathbb{Z}-basis$ from which we start the procedure.
 We observed that applying the inbuild Hash function to the sequence defined above, while giving a smaller hash, it often lead to collisions.
<pre><b> myHash</b>(S::AlgEtQOrd)->SeqEnum[RngInt]</pre>
Hash function for AlgEtQOrd.
<pre><b> 'eq'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
Checks equality of orders.
## Inclusion and coordinates
 Inclusion testing of elements, orders and ideals in a fixed order is perfomed by multiplying by the `inclusion matrix`. This matrix, which is stored in an attribute, is the inverse of the matrix with coefficients the $\mathbb{Z}$-basis of the order. If the output of the multiplication has integer coefficients then we have an inclusion.
 The same matrix can be used to obtain the coordinates of a sequence of elements with respect to the $\mathbb{Z}$-basis.
<pre><b> 'in'</b>(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt</pre>
Inclusion of an element in the order.
<pre><b> 'in'</b>(x::RngIntElt,O::AlgEtQOrd) -> BoolElt</pre>
Inclusion of elements.
<pre><b> 'in'</b>(x::FldRatElt,O::AlgEtQOrd) -> BoolElt</pre>
Inclusion of elements.
<pre><b> 'subset'</b>(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt</pre>
Checks if the first argument is inside the second.
<pre><b> AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum</pre>
 Returns the absolute coordinates of a sequence of elements  with respect to the stored $\mathbb{Z}$-basis of the given order.
## Special and random elements in an order
<pre><b> One</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
Returns the multiplicative unit of the order.
<pre><b> Zero</b>(S::AlgEtQOrd)->AlgEtQElt</pre>
Returns the additive neutral element of the order.
<pre><b> Random</b>(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
 Returns a random element of the given order with coefficients with respect to the stored $\mathbb{Z}$-basis bounded by `bd`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
<pre><b> Random</b>(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
 Returns a random element of the given order. The coefficients with respect to the stored $\mathbb{Z}$-basis can be bounded by setting the parameter `CoeffRange`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
<pre><b> Index</b>(T::AlgEtQOrd) -> FldRatElt</pre>
 Given an order $T$ returns the determinant of the change of basis from a $\mathbb{Z}-basis of $T$ to the basis of the parent algebra.
<pre><b> Index</b>(S::AlgEtQOrd, T::AlgEtQOrd) -> Any</pre>
 Given two orders $S$ and $T$ returns the ratio of their indices. If $T subset S$ the outputs is the index of the inclusion, that is, it equals $[S:T] = \#S/T$.}
## Product and intersection of orders
<pre><b> '*'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
Returns the order generated by the orders O1 and O2.
<pre><b> 'meet'</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd</pre>
Returns the order given by the intersection of two given orders.
## Orders which are direct products of orders
<pre><b> IsProductOfOrdersInComponents</b>(O::AlgEtQOrd)->BoolElt, Tup</pre><pre><b> IsProductOfOrders</b>(O::AlgEtQOrd)->BoolElt, Tup</pre>
Returns whether the argument is a product of orders in the components of its parent algebra. If so, it returns also a tuple containing these orders.
<pre><b> IsProductOfOrdersInFactorAlgebras</b>(S::AlgEtQOrd)->BoolElt,SeqEnum[AlgEtQElt],SeqEnum</pre>
Returns whether the given order is a product of orders living in some factor algebras of the parent algebra. 
This is equivalent for the given order to contain some idempotents of teh algebra other than 0 and 1. If this is the case, it returns also the idempotents.
# Example 3
 ```
 _<x>:=PolynomialRing(Integers());
 // We consider the following three number fields
 K1:=NumberField(x^2-2);
 K2:=NumberField(x^2-3);
 K3:=NumberField(x^2-5);
 // We define the product étale algebra A and the factor algebras B and C consisting of only the first two components and the last one, repsectively.
 B:=EtaleAlgebra([K1,K2]);
 C:=EtaleAlgebra([K3]);
 A,embs,projs:=DirectProduct([B,C]);
 // The maximal order of A is the product of the three ring of integers.
 OA:=MaximalOrder(A);
 IsProductOfOrdersInComponents(OA);
 // The equation order of A is not a product in any factor algebra
 EA:=EquationOrder(A);
 IsProductOfOrdersInFactorAlgebras(EA);
 // Now we construct an order that is a product of an order in B and one in C, but does not admit further splittings.
 a:=PrimitiveElement(A);
 e1:=A![1,0,0];
 e2:=A![0,1,1];
 R:=Order([a*e1,a*e2]);
 IsProductOfOrdersInFactorAlgebras(R);
 ```
# Ideals of orders in étale algebras over $\mathbb{Q}$
 Let $R$ be an order in an étale algebra $A$ over $\mathbb{Q}$. 
 A `fractional ideal` over $R$, or fractional $R$-ideal, is a sub-$R$-module of $A$ whose underlying additive group is free of rank equal to the absolute dimension of $A$.
 A fractional $R$-ideal is said to be `integral` if $I\subseteq R$.
 Note that every overorder of $R$ is a fractional $R$-ideal.
 
 In MAGMA, fractional ideals have type `AlgEtQIdl`. Elements will always be considered as elements of the algebra.
## Creation
<pre><b> Ideal</b>(S::AlgEtQOrd, gens::SeqEnum) -> AlgEtQIdl</pre>
 Given an order $S$ in an étale algebra and a sequence `gens` of elements coercible in the algebra, returns the fractional $S$-ideal generated by `gens`.}
<pre><b> Ideal</b>(S::AlgEtQOrd, idls::Tup) -> AlgEtQIdl</pre>
 Given an order $S$ which is a product of orders $S_i$ in the components of the parent étale algebra, and a tuple of fractional $S_i$-ideals $I_i$, returns the $S$-ideal corresponding to the direct sum of the $I_i$.}
<pre><b> Ideal</b>(S::AlgEtQOrd, gen::Any) -> AlgEtQIdl</pre><pre><b> '*'</b>(S::AlgEtQOrd, gen::AlgEtQElt) -> AlgEtQIdl</pre><pre><b> '*'</b>(gen::AlgEtQElt, S::AlgEtQOrd) -> AlgEtQIdl</pre><pre><b> '*'</b>(S::AlgEtQOrd, gen::RngIntElt) -> AlgEtQIdl</pre><pre><b> '*'</b>(gen::RngIntElt, S::AlgEtQOrd) -> AlgEtQIdl</pre><pre><b> '*'</b>(S::AlgEtQOrd, gen::FldRatElt) -> AlgEtQIdl</pre><pre><b> '*'</b>(gen::FldRatElt, S::AlgEtQOrd) -> AlgEtQIdl</pre>
 Returns the fractional $S$-ideal generated by `gen`.
## Coercion of fractional ideals
 Let $S$ and $T$ be orders in an étale algebra and let $I$ be a fractional $S$-ideal.
 If $T \subseteq S$ then $I$ is naturally also a fractional $T$-ideal.
 Otherwise, one can consider the fractional $T$-ideal $IT$ generated by $I$.
 This is done using the operator `!!`.
<pre><b> '!!'</b>(T::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
 Given an order $T$ and a fractional $S$-ideal $I$ for some order $S$ in the same étale algebra of $T$, returns the fractional $T$-ideal $IT$. If $T$ is a subset of $S$, then $IT$ and $I$ coincides as sets.
## Attributes
<pre><b> Algebra</b>(I::AlgEtQIdl) -> AlgEtQ</pre>
Returns the étale algebra of the fractional ideal.
<pre><b> Order</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
Returns the order of definition of the fractional ideal.
<pre><b> ZBasis</b>(I::AlgEtQIdl)->SeqEnum[AlgEtQElt]</pre>
 Returns a $\mathbb{Z}$-basis of the fractional ideal.
<pre><b> Generators</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre>
Returns a set of generators of the fractional ideal.
## Hashing and equality
 Let $I$ a fractional $S$-ideal in an étale algebra over $\mathbb{Q}$.
 Let $P$ be the upper triangular Hermite normal form of the integer square matrix $d\cdot M$ where $M$ is the matrix whose rows are the coefficients of a $\mathbb{Z}$-basis of $I$ and $d$ is the least common denominator of its entries.
 The `Hash` of $I$ is defined to be the sequence consisting of the least common denominator of $\frac{1}{d}\cdot P$ and the entries of the upper triangular part of $\frac{1}{d}\cdot P$. 
 This hashing method has no collisions and it is independent of the choice of $\mathbb{Z}-basis$ from which we start the procedure.
 We observed that applying the inbuild Hash function to the sequence defined above, while giving a smaller hash, it often lead to collisions.
<pre><b> myHash</b>(I::AlgEtQIdl)->RngInt</pre>
Hash function.
<pre><b> 'eq'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre><pre><b> 'eq'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt</pre><pre><b> 'eq'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
Equality testing.
<pre><b> 'ne'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> BoolElt</pre><pre><b> 'ne'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> BoolElt</pre><pre><b> 'ne'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre>
Disequality testing.
## Coordinates with respect to the $\mathbb{Z}$-basis and inclusion
 Inclusion testing of elements, orders and ideals in a given fractional ideal is perfomed by multiplying by the `inclusion matrix`. This matrix, which is stored in an attribute, is the inverse of the matrix with coefficients the $\mathbb{Z}$-basis of the order. If the output of the multiplication has integer coefficients then we have an inclusion.
 The same matrix can be used to obtain the coordinates of a sequence of elements with respect to the $\mathbb{Z}$-basis.
<pre><b> AbsoluteCoordinates</b>(seq::SeqEnum[AlgEtQElt],I::AlgEtQIdl) -> SeqEnum</pre>
 Given a sequence of elements and a fractional ideal $I$, returns the sequence of coordinates of the elements with respect to the stored $\mathbb{Z}$-basis of $I$.
<pre><b> 'in'</b>(x::AlgEtQElt , I::AlgEtQIdl ) -> BoolElt</pre><pre><b> 'in'</b>(x::RngIntElt , I::AlgEtQIdl ) -> BoolElt</pre><pre><b> 'in'</b>(x::FldRatElt , I::AlgEtQIdl ) -> BoolElt</pre>
Returns whether the element is in the given fractional ideal.
<pre><b> 'subset'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> BoolElt</pre><pre><b> 'subset'</b>(I::AlgEtQIdl,S::AlgEtQOrd) -> BoolElt</pre><pre><b> 'subset'</b>(I1 :: AlgEtQIdl, I2 :: AlgEtQIdl) -> BoolElt</pre>
Inclusion testing.
## Index
<pre><b> Index</b>(I::AlgEtQIdl) -> FldRatElt</pre>
 Given a fractional ideal $I$ returns the determinant of the change of basis from a $\mathbb{Z}$-basis of $I$ to the basis of the parent algebra.
<pre><b> Index</b>(J::AlgEtQIdl, I::AlgEtQIdl) -> Any</pre>
 Given fractional ideals $J$ and $I$ defined over the same order returns $[J:I] = [J:J \cap I]/[I : J \cap I]$.
<pre><b> Index</b>(S::AlgEtQOrd, I::AlgEtQIdl) -> Any</pre>
 Given an order $S$ and a fractional $S$-ideal $I$ returns $[S:I] = [S:S \cap I]/[I : S \cap I]$.
## Special ideals
<pre><b> OneIdeal</b>(S::AlgEtQOrd) -> AlgEtQIdl</pre>
Returns the ideal generated by the multiplicative unit of the given order.
<pre><b> Conductor</b>(O::AlgEtQOrd) ->AlgEtQOrdIdl</pre>
 Computes the conductor of an order $R$, defined as the colon ideal $(R:\mathcal{O})$, where $\mathcal{O}$ is the maximal order of the algebra.
## Operations on ideals
 The set of fractional $R$-ideals is closed under addition, multiplication, intersection and colon ideal.
 More precisely, given two fractional $R$-ideals $I$ and $J$ in the étale algebra $A$, we define
 - $I+J = \{ a + b : a \in I , b\in J \}$,
 - $I*J = \{ a*b : a \in I, b \in J \}$,
 - $I \cap J = \{ c \in A : c \in I \text{ and } c \in J \}$,
 - $(I : J) = \{ c \in A : c\cdot J \subseteq I\}$.
<pre><b> '+'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
Returns the sum of the two fractional ideals.
<pre><b> '*'</b>(I::AlgEtQIdl , J::AlgEtQIdl ) -> AlgEtQIdl</pre>
Returns the product of the two fractional ideals.
<pre><b> '*'</b>(I::AlgEtQIdl , x::AlgEtQElt ) -> AlgEtQIdl</pre><pre><b> '*'</b>(x::AlgEtQElt, I::AlgEtQIdl) -> AlgEtQIdl</pre><pre><b> '*'</b>(x::RngIntElt, I::AlgEtQIdl) -> AlgEtQIdl</pre><pre><b> '*'</b>(I::AlgEtQIdl, x::RngIntElt) -> AlgEtQIdl</pre><pre><b> '*'</b>(x::FldRatElt, I::AlgEtQIdl) -> AlgEtQIdl</pre><pre><b> '*'</b>(I::AlgEtQIdl, x::FldRatElt) -> AlgEtQIdl</pre>
 Returns the fractional ideal $x\cdot I$.
<pre><b> '^'</b>(I::AlgEtQIdl, n::RngIntElt) -> AlgEtQIdl</pre>
 Returns the $n$-th power of the fractional ideal.
<pre><b> 'meet'</b>(I::AlgEtQIdl, S::AlgEtQOrd) -> AlgEtQIdl</pre><pre><b> 'meet'</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> AlgEtQIdl</pre>
 Given a fractional ideal $I$ of the order $S$, returns $S \cap I$.
<pre><b> 'meet'</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> AlgEtQIdl</pre>
Given two fractional ideals over the same order, returns their intersection.
<pre><b> ColonIdeal</b>(I::AlgEtQIdl,J::AlgEtQIdl)->AlgEtQIdl</pre>
 Given fractional ideals $I$ and $J$ over the same order returns their colon ideal $(I:J)$.
<pre><b> ColonIdeal</b>(O::AlgEtQOrd,J::AlgEtQIdl)->AlgEtQIdl</pre>
 Given an order $\mathcal{O}$ and a fractional $\mathcal{O}$-ideal $J$, returns $(\mathcal{O}:J)$.
<pre><b> ColonIdeal</b>(I::AlgEtQIdl,O::AlgEtQOrd)->AlgEtQIdl</pre>
 Given an order $\mathcal{O}$ and a fractional $\mathcal{O}$-ideal $I$, returns $(I:\mathcal{O})$.
## Invertibility and multiplicator ring
 Let $I$ be a fractional $R$-ideal. 
 The `multiplicator ring` of $I$ is the biggest order $S$ such that $I$ is a fractional $S$-ideal.
 It is easy to verify that the multiplicator ring of $I$ is the overorder $(I:I)$ of $R$.
 We say that $I$ is `invertible` if $I\cdot (R:I) = R$.
 If $I$ is invertible, we say that $(R:I)$ is its `inverse`.
 This definition depends on the order of definition of $I$.
 In fact, if $I$ is invertible then $R=(I:I)$.
 If $R$ is the maximal order of the algebra, then every fractional $R$-ideal is invertible.
<pre><b> IsInvertible</b>(I::AlgEtQIdl) ->BoolElt</pre>
Checks if the ideal I is invertible in its order of definition.
<pre><b> Inverse</b>(I::AlgEtQIdl) ->AlgEtQIdl</pre>
Computes the inverse of a given invertible ideal.
<pre><b> MultiplicatorRing</b>(I::AlgEtQIdl) -> AlgEtQOrd</pre>
Computes the multiplicator ring of the given fractional ideal.
<pre><b> MultiplicatorRing</b>(R::AlgEtQOrd) -> AlgEtQOrd</pre>
Returns the multiplicator ring of the order, that is, the order itself.
<pre><b> Random</b>(I::AlgEtQIdl , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
 Returns a random element of the given fractonal ideal with coefficients with respect to the stored $\mathbb{Z}$-basis bounded by `bd`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
<pre><b> Random</b>(I::AlgEtQIdl : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt</pre>
 Returns a random element of the given fractional ideal. The coefficients with respect to the stored $\mathbb{Z}$-basis can be bounded by setting the parameter `CoeffRange`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
## Various properties of the ideals
<pre><b> IsProductOfIdeals</b>(I::AlgEtQIdl) -> BoolElt, Tup</pre>
 Returns if the argument is a product of ideals in the number fields defining the algebra. If so, it returns also the sequence of these ideals (in the appropriate orders). We require the order of definition of the ideal to be its multiplicator ring.
<pre><b> IsCoprime</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> BoolElt</pre>
 Given two integral ideals $I$ and $J$ over the same order $S$, returns whether $I+J=S$.}
<pre><b> IsIntegral</b>(I::AlgEtQIdl) -> BoolElt</pre>
Returns wheter the given fractional ideal is integral in its order of definition.
<pre><b> MakeIntegral</b>(I::AlgEtQIdl) -> AlgEtQIdl,RngIntElt</pre>
 Given a fractional $S$-ideal $I$, returns the ideal $d*I$, where $d$ is the smallest integer such that $d\cdot I$ is integral. Compare with the intrinsic `SmallRepresentative`.
<pre><b> MinimalInteger</b>(I::AlgEtQIdl) -> RngIntElt</pre>
Returns the smallest integer contained in the given fractional ideal, which must be integral.
<pre><b> CoprimeRepresentative</b>(I::AlgEtQIdl,J::AlgEtQIdl) -> AlgEtQElt,AlgEtQIdl</pre>
 Given an invertible fractional ideal $I$ and an integral ideal $J$, defined over the same order, returns a unit $x$ of the parent algebra such that $x\cdot I$ is an integral fractional ideal coprime with $J$, togheter with the product $x\cdot I$.
<pre><b> ZBasisLLL</b>(S::AlgEtQOrd)</pre>
A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.
<pre><b> ZBasisLLL</b>(S::AlgEtQIdl)</pre>
A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.
<pre><b> Quotient</b>(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
Given an ideal I and the ZBasis of an ideal or order J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.
<pre><b> Quotient</b>(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map</pre>
Given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J.
<pre><b> Quotient</b>(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map</pre>
Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order.
<pre><b> ResidueRing</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map</pre>
Given an integral ideal I of S, returns the abelian group S/I and the quotient map q:S -> S/I (with preimages). Important: the domain of q is the Algebra of S, since the elements of S are expressed as elements of A. We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.
<pre><b> ResidueField</b>(P::AlgEtQIdl) -> FldFin, Map</pre>
Given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.
<pre><b> PrimitiveElementResidueField</b>(P::AlgEtQIdl)->AlgEtQElt</pre>
Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^*.
<pre><b> QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).
<pre><b> QuotientVS</b>(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).
<pre><b> QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map</pre>
Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).
<pre><b> QuotientVS</b>(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map</pre>
Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with preimages).
<pre><b> IsMaximalAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt</pre>
Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.
<pre><b> MinimalOverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]</pre>
Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
<pre><b> MinimalOverOrders</b>(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]</pre>
Computes the minimal overorders of R.
<pre><b> OverOrdersAtPrime</b>(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]</pre>
Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
<pre><b> OverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]</pre>
We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.
<pre><b> FindOverOrders</b>(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]</pre>
We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.
<pre><b> GraphOverOrders</b>(R:AlgEtQOrd) -> GrphDir</pre>
Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].
<pre><b> Trace</b>(x::AlgEtQElt) -> Any</pre>
Returns the trace of the element x of an étale algebra.
<pre><b> Norm</b>(x::AlgEtQElt) -> Any</pre>
Returns the norm of the element x of an étale algebra.
<pre><b> AbsoluteTrace</b>(x::AlgEtQElt) -> Any</pre>
Returns the absolute trace of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Trace.
<pre><b> AbsoluteNorm</b>(x::AlgEtQElt) -> Any</pre>
Returns the absolute norm of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Norm.
<pre><b> TraceDualIdeal</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
Returns the trace dual ideal of the ideal I, that is, the set of elements x of the algebra such that Trace(x*I) is integer-valued.
<pre><b> TraceDualIdeal</b>(O::AlgEtQOrd) -> AlgEtQIdl</pre>
Returns the trace dual ideal of an order in an etale algebra, that is, the set of elements x of the algebra such that Trace(x*O) is integer-valued.
<pre><b> Completion</b>(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map</pre>
Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage.
<pre><b> HasComplexConjugate</b>(A::AlgEtQ) -> BoolElt</pre>
Returns if the algebra is the product of CM fields.
<pre><b> ComplexConjugate</b>(x::AlgEtQElt) -> AlgEtQElt</pre>
If A is a product of CM fields, it returns the complex conjugate of the argument.
<pre><b> IsConjugateStable</b>(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd</pre>
Given an order O in a CM-étale algebra, it returns wheter O is conjugate stable and the complex conjugate.
<pre><b> ComplexConjugate</b>(O::AlgEtQOrd) -> AlgEtQOrd</pre>
Given an order O in a CM-étale algebra, it returns the complex conjugate of the argument.
<pre><b> IsConjugateStable</b>(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl</pre>
Given a fractional ideal I in a CM-étale algebra, it returns wheter I is conjugate stable and the complex conjugate. Note: if the order of I is not conjugate stable, then the second output will be defined over the complex conjugate of the order.
<pre><b> ComplexConjugate</b>(I::AlgEtQIdl) -> AlgEtQIdl</pre>
If A is a product of CM fields, it returns the complex conjugate of the fractional ideal I. Note: if the order of I is not conjugate stable, then the output will be defined over the complex conjugate of the order.
<pre><b> CMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType
<pre><b> CreateCMType</b>(seq::SeqEnum[Map]) -> AlgEtQCMType</pre>
Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.
<pre><b> CMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive, that is, Im(phi(b))>0 for every phi in PHI.
<pre><b> CreateCMType</b>( b::AlgEtQElt  ) -> AlgEtQCMType</pre>
Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.
<pre><b> Print</b>( PHI :: AlgEtQCMType)</pre>
Print the AlgEtQCMType.
<pre><b> CMPositiveElement</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).
<pre><b> CMPosElt</b>( PHI::AlgEtQCMType )->AlgEtQElt</pre>
Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).
<pre><b> Homs</b>( PHI::AlgEtQCMType : prec:=Precision(GetDefaultRealField()) )->SeqEnum[Map]</pre>
Given a AlgEtQCMType PHI returns the sequence of maps to the complex field. The vararg prec (default value 30) determines the precision of the codomains of the maps.
<pre><b> 'eq'</b>(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=Precision(GetDefaultRealField()))->BoolElt</pre>
Returns whether two cm types are equal. This happens if and only if the quotient of (any) two CMPositiveElements is totally real and totally positive.
<pre><b> Precision</b>(PHI :: AlgEtQCMType)->RngIntElt </pre>
Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).
<pre><b> ChangePrecision</b>(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType</pre>
Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).
<pre><b> ChangePrecision</b>(~PHI :: AlgEtQCMType, prec::RngIntElt )</pre>
Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).
<pre><b> AllCMTypes</b>(A::AlgEtQ : Prec := Precision(GetDefaultRealField()) ) -> SeqEnum[AlgEtQCMType]</pre>
Returns all the AlgEtQCMTypes of A. The vararg Prec determined the precision of the codomain of the maps defining the CMTypes.
<pre><b> MinimalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals J subseteq I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subsetneq K subseteq I. Note J is never in the output.
<pre><b> IntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I. They are produced recursively from the minimal ones. The output includes I and J.
<pre><b> IntermediateIdealsWithPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I and (K:K)=S. They are produced recursively from the minimal ones. The output includes I, if (I:I)=S, and J, if (J:J)=S.
<pre><b> MaximalIntermediateIdeals</b>(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals J subseteq I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subseteq K subsetneq I. Note I is never in the output, while J is in the output if and only if the S-module I/J is simple, in which case the output consists only of J.
<pre><b> IntermediateIdealsWithTrivialExtension</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals I and J and an order O such that 
- S subset O,  
- J subseteq I, and 
- O subset (I:I),
returns all the fractional S-ideals K such that 
- J subseteq K subseteq I, and 
- O!!K = I. 
Note that the output always contains I. The output is produced by a recursiv use of MaximalIntermediateIdeals.
<pre><b> IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing</b>(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]</pre>
Given fractional S-ideals I and J and an order O such that 
- S subseteq O,  
- J subseteq I, and 
- O subseteq (I:I), 
returns all the fractional S-ideals K satisfying
- J subseteq K subseteq I, 
- O!!K = I, and 
- (K:K) eq S. 
In particular, the output contains I if and only if O = (I:I) = S. The output is produced by a recursive use of MaximalIntermediateIdeals.
<pre><b> IntermediateIdealsOfIndex</b>(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]</pre>
Given ideals J subseteq I over the same order, and a positive integer N, it returns all the ideals K such that 
- J subseteq K subseteq I, and 
- [I:K]=N. 
The output is produced by a recursive use of MaximalIntermediateIdeals.
<pre><b> IdealsOfIndex</b>(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.
<pre><b> IdealsOfIndex</b>(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]</pre>
Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.
<pre><b> IdealsOfIndex</b>(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]</pre>
Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.
<pre><b> IdealsOfIndex</b>(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".
<pre><b> IdealsOfIndex</b>(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]</pre>
Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".
<pre><b> ShortElement</b>(I::AlgEtQIdl) ->AlgEtQElt</pre>
Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by randomly picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).
<pre><b> SmallRepresentative</b>(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt</pre>
Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.
<pre><b> TwoGeneratingSet</b>(I::AlgEtQIdl)</pre>
A procedure that given an invertible ideal I put in the attibute Generators of I two non-zerodivisors in I that generate I. If I is known to be principal nothing is done.
<pre><b> ChineseRemainderTheorem</b>(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt</pre>
Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.
<pre><b> ChineseRemainderTheorem</b>(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt</pre>
Given a sequence `Is` of ideals of S, pairwise coprime, and a sequence `as` of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.
<pre><b> ChineseRemainderTheoremFunctions</b>(Is::SeqEnum[AlgEtQIdl])-> Map,Map</pre>
Given a sequence `Is` of N ideals of S, pairwise coprime, returns a function S->S^N representing the natural isomorphism S/&*(Is) -> \prod_(I in Is) S/I and a function S^N->S representing the inverse.
<pre><b> ResidueRingUnits</b>(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map</pre>
Returns the group (S/I)^* and a map (S/I)^* -> S. The MultiplicatorRing(I) must be the maximal order.
<pre><b> ResidueRingUnits</b>(I::AlgEtQIdl) -> GrpAb,Map</pre>
Returns the group (S/I)^* and a map (S/I)^* -> S, where S=Order(I) and the multiplicator ring of I is maximal.
<pre><b> ResidueRingUnitsSubgroupGenerators</b>(F::AlgEtQIdl) -> SeqEnum[AlgEtQElt]</pre>
Returns generators of (S/F)^* where F is an ideal of the order S.
<pre><b> IsPrincipal</b>(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgEtQElt</pre>
Return if the argument is a principal ideal; if so the function returns also the generator. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".
<pre><b> PicardGroup</b>( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes. The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".
<pre><b> ExtensionHomPicardGroups</b>(S::AlgEtQOrd,T::AlgEtQOrd : GRH:="false")->Map</pre>
Given orders S subseteq T, it returns the surjective extension map from PicardGroup(S) to PicardGroup(T). The vararg GRH, default false, is passed to PicardGroup.
<pre><b> UnitGroup</b>(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map</pre>
Return the unit group of a order in a etale algebra. The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".
<pre><b> IsIsomorphic</b>(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgEtQElt</pre>
Checks if I=x*J, for some x. If so, also x is returned. The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".
<pre><b> Factorization</b>(I::AlgEtQIdl) -> Tup</pre>
Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.
<pre><b> PrimesAbove</b>(I::AlgEtQIdl) -> SeqEnum[AlgEtQOrdIdl]</pre>
Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.
<pre><b> PlacesAboveRationalPrime</b>(E::AlgEtQ,p::RngIntElt)->SeqEnum[AlgEtQIdl]</pre>
Returns the maximal ideals of maximal order of the algebra E above the rational prime p.
<pre><b> SingularPrimes</b>(R::AlgEtQOrd) -> SeqEnum[AlgEtQOrdIdl]</pre>
Returns the non-invertible primes of the order.
<pre><b> NonInvertiblePrimes</b>(R::AlgEtQOrd) -> SetIndx</pre>
Returns the non-invertible primes of the order.
<pre><b> IsPrime</b>(I::AlgEtQIdl) -> BoolElt</pre>
Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.
<pre><b> IsMaximal</b>(I::AlgEtQIdl) -> BoolElt</pre>
Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.
<pre><b> IsMaximalIdeal</b>(I::AlgEtQIdl) -> BoolElt</pre>
Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is, a maximal S ideal.
<pre><b> IsBassAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
Check if the order is Bass at the prime ideal P, that is, if every overorder of S is Gorenstein at the primes above P.
<pre><b> IsBass</b>(S::AlgEtQOrd) -> BoolElt</pre>
Check if the order is Bass, that is, if every overorder of S is Gorenstein.
<pre><b> IsGorensteinAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt</pre>
Check if the order is Gorenstein at the prime ideal P, that is, if every fractional ideal I with (I:I)=S is locally principal at P.
<pre><b> IsGorenstein</b>(O::AlgEtQOrd)->BoolElt</pre>
Checks if the order O is Gorenstein, that is if the TraceDualIdeal of O is invertible, or equivalently, if all fractional ideals I with (I:I)=O are invertible.
<pre><b> Valuation</b>(x::AlgEtQElt,P::AlgEtQIdl)->RngIntElt</pre>
Given an element x and a maximal ideal P of the maximal order, it returns the valuation of x at P.
<pre><b> Valuation</b>(I::AlgEtQIdl,P::AlgEtQIdl)->RngIntElt</pre>
Given a fractional ideal I and a maximal ideal P, both of the maximal order, it returns the valuation of I at P.
<pre><b> InertiaDegree</b>(P::AlgEtQIdl)->RngIntElt</pre>
Given a maximal ideal P of the maximal order O above the rational prime p, it returns the inertia degree of P, that is, the index of the finite field extension GF(p)->O/P.
<pre><b> RamificationIndex</b>(P::AlgEtQIdl)->RngIntElt</pre>
Given a maximal ideal P of the maximal order O, it returns the reamification index of P.
<pre><b> Uniformizers</b>(PPs::SeqEnum[AlgEtQIdl])->SeqEnum</pre>
Given a sequence of maximal ideals P of the maximal order, it returns a sequence of elements t_P such that t_P is a uniformizer of P and a unit at every other prime in the sequence.
<pre><b> NonGorensteinPrimes</b>(S::AlgEtQOrd)->SeqEnum,SeqEnum</pre>
Given an order S it returns two sequences: the first containis the primes at which S is locally not Gorenstein; the second containis the CohenMacaulay types of S at these primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.
<pre><b> CohenMacaulayTypeAtPrime</b>(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt</pre>
Given an order S and a prime ideal P, it returns its Cohen-Macaulay Type at P. This integer equals the dimension of S^t/P*S^t where S^t is the trace dual of S.
<pre><b> CohenMacaulayType</b>(S::AlgEtQOrd)->RngIntElt</pre>
Given an order S returns its Cohen-Macaulay Type. This integer equals the max dimension of S^t/P*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.
<pre><b> WKICM_bar</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]</pre>
Returns representatives I of all weak equivalence classes, such that (I:I)=S. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.
<pre><b> WeakEquivalenceClassesWithPrescribedMultiplicatorRing</b>(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]</pre>
Returns representatives I of all weak equivalence classes, such that (I:I)=S. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.
<pre><b> WKICM</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum[AlgEtQIdl]</pre>
Returns a set of representatives of the weak equivalence class monoid of E. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.
<pre><b> WeakEquivalenceClassMonoid</b>(E::AlgEtQOrd : Method:="Auto")->SeqEnum[AlgEtQIdl]</pre>
Returns a set of representatives of the weak equivalence class monoid of E. The VarArg Method (default "Auto") is not used and kept for retrocompatibility.
<pre><b> IsWeakEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt</pre>
Checks if I and J are weakly equivalent, that is, if 1 \in (I:J)*(J:I), or equivalently, if I and J are locally equivalent at all prime of their common multiplicator ring. This function does not require that the ideals are defined over the same order.
<pre><b> IsWeakEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
Check if the two orders are weakly equivalent, that is equal.
<pre><b> IsWeakEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt</pre>
Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.
<pre><b> IsWeakEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt</pre>
Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.
<pre><b> IsWeaklyEquivalent</b>(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt</pre>
Checks if I and J are weakly equivalent, that is, if 1 \in (I:J)*(J:I), or equivalently, if I and J are locally equivalent at all prime of their common multiplicator ring. This function does not require that the ideals are defined over the same order.
<pre><b> IsWeaklyEquivalent</b>(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt</pre>
Check if the two orders are weakly equivalent, that is equal.
<pre><b> IsWeaklyEquivalent</b>(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt</pre>
Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.
<pre><b> IsWeaklyEquivalent</b>(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt</pre>
Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.
<pre><b> Print</b>(x::AlgEtQWECMElt)</pre>
Print the element.
<pre><b> IsCoercible</b>(W::AlgEtQWECM, x::.) -> BoolElt, .</pre>
Return whether the element is coercible into W and the result of the coercion if so.
<pre><b> '!'</b>(W::AlgEtQWECM, x::.) -> AlgEtQWECMElt</pre>
Coerce x into W, when possible.
<pre><b> 'in'</b>(x::AlgEtQWECMElt,W::AlgEtQWECM) -> BoolElt</pre>
Returns whether x is in W.
<pre><b> Parent</b>(x::AlgEtQWECMElt)->AlgEtQWECM</pre>
Returns the parent of x, that is, the weak equivalence class monoid it belongs to.
<pre><b> Ideal</b>(x::AlgEtQWECMElt)->AlgEtQIdl</pre>
Returns the ideal representing x.
<pre><b> MultiplicatorRing</b>(x::AlgEtQWECMElt)->AlgEtQOrd</pre>
Returns the multiplicator ring of x.
<pre><b> 'eq'</b>(x::AlgEtQWECMElt,y::AlgEtQWECMElt)->BoolElt</pre>
Returns whether the two elements define the same class.
<pre><b> SetRepresentative</b>(x::AlgEtQWECMElt,J::AlgEtQIdl)</pre>
It changes the Ideal attribute of x which consisits of the ideal representing the weak equivalence class to J, which must be weakly equivalent to the one of x.
<pre><b> MultiplicationTable</b>(W::AlgEtQWECM)->Assoc</pre>
Computes the multiplication table of W. It is returned as an associative array where, given two classes x and y of W, the value at the key \{x,y\} is x*y.
<pre><b> '*'</b>(x::AlgEtQWECMElt,y::AlgEtQWECMElt)->AlgEtQWECMElt</pre>
Returns weak equivalence class corresponding to the product.
<pre><b> '^'</b>(x::AlgEtQWECMElt,n::RngIntElt)->AlgEtQWECMElt</pre>
Given a weak equivalence class x and a non-negative integer n, returns the weak equivalence class x^n.
<pre><b> IsOne</b>(x::AlgEtQWECMElt)->BoolElt</pre>
Returns whether x is the neutral element of the weak equivalence class monoid it belongs to.
<pre><b> IsIdempotent</b>(x::AlgEtQWECMElt)->BoolElt</pre>
Returns whether x is an idempotent of weak equivalence class monoid it belongs to.
<pre><b> WeakEquivalenceClassMonoidAbstract</b>(R::AlgEtQOrd : Method:="Auto") -> AlgEtQWECM,Map</pre>
Returns the weak equivalence class monoid W of R together with a map w (with preimages) sending each class of W to a representative (determined by WeakEquivalenceClassMonoid). The vararg Methods is passed to WeakEquivalenceClassMonoid.
<pre><b> Print</b>(W::AlgEtQWECM)</pre>
Print the weak equivalence class monoid.
<pre><b> Order</b>(W::AlgEtQWECM)->AlgEtQOrd</pre>
Returns the order of W.
<pre><b> Array</b>(W::AlgEtQWECM)->Assoc</pre>
Return the underlying associative array of W, which is indexed by the overorder of R and values given by the weak equivalence classes with prescribed multiplicator ring.
<pre><b> Map</b>(W::AlgEtQWECM)->Map</pre>
Returns the map from W to the set of ideals which returns the representative of each class.
<pre><b> 'eq'</b>(W::AlgEtQWECM,WW::AlgEtQWECM)->BoolElt</pre>
Returns whether the two weak equivalence class monoid are the same, that is, if the underlying orders are.
<pre><b> '#'</b>(W::AlgEtQWECM)->RngInt</pre>
Returns the size of W.
<pre><b> Classes</b>(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]</pre>
Returns the sequence of the classes in W.
<pre><b> Representatives</b>(W::AlgEtQWECM)->SeqEnum[AlgEtQIdl]</pre>
Returns the sequence of representatives of the classes in W.
<pre><b> One</b>(W::AlgEtQWECM)->AlgEtQWECMElt</pre>
Returns the neutral element of W.
<pre><b> Random</b>(W::AlgEtQWECM)->AlgEtQWECMElt</pre>
Returns a random element of W.
<pre><b> Idempotents</b>(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]</pre>
Returns the sequence of the classes in W which are idempotent.
<pre><b> Localization</b>(W::AlgEtQWECM,P::AlgEtQIdl)->AlgEtQWECM</pre>
Given the weak equivalence class monoid of an order R in the etale algebra K and a maximal ideal P of R, returns the weak equivalance class monoid of the unique overorder of R which is locally equal to R at P and locally maximal at every other maximal ideal. This order is R+P^kO, where O is the maximal order of K and k is a non-negative integer big enough.
<pre><b> IsGeneratingSet</b>(W::AlgEtQWECM,seq::SeqEnum[AlgEtQWECMElt])->BoolElt</pre>
Given the weak equivalence class monoid of R and a sequence of classes in W, returns whether the sequence is a generating set of W.
<pre><b> IsGeneratingSet</b>(W::AlgEtQWECM,seq::SeqEnum[AlgEtQIdl])->BoolElt</pre>
Given the weak equivalence class monoid of R and a sequence of fractional R-ideals, returns whether the sequence represents a generating set of W.
<pre><b> ICM_bar</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
Returns the ideal classes of fractional S-ideals having MultiplicatorRing equal to S. This is the same as the orbit of the action of PicardGroup(S) on WKICM_bar(S).
<pre><b> ICM</b>(S::AlgEtQOrd : GRH:=false ) -> SeqEnum</pre>
Returns the ideal class monoid of the order, that is, a set of representatives for the isomorphism classes of the fractional S-ideals.
<pre><b> Print</b>(x::AlgEtQICMElt)</pre>
Print the element.
<pre><b> IsCoercible</b>(icm::AlgEtQICM, x::.) -> BoolElt, .</pre>
Return whether the element is coercible into icm and the result of the coercion if so.
<pre><b> '!'</b>(icm::AlgEtQICM, x::.) -> AlgEtQICMElt</pre>
Coerce x into icm, when possible.
<pre><b> 'in'</b>(x::AlgEtQICMElt,icm::AlgEtQICM) -> BoolElt</pre>
Returns whether x is in icm.
<pre><b> Parent</b>(x::AlgEtQICMElt)->AlgEtQICM</pre>
Returns the parent of x, that is, the ideal class monoid it belongs to.
<pre><b> WEClass</b>(x::AlgEtQICMElt)->AlgEtQWECMElt</pre>
Returns the weak equivalence class of x.
<pre><b> PicClass</b>(x::AlgEtQICMElt)->GrpAbElt,Map</pre>
Returns the elements of the PicardGroup Pic(T) of the multiplicator ring T of x together with the map from Pic(T) to the set of ideals corresponding to x.
<pre><b> Ideal</b>(x::AlgEtQICMElt)->AlgEtQIdl</pre>
Returns a deterministically computed ideal representing the ideal class. This is not stored in an attribute to save memory.
<pre><b> MultiplicatorRing</b>(x::AlgEtQICMElt)->AlgEtQOrd</pre>
Returns the multiplicator ring of x.
<pre><b> 'eq'</b>(x::AlgEtQICMElt,y::AlgEtQICMElt)->BoolElt</pre>
Returns whether the two elements define the same class.
<pre><b> '*'</b>(x::AlgEtQICMElt,y::AlgEtQICMElt)->AlgEtQICMElt</pre>
Returns ideal class corresponding to the product.
<pre><b> '^'</b>(x::AlgEtQICMElt,n::RngIntElt)->AlgEtQICMElt</pre>
Given an ideal class x and a non-negative integer n, returns the ideal class x^n.
<pre><b> IdealClassMonoidAbstract</b>(R::AlgEtQOrd : Method:="Auto") -> AlgEtQICM,Map</pre>
Returns the ideal class monoid icm of R together with a map m (with preimages) sending each class of icm to a representative (determined by WeakEquivalenceClassMonoid and PicardGroup). The vararg Methods is passed to WeakEquivalenceClassMonoid.
<pre><b> IsOne</b>(x::AlgEtQICMElt)->BoolElt</pre>
Returns whether x is the neutral element of the ideal class monoid it belongs to.
<pre><b> IsInvertibleInMultiplicatorRing</b>(x::AlgEtQICMElt)->BoolElt</pre>
Returns whether x the ideal class is invertible in its own multiplicator ring.
<pre><b> Print</b>(W::AlgEtQICM)</pre>
Print the ideal class monoid.
<pre><b> Order</b>(icm::AlgEtQICM)->AlgEtQOrd</pre>
Returns the order of icm.
<pre><b> Map</b>(W::AlgEtQICM)->Map</pre>
Returns the map from W to the set of ideals which returns the representative of each class.
<pre><b> 'eq'</b>(icm1::AlgEtQICM,icm2::AlgEtQICM)->BoolElt</pre>
Returns whether the two ideal class monoid are the same, that is, if the underlying orders are.
<pre><b> '#'</b>(icm::AlgEtQICM)->RngInt</pre>
Returns the size of W.
<pre><b> Classes</b>(icm::AlgEtQICM)->SeqEnum[AlgEtQICMElt]</pre>
Returns the sequence of the classes in icm.
<pre><b> Representatives</b>(icm::AlgEtQICM)->SeqEnum[AlgEtQIdl]</pre>
Returns the sequence of representatives of the classes in icm.
<pre><b> One</b>(icm::AlgEtQICM)->AlgEtQICMElt</pre>
Returns the neutral element of W.
<pre><b> Random</b>(icm::AlgEtQICM)->AlgEtQICMElt</pre>
Returns a random element of W.
<pre><b> PrintSeqAlgEtQElt</b>(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt</pre>
Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it returns a string that can be printed to file.
<pre><b> PrintWKICM</b>(R::AlgEtQOrd) -> MonStgElt</pre>
Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered from this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.
<pre><b> LoadWKICM</b>(str::MonStgElt) -> AlgEtQOrd</pre>
Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. These can be recovered with the approriate intrinsics.
