/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtQElt, 3;

declare attributes AlgEtQElt : Algebra, // AlgEtQ
                              AbsoluteCoordinates,
                              Components; // Tup

declare attributes AlgEtQ : Basis,
                           AbsoluteBasis,
                           PrimitiveElement,
                           PowerBasis,
                           OrthogonalIdempotents,
                           Idempotents;

//------------
// Printing 
//------------

intrinsic Print(x::AlgEtQElt)
{Print the element.}
    printf "%o", Components(x);
end intrinsic;

//------------
// Access attributes
//------------

intrinsic Parent(x::AlgEtQElt) -> AlgEtQ
{Returns the algebra to which the elemenet belongs to.}
  return x`Algebra;
end intrinsic;

intrinsic Algebra(x::AlgEtQElt) -> AlgEtQ
{Returns the algebra to which the elemenet belongs to.}
  return x`Algebra;
end intrinsic;

intrinsic Components(x::AlgEtQElt) -> SeqEnum
{Given an element, returns its components, which are elements of number fields.}
  return x`Components;
end intrinsic;

intrinsic AbsoluteCoordinates(x::AlgEtQElt) -> SeqEnum
{Given an element, returns the coordinates relative to the absolute basis, which are elements of the prime rational field.}
    if not assigned x`AbsoluteCoordinates then
        x`AbsoluteCoordinates:=&cat[ Flat(c) : c in Components(x) ];
    end if;
    return x`AbsoluteCoordinates;
end intrinsic;

//------------
// Coercion
//------------

function CreateAlgEtQElt(A,comp)
// given an AlgEtQ A and a Tup comp containing the components in the number fields of A, returns the corresponding elemnt
    x1:=New(AlgEtQElt);
    x1`Algebra:=A;
    x1`Components:=comp;
    return x1;
end function;

intrinsic IsCoercible(A::AlgEtQ, x::.) -> BoolElt, .
{Return whether the element is coercible into A and the result of the coercion if so.}
    if Parent(x) cmpeq A then
        return true,x;
    elif Type(x) in {List,Tup,SeqEnum} then
        comp:=<>;
        for i in [1..#x] do
            t,xi:=IsCoercible(Components(A)[i],x[i]); 
            if not t then
                return false,""; // early exit if not coercible
            else
                Append(~comp,xi);
            end if;
        end for;
        // now we know the elt is coercible in A
        x1:=CreateAlgEtQElt(A,comp);
        return true,x1;
    elif Type(x) eq RngIntElt or Type(x) eq FldRatElt then
        coordinates:=<>;
        nf:=Components(A);
        comp:=<nf[i]!x : i in [1..#nf]>; //diagonal embedding
        x1:=CreateAlgEtQElt(A,comp);
        return true,x1;
    else 
        return false,"";
    end if;
end intrinsic;

intrinsic '!'(A::AlgEtQ, x::.) -> AlgEtQElt
{Coerce x into A.}
    bool,x:=IsCoercible(A,x);
    require bool : "The element cannot be coerced in the algebra.";
    return x;
end intrinsic;

//------------
// Basic and Random elements
//------------

intrinsic One(A::AlgEtQ) -> AlgEtQElt
{The multiplicative neutral element of A.}   
    return A![1 : i in [1..#Components(A)]];
end intrinsic;

intrinsic Zero(A::AlgEtQ) -> AlgEtQElt
{The additive neutral element of A.}   
    return A![0 : i in [1..#Components(A)]];
end intrinsic;

intrinsic IsUnit(x::AlgEtQElt) -> BoolElt
{Returns wheter x is a unit in A.}   
    return forall{ c : c in Components(x) | c ne 0};
end intrinsic;

intrinsic IsZeroDivisor(x::AlgEtQElt) -> BoolElt
{Returns wheter x is a not unit in A.}   
    return not IsUnit(x);
end intrinsic;

intrinsic Random(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt
{Random element of A. The Coefficients are bounded by the positive integer bd.}   
    require bd gt 0 : "The bound needs to be a positive integer.";
    abs:=AbsoluteBasis(A);
    coeff:=[ Random(-bd,bd) : i in [1..#abs] ];
    num:=&+[abs[i]*coeff[i] : i in [1..#abs]];
    repeat
        abs:=AbsoluteBasis(A);
        coeff:=[ Random(-bd,bd) : i in [1..#abs] ];
        den:=&+[abs[i]*coeff[i] : i in [1..#abs]];
    until IsUnit(den);
    return num/den;
end intrinsic;

intrinsic Random(A::AlgEtQ : bd:=3) -> AlgEtQElt
{Random element of A. The Coefficients are bounded by the VarArg bd (default 3).}   
    return Random(A,bd);
end intrinsic;

intrinsic RandomUnit(A::AlgEtQ , bd::RngIntElt) -> AlgEtQElt
{Random unit of A. The Coefficients are bounded by the positive integer bd.}   
    require bd gt 0 : "The bound needs to be a positive integer.";
    nf:=Components(A);
    require BaseField(A) eq Rationals() : "The function works only for product of absolute number fiedls'" ;
    repeat
        num:=A!< E!Random(EquationOrder(E),bd) : E in nf>;
        den:=A!< E!Random(EquationOrder(E),bd) : E in nf>;
    until IsUnit(den) and IsUnit(num);
    return num/den;
end intrinsic;

//------------
// Equality and Operations
//------------

// eq 
intrinsic 'eq'(x1::AlgEtQElt,x2::AlgEtQElt) -> BoolElt
{Is x1=x2 ?}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    return (Components(x1) eq Components(x2)); // first part is faster for negatives 
                                                                        // and second is safer for positive
end intrinsic;

intrinsic 'eq'(x1::RngIntElt,x2::AlgEtQElt) -> BoolElt
{Is x1=x2 ?}
    return (Algebra(x2)!x1) eq x2;
end intrinsic;

intrinsic 'eq'(x1::FldRatElt,x2::AlgEtQElt) -> BoolElt
{Is x1=x2 ?}
    return (Algebra(x2)!x1) eq x2;
end intrinsic;

intrinsic 'eq'(x1::AlgEtQElt,x2::RngIntElt) -> BoolElt
{Is x1=x2 ?}
    return x1 eq (Algebra(x1)!x2);
end intrinsic;

intrinsic 'eq'(x1::AlgEtQElt,x2::FldRatElt) -> BoolElt
{Is x1=x2 ?}
    return x1 eq (Algebra(x1)!x2);
end intrinsic;

// +
intrinsic '+'(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
{x1+x2.}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    x3:=CreateAlgEtQElt(A,< Components(x1)[i] + Components(x2)[i] : i in [1..#Components(A)] >);
    return x3;
end intrinsic;

intrinsic '+'(x1::.,x2::AlgEtQElt) -> AlgEtQElt
{x1+x2.}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 + x2;
end intrinsic;

intrinsic '+'(x1::AlgEtQElt,x2::.) -> AlgEtQElt
{x1+x2.}
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 + x2;
end intrinsic;

intrinsic '+'(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
{x1+x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1+c : c in Components(x2)>);
end intrinsic;

intrinsic '+'(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
{x1+x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1+c : c in Components(x2)>);
end intrinsic;

intrinsic '+'(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
{x1+x2.}
    return CreateAlgEtQElt(Algebra(x1),<c+x2 : c in Components(x1)>);
end intrinsic;

intrinsic '+'(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
{x1+x2.}
    return CreateAlgEtQElt(Algebra(x1),<c+x2 : c in Components(x1)>);
end intrinsic;

// -
intrinsic '-'(x::AlgEtQElt) -> AlgEtQElt
{-x.}
    A:=Parent(x);
    comp:=< - Components(x)[i] : i in [1..#Components(A)] >;
    y:=CreateAlgEtQElt(A,comp);
    return y;
end intrinsic;

intrinsic '-'(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
{x1-x2.}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    comp:=< Components(x1)[i] - Components(x2)[i] : i in [1..#Components(A)] >;
    x3:=CreateAlgEtQElt(A,comp);
    return x3;
end intrinsic;

intrinsic '-'(x1::.,x2::AlgEtQElt) -> AlgEtQElt
{x1-x2.}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 - x2;
end intrinsic;

intrinsic '-'(x1::AlgEtQElt,x2::.) -> AlgEtQElt
{x1-x2.}
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 - x2;
end intrinsic;

intrinsic '-'(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
{x1-x2.}
    return Algebra(x2)! <x1-c : c in Components(x2)>;
end intrinsic;

intrinsic '-'(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
{x1-x2.}
    return Algebra(x2)! <x1-c : c in Components(x2)>;
end intrinsic;

intrinsic '-'(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
{x1-x2.}
    return Algebra(x1)! <c-x2 : c in Components(x1)>;
end intrinsic;

intrinsic '-'(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
{x1-x2.}
    return Algebra(x1)! <c-x2 : c in Components(x1)>;
end intrinsic;

// *
intrinsic '*'(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
{x1*x2.}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    comp:=< Components(x1)[i] * Components(x2)[i] : i in [1..#Components(A)] >;
    x3:=CreateAlgEtQElt(A,comp);
    return x3;
end intrinsic;

intrinsic '*'(x1::.,x2::AlgEtQElt) -> AlgEtQElt
{x1*x2.}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 * x2;
end intrinsic;

intrinsic '*'(x1::AlgEtQElt,x2::.) -> AlgEtQElt
{x1*x2.}
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 * x2;
end intrinsic;

intrinsic '*'(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
{x1*x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1*c : c in Components(x2)>);
end intrinsic;

intrinsic '*'(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
{x1*x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1*c : c in Components(x2)>);
end intrinsic;

intrinsic '*'(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
{x1*x2.}
    return CreateAlgEtQElt(Algebra(x1),<c*x2 : c in Components(x1)>);
end intrinsic;

intrinsic '*'(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
{x1*x2.}
    return CreateAlgEtQElt(Algebra(x1),<c*x2 : c in Components(x1)>);
end intrinsic;

// inverse and power
intrinsic Inverse(x::AlgEtQElt) -> AlgEtQElt
{1/x.}
    require IsUnit(x) : "The element is not invertible.";
    A:=Parent(x);
    comp:=< 1/(Components(x)[i]) : i in [1..#Components(A)] >;
    y:=CreateAlgEtQElt(A,comp);
    return y;
end intrinsic;

intrinsic '^'(x::AlgEtQElt,n::RngIntElt) -> AlgEtQElt
{x^n.}
    A:=Algebra(x);
    if n eq 0 then
        return One(A);
    elif n gt 0 then
        comp:=< Components(x)[i]^n : i in [1..#Components(A)] >;
        y:=CreateAlgEtQElt(A,comp);
        return y;
    elif n lt 0 then
        require IsUnit(x) : "The element is not invertible.";
        comp:=< Components(x)[i]^n : i in [1..#Components(A)] >;
        y:=CreateAlgEtQElt(A,comp);
        return y;
    end if;
end intrinsic;

// /
intrinsic '/'(x1::AlgEtQElt,x2::AlgEtQElt) -> AlgEtQElt
{x1/x2.}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    require IsUnit(x2) : "The denominator is not invertible.";
    comp:=< Components(x1)[i]/Components(x2)[i] : i in [1..#Components(A)] >;
    y:=CreateAlgEtQElt(A,comp);
    return y;
end intrinsic;

intrinsic '/'(x1::.,x2::AlgEtQElt) -> AlgEtQElt
{x1/x2.}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 / x2;
end intrinsic;

intrinsic '/'(x1::AlgEtQElt,x2::.) -> AlgEtQElt
{x1/x2.}
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 / x2;
end intrinsic;

intrinsic '/'(x1::RngIntElt,x2::AlgEtQElt) -> AlgEtQElt
{x1/x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1/c : c in Components(x2)>);
end intrinsic;

intrinsic '/'(x1::FldRatElt,x2::AlgEtQElt) -> AlgEtQElt
{x1/x2.}
    return CreateAlgEtQElt(Algebra(x2),<x1/c : c in Components(x2)>);
end intrinsic;

intrinsic '/'(x1::AlgEtQElt,x2::RngIntElt) -> AlgEtQElt
{x1/x2.}
    return CreateAlgEtQElt(Algebra(x1),<c/x2 : c in Components(x1)>);
end intrinsic;

intrinsic '/'(x1::AlgEtQElt,x2::FldRatElt) -> AlgEtQElt
{x1/x2.}
    return CreateAlgEtQElt(Algebra(x1),<c/x2 : c in Components(x1)>);
end intrinsic;

//------------
// Sums and Products of sequences
//------------

intrinsic '&+'(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt
{Given a sequence of AlgEtQElt returns the sum of the entries.}
    if #seq eq 0 then
        out:=Zero(Universe(seq));
    elif #seq eq 1 then
        out:=seq[1];
    else
        comps:=<>;
        N:=#Components(seq[1]);
        comps:=< &+[ Components(s)[i] : s in seq ] : i in [1..N] >;
        out:=Algebra(seq[1])!comps;
    end if;
    return out;
end intrinsic;

intrinsic '&*'(seq::SeqEnum[AlgEtQElt]) -> AlgEtQElt
{Given a sequence of AlgEtQElt returns the product of the entries.}
    if #seq eq 0 then
        out:=One(Universe(seq));
    elif #seq eq 1 then
        out:=seq[1];
    else
        comps:=<>;
        N:=#Components(seq[1]);
        comps:=< &*[ Components(s)[i] : s in seq ] : i in [1..N] >;
        out:=Algebra(seq[1])!comps;
    end if;
    return out;
end intrinsic;

intrinsic SumOfProducts(as::SeqEnum[AlgEtQElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
{Given sequences as and bs, such that #as eq #bs, returns `&+[as[i]*bs[i] : i in [1..#as]]`. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.}
    A:=Universe(bs);
    N:=#as;
    require N eq #bs : "The sequences don't have the same size";
    as:=[ Components(ai) : ai in as ];
    bs:=[ Components(bi) : bi in bs ];
    N_comp:=#bs[1];
    return A ! < &+[ as[i][c]*bs[i][c] : i in [1..N] ] : c in [1..N_comp]>;
end intrinsic;

intrinsic SumOfProducts(as::SeqEnum[RngIntElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
{Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.}
    A:=Universe(bs);
    N:=#as;
    require N eq #bs : "The sequences don't have the same size";
    bs:=[ Components(bi) : bi in bs ];
    N_comp:=#bs[1];
    return CreateAlgEtQElt(A,< &+[ as[i]*bs[i][c] : i in [1..N] ] : c in [1..N_comp]>);
end intrinsic;

intrinsic SumOfProducts(as::SeqEnum[FldRatElt],bs::SeqEnum[AlgEtQElt]) -> AlgEtQElt
{Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.}
    A:=Universe(bs);
    N:=#as;
    require N eq #bs : "The sequences don't have the same size";
    bs:=[ Components(bi) : bi in bs ];
    N_comp:=#bs[1];
    return CreateAlgEtQElt(A,< &+[ as[i]*bs[i][c] : i in [1..N] ] : c in [1..N_comp]>);
end intrinsic;

intrinsic SumOfProducts(as::SeqEnum[AlgEtQElt],bs::SeqEnum[RngIntElt]) -> AlgEtQElt
{Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.}
    A:=Universe(as);
    N:=#as;
    require N eq #bs : "The sequences don't have the same size";
    as:=[ Components(ai) : ai in as ];
    N_comp:=#as[1];
    return CreateAlgEtQElt(A,< &+[ as[i][c]*bs[i] : i in [1..N] ] : c in [1..N_comp]>);
end intrinsic;

intrinsic SumOfProducts(as::SeqEnum[AlgEtQElt],bs::SeqEnum[FldRatElt]) -> AlgEtQElt
{Given sequences as and bs, such that #as eq #bs, returns &+[as[i]*bs[i] : i in [1..#as]]. This intrinsic is included to obviate to the loss of efficiency due to the many calls of IsCoercible.}
    A:=Universe(as);
    N:=#as;
    require N eq #bs : "The sequences don't have the same size";
    as:=[ Components(ai) : ai in as ];
    N_comp:=#as[1];
    return CreateAlgEtQElt(A,< &+[ as[i][c]*bs[i] : i in [1..N] ] : c in [1..N_comp]>);
end intrinsic;

//------------
// Minimal polynomials and integrability
//------------

intrinsic MinimalPolynomial(x::AlgEtQElt) -> RngUPolElt
{Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.}
    A:=Algebra(x);
    require HasBaseField(A) : "The number fields of A shoud all be defined over the same base ring.";
    m:=LCM( [ MinimalPolynomial(c) : c in Components(x) ] );
    assert2 (0 eq Evaluate(m,x));
    return m;
end intrinsic;

intrinsic MinimalPolynomial(x::AlgEtQElt, F::Rng) -> RngUPolElt
{Returns the minimal polynommial over the ring F of the element x.}
    m:=LCM( [ MinimalPolynomial(c,F) : c in Components(x) ] );
    assert2 (0 eq Evaluate(m,x));
    return m;
end intrinsic;

intrinsic AbsoluteMinimalPolynomial(x::AlgEtQElt) -> RngUPolElt
{Returns the minimal polynommial over the prime field of the element x.}
    return MinimalPolynomial(x,PrimeField(Algebra(x)));
end intrinsic;

intrinsic IsIntegral(x::AlgEtQElt) -> BoolElt
{Returns whether the element x is integral (over the integers).}
    m:=AbsoluteMinimalPolynomial(x);
    return IsMonic(m) and IsCoercible(PolynomialRing(Integers()),m);
end intrinsic;

intrinsic Evaluate(f::RngUPolElt,a::AlgEtQElt) -> AlgEtQElt
{Evaluate the polynomial f at the element a.}
    A:=Algebra(a);
    coeff:=Coefficients(f); 
    pow_a:=[ i gt 1 select Self(i-1)*a else One(A) : i in [1..#coeff]];
    //Self makes it much more efficient
    return A!&+[ coeff[i]*pow_a[i] : i in [1..#coeff] ];
end intrinsic;

//------------
// Primitive Element and Power Basis
//------------

intrinsic PrimitiveElement(A::AlgEtQ) -> AlgEtQElt
{Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.}
    if not assigned A`PrimitiveElement then
        nf,embs:=Components(A);
        require #Seqset(nf) eq #nf : "The number fields defining the algebra are not distinct. Hence the algebra does not have a primitive element";
        A`PrimitiveElement:= CreateAlgEtQElt(A,<PrimitiveElement(F) : F in nf>);
    end if;
    return A`PrimitiveElement;
end intrinsic;

intrinsic PowerBasis(A::AlgEtQ) -> SeqEnum[AlgEtQElt]
{Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A.}
    if not assigned A`PowerBasis then
        a:=PrimitiveElement(A);
        pow_a:=[ i gt 1 select Self(i-1)*a else One(A) : i in [1..Dimension(A)]];
        //Self makes it much more efficient
        A`PowerBasis:=pow_a;
    end if;
    return A`PowerBasis;
end intrinsic;

//------------
// Basis and Vector Representations
//------------

intrinsic Basis(A::AlgEtQ) -> SeqEnum
{Returns a basis of the algebra over the common base field.}
    if not assigned A`Basis then
        nf,embs:=Components(A);
        require HasBaseField(A) : "The number fields do not have a common base field.";
        cc:= &cat[ [embs[i](b) : b in Basis(nf[i])] : i in [1..#nf]]; 
        A`Basis := [ A ! c : c in cc ];
    end if;
    return A`Basis;
end intrinsic;

intrinsic AbsoluteBasis(A::AlgEtQ) -> SeqEnum
{Returns a basis of the algebra over the prime field.}
    if not assigned A`AbsoluteBasis then
        nf,embs:=Components(A);
        cc:= &cat[ [embs[i](b) : b in AbsoluteBasis(nf[i])] : i in [1..#nf]]; 
        A`AbsoluteBasis := [ A ! c : c in cc ];
    end if;
    return A`AbsoluteBasis;
end intrinsic;

intrinsic AbsoluteCoordinates(seq::SeqEnum[AlgEtQElt] , basis::SeqEnum[AlgEtQElt]) -> SeqEnum
{Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.}
    A:=Algebra(seq[1]);
    require #basis eq AbsoluteDimension(A) : "the elements given do not form a basis";
    Mbasis:=Matrix([ AbsoluteCoordinates(b) : b in basis ]);
    Mseq:=[ Matrix([AbsoluteCoordinates(b)]) : b in seq ];
    Mcoeff:=[ Solution(Mbasis,S) : S in Mseq ];
    out:=[ Eltseq(C) : C in Mcoeff ];
    vprintf AlgEtQElt,3 : "Mbasis=\n%o\nMseq=\n%o\nMcoeff=\n%o\nout=\n%o\n",Mbasis,Mseq,Mcoeff,out;
    assert2 forall{i : i in [1..#seq] | seq[i] eq &+[out[i][j]*basis[j] : j in [1..#basis]]};
    return out;
end intrinsic;

//------------
// Idempotents
//------------

intrinsic OrthogonalIdempotents(A::AlgEtQ) -> SeqEnum
{Returns the orthogonal ideampotent element of the étale algebra A.}
    if not assigned A`OrthogonalIdempotents then
        nf,embs:=Components(A);
        A`OrthogonalIdempotents := [embs[i](One(nf[i])) : i in [1..#nf]];
    end if;
    return A`OrthogonalIdempotents;
end intrinsic;

intrinsic Idempotents(A::AlgEtQ) -> SeqEnum
{Returns the ideampotent element of the étale algebra A.}
    if not assigned A`Idempotents then
        ortid:=Seqset(OrthogonalIdempotents(A));
        ss:=Subsets(ortid);
        id:=[ Zero(A) + &+(Setseq(S)) : S in ss ];
        assert One(A) in id;
        A`Idempotents := id;
    end if;
    return A`Idempotents;
end intrinsic;

/* TESTS

    printf "### Testing Elements:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQElt",2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq f;
    _:=A!1;
    _:=A!(1/2);
    _:=Random(A,3)+Random(A,3);
    _:=Random(A,3)-Random(A,3);
    _:=Random(A,3)*Random(A,3);
    _:=Random(A,3)/RandomUnit(A,3);
    _:=Random(A);
    assert not IsIntegral(A!1/2);
    assert IsIntegral(A!2);

    ort:=OrthogonalIdempotents(A);
    idem:=Idempotents(A);
    assert ort subset idem;
    assert One(A) in idem;
    assert Zero(A) in idem;
    assert forall{i : i in idem | i^2 eq i};
    assert forall{i : i,j in [1..#ort] | ort[i]*ort[j] eq ort[i]*KroneckerDelta(i,j)};

    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
        printf ".";
    end for;

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq DefiningPolynomial(A);
    _:=A!1;
    _:=A!(1/2);
    _:=A!<seq[1]!1,seq[2]!(1/2)>;
    _:=A![1,2];
    assert A!1 + A!(1/2) eq A!3/2;
    e:=A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    assert not IsIntegral(e);
    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
        printf ".";
    end for;

    // testing sequences
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    seq:=[Random(A) : i in [1..10^3]];
    c0:=&+seq;
    c:=seq[1];
    for i in [2..#seq] do
        c+:=seq[i];
    end for;
    assert c eq c0;
    c0:=&*seq;
    c:=seq[1];
    for i in [2..#seq] do
        c*:=seq[i];
    end for;
    assert c eq c0;

    s1:=&+[seq[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq,seq);
    assert s1 eq s2;
    seq0:=[Random(-10,10) : i in [1..#seq]];
    s1:=&+[seq0[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq0,seq);
    s3:=&+[seq[i]*seq0[i] : i in [1..#seq]];
    s4:=SumOfProducts(seq,seq0);
    assert s1 eq s2 and s1 eq s3 and s1 eq s4;
    seq1:=[Random(-10,10)/Random(1,20) : i in [1..#seq]];
    s1:=&+[seq1[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq1,seq);
    s3:=&+[seq[i]*seq1[i] : i in [1..#seq]];
    s4:=SumOfProducts(seq,seq1);
    assert s1 eq s2 and s1 eq s3 and s1 eq s4;

    printf " all good!\n"; 

*/
