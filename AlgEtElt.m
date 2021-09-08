/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtElt, 1;

/*TODO:

*/

declare attributes AlgEtElt : Algebra, // AlgEt
                              AbsoluteCoordinates,
                              Coordinates; // Tup

declare attributes AlgEt : Basis,
                           AbsoluteBasis;

//------------
// Printing for AlgEtElt
//------------

intrinsic Print(x::AlgEtElt)
{Print the AlgEtElt.}
    printf "%o", Coordinates(x);
end intrinsic;

//------------
// Access attributes
//------------

intrinsic Parent(x::AlgEtElt) -> AlgEt
{Returns the algebra to which the elemenet belongs to.}
  return x`Algebra;
end intrinsic;

intrinsic Algebra(x::AlgEtElt) -> AlgEt
{Returns the algebra to which the elemenet belongs to.}
  return x`Algebra;
end intrinsic;

intrinsic Coordinates(x::AlgEtElt) -> SeqEnum
{Given an element x returns its coordinates, which are elements of number fields.}
  return x`Coordinates;
end intrinsic;

intrinsic AbsoluteCoordinates(x::AlgEtElt) -> SeqEnum
{Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.}
    if not assigned x`AbsoluteCoordinates then
        x`AbsoluteCoordinates:=&cat[ Flat(c) : c in Coordinates(x) ];
    end if;
    return x`AbsoluteCoordinates;
end intrinsic;

//------------
// Coercion
//------------

intrinsic IsCoercible(A::AlgEt, x::.) -> BoolElt, .
{Return whether x is coercible into A and the result of the coercion if so.}
    if Parent(x) cmpeq A then
        return true,x;
    elif Type(x) in {List,Tup,SeqEnum} then
        coordinates:=<>;
        for i in [1..#x] do
            t,xi:=IsCoercible(NumberFields(A)[i],x[i]); 
            if not t then
                return false,_; // early exit if not coercible
            else
                Append(~coordinates,xi);
            end if;
        end for;
        // now we know the elt is coercible in A
        x1:=New(AlgEtElt);
        x1`Algebra:=A;
        x1`Coordinates:=coordinates;
        return true,x1;
    elif Type(x) eq RngIntElt or Type(x) eq FldRatElt then
        coordinates:=<>;
        nf:=NumberFields(A);
        x1:=New(AlgEtElt);
        x1`Algebra:=A;
        x1`Coordinates:=<nf[i]!x : i in [1..#nf]>; //diagonal embedding
        return true,x1;
    else 
        return false,_;
    end if;
end intrinsic;

intrinsic '!'(A::AlgEt, x::.) -> AlgEtElt
{Coerce x into A.}
    bool,x:=IsCoercible(A,x);
    require bool : "The element cannot be coerced in the algebra.";
    return x;
end intrinsic;

//------------
// Basic and Random elements
//------------

intrinsic One(A::AlgEt) -> AlgEtElt
{The multiplicative neutral element of A.}   
    return A![1 : i in [1..#NumberFields(A)]];
end intrinsic;

intrinsic Zero(A::AlgEt) -> AlgEtElt
{The additive neutral element of A.}   
    return A![0 : i in [1..#NumberFields(A)]];
end intrinsic;

intrinsic IsUnit(x::AlgEtElt) -> BoolElt
{Returns wheter x is a unit in A.}   
    return forall{ c : c in Coordinates(x) | c ne 0};
end intrinsic;

intrinsic IsZeroDivisor(x::AlgEtElt) -> BoolElt
{Returns wheter x is a not unit in A.}   
    return not IsUnit(x);
end intrinsic;

//kept for retro-compatbility
intrinsic IsZeroDivisor2(x::AlgEtElt) -> BoolElt
{Returns wheter x is a not unit in A.}   
    return not IsUnit(x);
end intrinsic;

intrinsic Random(A::AlgEt , bd::RngIntElt) -> AlgEtElt
{Random element of A. The Coefficients are bounded by the positive integer bd.}   
    require bd gt 0 : "The bound needs to be a positive integer.";
    nf:=NumberFields(A);
    require BaseField(A) eq Rationals() : "The function works only for product of absolute number fiedls'" ;
    num:=A!< E!Random(EquationOrder(E),bd) : E in nf>;
    repeat
        den:=A!< E!Random(EquationOrder(E),bd) : E in nf>;
    until IsUnit(den);
    return num/den;
end intrinsic;

intrinsic RandomUnit(A::AlgEt , bd::RngIntElt) -> AlgEtElt
{Random unit of A. The Coefficients are bounded by the positive integer bd.}   
    require bd gt 0 : "The bound needs to be a positive integer.";
    nf:=NumberFields(A);
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

intrinsic 'eq'(x1::AlgEtElt,x2::AlgEtElt) -> BoolElt
{Is x1=x2 ?}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    return Coordinates(x1) eq Coordinates(x2);
end intrinsic;

intrinsic 'eq'(x1::.,x2::AlgEtElt) -> BoolElt
{Is x1=x2 ?}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 = x2;
    end intrinsic;

intrinsic 'eq'(x1::AlgEtElt,x2::.) -> BoolElt
{Is x1=x2 ?}
    bool,x2:=IsCoercible(Algebra(x2),x1);
    require bool : "x2 not coercible";
    return x1 = x2;
end intrinsic;

intrinsic '+'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1+x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    x3:=A!< Coordinates(x1)[i] + Coordinates(x2)[i] : i in [1..#NumberFields(A)] >;
    return x3;
end intrinsic;

intrinsic '+'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1+x2}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 + x2;
end intrinsic;

intrinsic '+'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1+x2}
    bool,x2:=IsCoercible(Algebra(x2),x1);
    require bool : "x2 not coercible";
    return x1 + x2;
end intrinsic;

intrinsic '-'(x::AlgEtElt) -> AlgEtElt
{-x}
    A:=Parent(x);
    y:=A!< - Coordinates(x)[i] : i in [1..#NumberFields(A)] >;
    return y;
end intrinsic;

intrinsic '-'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1-x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    x3:=A!< Coordinates(x1)[i] - Coordinates(x2)[i] : i in [1..#NumberFields(A)] >;
    return x3;
end intrinsic;

intrinsic '-'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1-x2}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 - x2;
end intrinsic;

intrinsic '-'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1-x2}
    bool,x2:=IsCoercible(Algebra(x2),x1);
    require bool : "x2 not coercible";
    return x1 - x2;
end intrinsic;

intrinsic '*'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1*x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    x3:=A!< Coordinates(x1)[i] * Coordinates(x2)[i] : i in [1..#NumberFields(A)] >;
    return x3;
end intrinsic;

intrinsic '*'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1*x2}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 * x2;
end intrinsic;

intrinsic '*'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1*x2}
    bool,x2:=IsCoercible(Algebra(x2),x1);
    require bool : "x2 not coercible";
    return x1 * x2;
end intrinsic;

intrinsic Inverse(x::AlgEtElt) -> AlgEtElt
{1/x}
    require IsUnit(x) : "The element is not invertible.";
    A:=Parent(x);
    y:=A!< 1/(Coordinates(x)[i]) : i in [1..#NumberFields(A)] >;
    return y;
end intrinsic;

intrinsic '^'(x::AlgEtElt,n::RngIntElt) -> AlgEtElt
{x^n}
    A:=Algebra(x);
    if n eq 0 then
        return One(A);
    elif n gt 0 then
        return A!< Coordinates(x)[i]^n : i in [1..#NumberFields(A)] >;
    elif n lt 0 then
        require IsUnit(x) : "The element is not invertible.";
        return A!< Coordinates(x)[i]^n : i in [1..#NumberFields(A)] >;
    end if;
end intrinsic;

intrinsic '/'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1/x2}
    require IsUnit(x2) : "The denominator is not invertible.";
    return x1*Inverse(x2);
end intrinsic;

intrinsic '/'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1/x2}
    bool,x1:=IsCoercible(Algebra(x2),x1);
    require bool : "x1 not coercible";
    return x1 / x2;
end intrinsic;

intrinsic '/'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1/x2}
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 / x2;
end intrinsic;


//------------
// Minimal polynomials and integrability
//------------

intrinsic MinimalPolynomial(x::AlgEtElt) -> RngUPolElt
{Returns the minimal polynommial over the common base ring of the number fields defining A of the element x.}
    A:=Algebra(x);
    nf:=NumberFields(A);
    require HasBaseField(A) : "The number fields of A shoud all be defined over the same base ring.";
    m:=LCM( [ MinimalPolynomial(c) : c in Coordinates(x) ] );
    return m;
end intrinsic;

intrinsic MinimalPolynomial(x::AlgEtElt, F::Rng) -> RngUPolElt
{Returns the minimal polynommial over the ring F of the element x.}
    m:=LCM( [ MinimalPolynomial(c,F) : c in Coordinates(x) ] );
    return m;
end intrinsic;

intrinsic AbsoluteMinimalPolynomial(x::AlgEtElt) -> RngUPolElt
{Returns the minimal polynommial over the prime field of the element x.}
    return MinimalPolynomial(x,PrimeField(Algebra(x)));
end intrinsic;

intrinsic IsIntegral(x::AlgEtElt) -> BoolElt
{Returns whether the element x is integral (over the integers).}
    m:=AbsoluteMinimalPolynomial(x);
    return IsMonic(m) and IsCoercible(PolynomialRing(Integers()),m);
end intrinsic;

//------------
// Basis and Vector Representations
//------------

intrinsic Basis(A::AlgEt) -> SeqEnum
{Returns a basis of the algebra over the common base field.}
    if not assigned A`Basis then
        nf,embs:=NumberFields(A);
        require HasBaseField(A) : "The number fields do not have a common base field.";
        cc:= &cat[ [embs[i](b) : b in Basis(nf[i])] : i in [1..#nf]]; 
        A`Basis := [ A ! c : c in cc ];
    end if;
    return A`Basis;
end intrinsic;

intrinsic AbsoluteBasis(A::AlgEt) -> SeqEnum
{Returns a basis of the algebra over the prime field.}
    if not assigned A`AbsoluteBasis then
        nf,embs:=NumberFields(A);
        cc:= &cat[ [embs[i](b) : b in AbsoluteBasis(nf[i])] : i in [1..#nf]]; 
        A`AbsoluteBasis := [ A ! c : c in cc ];
    end if;
    return A`AbsoluteBasis;
end intrinsic;

/* CONTINUE FROM HERE

*/


/*TO MOVE OR REMOVE

intrinsic '.'(A::AlgEt,i::RngIntElt) -> AlgEtElt
{A.i returns the ith generator of the Algebra A}
   alg:=AssociativeAlgebra(A);
   y:=New(AlgEtElt);
   y`AlgAssElt:=alg.i;
   y`Algebra:=A;
   return y;
end intrinsic;


*/

/* TEST

    Attach("~/packages_github/AlgEt/AlgEt.m");
    Attach("~/packages_github/AlgEt/AlgEtElt.m");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    A!1,A!(1/2);
    Random(A,3)+Random(A,3);
    Random(A,3)-Random(A,3);
    Random(A,3)*Random(A,3);
    Random(A,3)/RandomUnit(A,3);
    IsIntegral(A!1/2);
    IsIntegral(A!2);

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    A!1,A!(1/2),A!<seq[1]!1,seq[2]!(1/2)>,A![1,2];
    A!1 + A!(1/2);
    e:=A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    IsIntegral(e);

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    A!1 - A!(1/2);
    A!1 / A!(1/2);
    
    A:=EtaleAlgebra([K,K]);
    e:=A!<K.1,K.1^2>;
    MinimalPolynomial(e);
    IsIntegral(e);
    e2:=A!<K.1^2,K.1^2>;
    assert MinimalPolynomial(e2) eq MinimalPolynomial(K.1^2);
    _,embs,projs:=NumberFields(A);
    assert embs[1](K.1)+embs[2](K.1^2) eq e;
    assert projs[1](e) eq K.1;
    assert projs[2](e) eq K.1^2;

    E:=NumberField(p);
    seq:=[E,K];
    A:=EtaleAlgebra(seq);
    A!1 + A!(1/2);
    A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    assert #AbsoluteBasis(A) eq AbsoluteDimension(A);
    e:=A!<K.1+E.1,K.1^2>;
    assert #AbsoluteCoordinates(e) eq AbsoluteDimension(A);
    assert e eq &+[AbsoluteCoordinates(e)[i]*AbsoluteBasis(A)[i] : i in [1..AbsoluteDimension(A)]];

    A:=EtaleAlgebra([E,E]);
    assert #Basis(A) eq Dimension(A);
    assert #AbsoluteBasis(A) eq AbsoluteDimension(A);

*/
