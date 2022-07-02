/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtElt, 3;

/*TODO:

*/

declare attributes AlgEtElt : Algebra, // AlgEt
                              AbsoluteCoordinates,
                              Hash,
                              Components; // Tup

declare attributes AlgEt : Basis,
                           AbsoluteBasis,
                           PrimitiveElement,
                           PowerBasis,
                           OrthogonalIdempotents,
                           Idempotents;

//------------
// Printing for AlgEtElt
//------------

intrinsic Print(x::AlgEtElt)
{Print the AlgEtElt.}
    printf "%o", Components(x);
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

intrinsic Components(x::AlgEtElt) -> SeqEnum
{Given an element x returns its components, which are elements of number fields.}
  return x`Components;
end intrinsic;

intrinsic AbsoluteCoordinates(x::AlgEtElt) -> SeqEnum
{Given an element x returns the coordinates relative to the absolute basis, which are elements of the prime field.}
    if not assigned x`AbsoluteCoordinates then
        x`AbsoluteCoordinates:=&cat[ Flat(c) : c in Components(x) ];
    end if;
    return x`AbsoluteCoordinates;
end intrinsic;

//------------
// Coercion
//------------

function CreateAlgEtElt(A,comp)
// given an AlgEt A and a Tup comp containing the components in the number fields of A, returns the corresponding elemnt
    x1:=New(AlgEtElt);
    x1`Algebra:=A;
    x1`Components:=comp;
    x1`Hash:=Hash(comp);
    return x1;
end function;

intrinsic IsCoercible(A::AlgEt, x::.) -> BoolElt, .
{Return whether x is coercible into A and the result of the coercion if so.}
    if Parent(x) cmpeq A then
        return true,x;
    elif Type(x) in {List,Tup,SeqEnum} then
        comp:=<>;
        for i in [1..#x] do
            t,xi:=IsCoercible(NumberFields(A)[i],x[i]); 
            if not t then
                return false,""; // early exit if not coercible
            else
                Append(~comp,xi);
            end if;
        end for;
        // now we know the elt is coercible in A
        x1:=CreateAlgEtElt(A,comp);
        return true,x1;
    elif Type(x) eq RngIntElt or Type(x) eq FldRatElt then
        coordinates:=<>;
        nf:=NumberFields(A);
        comp:=<nf[i]!x : i in [1..#nf]>; //diagonal embedding
        x1:=CreateAlgEtElt(A,comp);
        return true,x1;
    else 
        return false,"";
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
    return forall{ c : c in Components(x) | c ne 0};
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

intrinsic Random(A::AlgEt : bd:=3) -> AlgEtElt
{Random element of A. The Coefficients are bounded by VarArg bd (default 3).}   
    return Random(A,bd);
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
    return (x1`Hash eq x2`Hash) and (Components(x1) eq Components(x2)); // first part is faster for negatives 
                                                                        // and second is safer for positive
end intrinsic;

intrinsic 'eq'(x1::RngIntElt,x2::AlgEtElt) -> BoolElt
{Is x1=x2 ?}
    return (Algebra(x2)!x1) eq x2;
end intrinsic;

intrinsic 'eq'(x1::FldRatElt,x2::AlgEtElt) -> BoolElt
{Is x1=x2 ?}
    return (Algebra(x2)!x1) eq x2;
end intrinsic;

intrinsic 'eq'(x1::AlgEtElt,x2::RngIntElt) -> BoolElt
{Is x1=x2 ?}
    return x1 eq (Algebra(x1)!x2);
end intrinsic;

intrinsic 'eq'(x1::AlgEtElt,x2::FldRatElt) -> BoolElt
{Is x1=x2 ?}
    return x1 eq (Algebra(x1)!x2);
end intrinsic;

intrinsic '+'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1+x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    x3:=CreateAlgEtElt(A,< Components(x1)[i] + Components(x2)[i] : i in [1..#NumberFields(A)] >);
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
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 + x2;
end intrinsic;

intrinsic '-'(x::AlgEtElt) -> AlgEtElt
{-x}
    A:=Parent(x);
    comp:=< - Components(x)[i] : i in [1..#NumberFields(A)] >;
    y:=CreateAlgEtElt(A,comp);
    return y;
end intrinsic;

intrinsic '-'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1-x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    comp:=< Components(x1)[i] - Components(x2)[i] : i in [1..#NumberFields(A)] >;
    x3:=CreateAlgEtElt(A,comp);
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
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 - x2;
end intrinsic;

intrinsic '*'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1*x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    comp:=< Components(x1)[i] * Components(x2)[i] : i in [1..#NumberFields(A)] >;
    x3:=CreateAlgEtElt(A,comp);
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
    bool,x2:=IsCoercible(Algebra(x1),x2);
    require bool : "x2 not coercible";
    return x1 * x2;
end intrinsic;

intrinsic Inverse(x::AlgEtElt) -> AlgEtElt
{1/x}
    require IsUnit(x) : "The element is not invertible.";
    A:=Parent(x);
    comp:=< 1/(Components(x)[i]) : i in [1..#NumberFields(A)] >;
    y:=CreateAlgEtElt(A,comp);
    return y;
end intrinsic;

intrinsic '^'(x::AlgEtElt,n::RngIntElt) -> AlgEtElt
{x^n}
    A:=Algebra(x);
    if n eq 0 then
        return One(A);
    elif n gt 0 then
        comp:=< Components(x)[i]^n : i in [1..#NumberFields(A)] >;
        y:=CreateAlgEtElt(A,comp);
        return y;
    elif n lt 0 then
        require IsUnit(x) : "The element is not invertible.";
        comp:=< Components(x)[i]^n : i in [1..#NumberFields(A)] >;
        y:=CreateAlgEtElt(A,comp);
        return y;
    end if;
end intrinsic;

intrinsic '/'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1/x2}
    A:=Parent(x1);
    require A cmpeq Parent(x2): "The elements must belong to the same algebra.";
    require IsUnit(x2) : "The denominator is not invertible.";
    comp:=< Components(x1)[i]/Components(x2)[i] : i in [1..#NumberFields(A)] >;
    y:=CreateAlgEtElt(A,comp);
    return y;
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
    require HasBaseField(A) : "The number fields of A shoud all be defined over the same base ring.";
    m:=LCM( [ MinimalPolynomial(c) : c in Components(x) ] );
    assert2 (0 eq Evaluate(m,x));
    return m;
end intrinsic;

intrinsic MinimalPolynomial(x::AlgEtElt, F::Rng) -> RngUPolElt
{Returns the minimal polynommial over the ring F of the element x.}
    m:=LCM( [ MinimalPolynomial(c,F) : c in Components(x) ] );
    assert2 (0 eq Evaluate(m,x));
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

intrinsic Evaluate(f::RngUPolElt,a::AlgEtElt) -> AlgEtElt
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

intrinsic PrimitiveElement(A::AlgEt) -> AlgEtElt
{Returns the primitive element of the étale algebra A. Note that A has a primitive element only if it is the product of distinct number fields.}
    if not assigned A`PrimitiveElement then
        nf,embs:=NumberFields(A);
        require #Seqset(nf) eq #nf : "The number fields defining the algebra are not distinct. Hence the algebra does not have a primitive element";
        A`PrimitiveElement:= CreateAlgEtElt(A,<PrimitiveElement(F) : F in nf>);
    end if;
    return A`PrimitiveElement;
end intrinsic;

intrinsic PowerBasis(A::AlgEt) -> SeqEnum[AlgEtElt]
{Returns the power basis of the étale algebra A, consisting of powers of the PrimitiveElement of A}
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

intrinsic AbsoluteCoordinates(seq::SeqEnum[AlgEtElt] , basis::SeqEnum[AlgEtElt]) -> SeqEnum
{Given a sequence of elements and a basis over the PrimeField returns a sequence whose entries are the coordinates in the PrimeField with respect to the given basis.}
    A:=Algebra(seq[1]);
    require #basis eq AbsoluteDimension(A) : "the elements given do not form a basis";
    Mbasis:=Matrix([ AbsoluteCoordinates(b) : b in basis ]);
    Mseq:=[ Matrix([AbsoluteCoordinates(b)]) : b in seq ];
    Mcoeff:=[ Solution(Mbasis,S) : S in Mseq ];
    out:=[ Eltseq(C) : C in Mcoeff ];
    vprintf AlgEtElt,3 : "Mbasis=\n%o\nMseq=\n%o\nMcoeff=\n%o\nout=\n%o\n",Mbasis,Mseq,Mcoeff,out;
    assert2 forall{i : i in [1..#seq] | seq[i] eq &+[out[i][j]*basis[j] : j in [1..#basis]]};
    return out;
end intrinsic;


//------------
// Idempotents
//------------

intrinsic OrthogonalIdempotents(A::AlgEt) -> SeqEnum
{Returns the orthogonal ideampotent element of the étale algebra A}
    if not assigned A`OrthogonalIdempotents then
        nf,embs:=NumberFields(A);
        A`OrthogonalIdempotents := [embs[i](One(nf[i])) : i in [1..#nf]];
    end if;
    return A`OrthogonalIdempotents;
end intrinsic;

intrinsic Idempotents(A::AlgEt) -> SeqEnum
{Returns the ideampotent element of the étale algebra A}
    if not assigned A`Idempotents then
        ortid:=Seqset(OrthogonalIdempotents(A));
        ss:=Subsets(ortid);
        id:=[ Zero(A) + &+(Setseq(S)) : S in ss ];
        assert One(A) in id;
        A`Idempotents := id;
    end if;
    return A`Idempotents;
end intrinsic;

/* CONTINUE FROM HERE

intrinsic CoordinatesOverBaseRing(seq::SeqEnum[AlgEtElt] , basis::SeqEnum[AlgEtElt]) -> SeqEnum
{//TODO I am not sure what this function should do}
    return ;
end intrinsic;

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

    AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtElt",2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq f;
    A!1,A!(1/2);
    Random(A,3)+Random(A,3);
    Random(A,3)-Random(A,3);
    Random(A,3)*Random(A,3);
    Random(A,3)/RandomUnit(A,3);
    Random(A);
    IsIntegral(A!1/2);
    IsIntegral(A!2);

    OrthogonalIdempotents(A);
    Idempotents(A);

    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
    end for;

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq DefiningPolynomial(A);
    A!1,A!(1/2),A!<seq[1]!1,seq[2]!(1/2)>,A![1,2];
    A!1 + A!(1/2);
    e:=A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    IsIntegral(e);
    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
    end for;

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq p;
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
    OrthogonalIdempotents(A);
    Idempotents(A);

    A:=EtaleAlgebra([E,E]);
    assert #Basis(A) eq Dimension(A);
    assert #AbsoluteBasis(A) eq AbsoluteDimension(A);
    OrthogonalIdempotents(A);
    Idempotents(A);

*/
