/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEt, 1;

/*TODO:

*/

declare type AlgEt[AlgEtElt];

declare attributes AlgEt : DefiningPolynomial, 
                           // ass_algebra, 
                           NumberFields;

//------------
// Creation and Printig for AlgEt
//------------

intrinsic Print(A::AlgEt)
{Prints the defining polynomial or the number fields defining A.}
    if assigned A`DefiningPolynomial then
        printf "Etale Algebra over QQ defined by %o", DefiningPolynomial(A);
    else
        printf "Etale Algebra product of %o", NumberFields(A);
    end if;
end intrinsic;

intrinsic EtaleAlgebra(seq::SeqEnum[FldNum]) -> AlgEt
{Given a sequence of number fields returns the étale algebra corresponding to the direct product.}
    A:=New(AlgEt);
    A`NumberFields:=[seq];
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt) -> AlgEt
{Given a squarefree polynomial returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    return EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
end intrinsic;

//------------
// Access attributes
//------------

intrinsic DefiningPolynomial(A::AlgEt) -> RngUPolElt
{Returns the defining polynomial of A, if the corresponding number fields are distinct.}
    if not assigned A`DefiningPolynomial then
        polys:=[DefiningPolynomial(L) : L in NumberFields(A)];
        require #polys eq #Seqset(polys) : "The number fields defining A are not distinct.";
        A`DefiningPolynomial:=&*(polys);
    end if;
    return A`DefiningPolynomial;
end intrinsic;

intrinsic NumberFields(A::AlgEt) -> SeqEnum
{Returns the number fields of which A is a product of.}
    return A`NumberFields;
end intrinsic;




/*

// MOVE OR REMOVE THE NEXT
intrinsic HomsToC(A::AlgEt : Precision:=30)->SeqEnum[Map]
{returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30}
    CC:=ComplexField(Precision);
    images:=function(x)
        return &cat[[CC ! z : z in Conjugates(y : Precision:=Precision)] :y in Components(x)];
    end function;
    maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Degree(A)] >;
    f:=&*[DefiningPolynomial(L[1]) : L in A`NumberFields];
    assert &and [ Abs(Evaluate(f,g(PrimitiveElement(A)))) lt 10^-10 : g in maps]; //the precision here is quite arbitrary...
    return maps;
end intrinsic;

intrinsic OrthogonalIdempotents(A::AlgEt)->SeqEnum
{returns a sequence containing the orthogonal ideampotents of the algebra, that is the image of the units of the number fields the algebra is product of}
    return [L[2](One(L[1])) : L in A`NumberFields ];
end intrinsic;

intrinsic Dimension(A::AlgEt)->RngInt
{Dimenison of A}    
    return Dimension(AssociativeAlgebra(A));
end intrinsic;

intrinsic Degree(A::AlgEt)->RngInt
{Dimenison of A}    
    return Degree(AssociativeAlgebra(A));
end intrinsic;

intrinsic PrimitiveElement(A::AlgEt) -> AlgEtElt
{ it returns an element which corresponds to the class of X in Q[X]/(f(X)) }
    assert IsSquarefree(DefiningPolynomial(A));
    return &+[L[2](PrimitiveElement(L[1])) : L in A`NumberFields];
end intrinsic;

intrinsic Components(x::AlgEtElt) -> SeqEnum
{returns the components of the element as in the product of number fields}
    A:=Parent(x);
    idem:=OrthogonalIdempotents(A);
    x_asProd:=[**];
    for i in [1..#A`NumberFields] do
        L:=A`NumberFields[i];
        x_L:=(x*idem[i])@@L[2];
        Append(~x_asProd,x_L);
    end for;
    return x_asProd;
end intrinsic;

intrinsic Idempotents(A::AlgEt)->SeqEnum
{returns a sequence containing the ideampotents of the algebra, zero included}
    ort_idem:=OrthogonalIdempotents(A);
    cc:=CartesianProduct([[A!0,oi] : oi in ort_idem]);
    idem:=[&+([cj : cj in c]) : c in cc];
    return idem;
end intrinsic;

intrinsic Coordinates(seq::SeqEnum[AlgEtElt],basis::SeqEnum[AlgEtElt]) -> SeqEnum
{ the coordinates of the sequence S of elements in an etale algebra A, relative to the given basis of A over the rationals. }
    vprintf et_algebras: "Coordinates";
    require Universe(seq) eq Universe(basis) : "the sequences must be defined over the same algebra";
    A:=AssAlgebra(Universe(seq));
    seq:=[A ! x : x in seq ];
    basis:=[A ! x : x in basis ];
    coord:=Coordinates(seq,basis);
    return coord;
end intrinsic;

intrinsic IsZeroDivisor(x::AlgEtElt) -> BoolElt
{returns if the element x is a zero divisor}
    return exists{c : c in Components(x)|c eq Zero(Parent(c))};
end intrinsic;

intrinsic IsZeroDivisor2(x::AlgEtElt) -> BoolElt
{returns if the element x is a zero divisor}
    return IsZeroDivisor(x);
end intrinsic;

intrinsic Evaluate(p::RngUPolElt, y::AlgEtElt) -> AlgEtElt
{evaluates a polynomial at y}
    evl:=Evaluate(p,AssElt(y));
    return Algebra(y)!evl;
end intrinsic;

*/

/* TEST

    Attach("~/packages_github/AlgEt/AlgEt.m");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    A;

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    A;

    A:=EtaleAlgebra([NumberField(p),NumberField(x^2-5)]);
    A;


*/
