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
                           Dimension,
                           AbsoluteDimension,
                           BaseField, //a tup : <F,m> where F is the Base field and m is the diagonal embedding into A
                           HasBaseField, //a boolean
                           PrimeField,
                           Components; //a tup of 3 sequences: the first are the NF, 
                                         //the second are embeddings and the third are projections

//------------
// Creation and Printing for AlgEt
//------------

intrinsic Print(A::AlgEt)
{Prints the defining polynomial or the number fields defining A.}
    if assigned A`DefiningPolynomial then
        f:=DefiningPolynomial(A);
        printf "Etale Algebra over (%o) defined by %o", FieldOfFractions(CoefficientRing(f)),f;    
    else
        printf "Etale Algebra product of %o", Components(A);
    end if;
end intrinsic;

intrinsic EtaleAlgebra(seq::SeqEnum[FldNum]) -> AlgEt
{Given a sequence of number fields returns the étale algebra corresponding to the direct product.}
    A:=New(AlgEt);
    embs:=[ map< seq[i]->A | x:-> A! (<seq[j]!0 : j in [1..i-1]> cat <x> cat <seq[j]!0 : j in [i+1..#seq]>)  > : i in [1..#seq] ];
    projs:=[ map< A->seq[i] | y:-> Components(y)[i] > : i in [1..#seq] ];
    A`Components:=<seq,embs,projs>;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt) -> AlgEt
{Given a squarefree polynomial returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

//------------
// Access attributes
//------------

intrinsic DefiningPolynomial(A::AlgEt) -> RngUPolElt
{Returns the defining polynomial of A, if the corresponding number fields are distinct.}
    if not assigned A`DefiningPolynomial then
        polys:=[DefiningPolynomial(L) : L in Components(A)];
        require #polys eq #Seqset(polys) : "The number fields defining A are not distinct.";
        A`DefiningPolynomial:=&*(polys);
    end if;
    return A`DefiningPolynomial;
end intrinsic;

intrinsic Components(A::AlgEt) -> SeqEnum
{Returns the number fields of which A is a product of,together with embeddings and projections}
    return Explode(A`Components);
end intrinsic;

intrinsic Dimension(A::AlgEt)->RngInt
{Dimension of A}    
    if not assigned A`Dimension then
        nf:=Components(A);
        require HasBaseField(A) : "The number fields of A shoud all be defined over the same base ring.";
        A`Dimension:=&+[Degree(E) : E in nf];
    end if;
    return A`Dimension;
end intrinsic;

intrinsic AbsoluteDimension(A::AlgEt)->RngInt
{Dimension of A over the prime field}    
    if not assigned A`AbsoluteDimension then
        A`AbsoluteDimension:=&+[AbsoluteDegree(E) : E in Components(A)];
    end if;
    return A`AbsoluteDimension;
end intrinsic;

intrinsic HasBaseField(A::AlgEt) -> BoolElt,FldNum
{Returns whether A has common base field. If this is the case it returns it.}
    if not assigned A`HasBaseField then
        nf,embs:=Components(A);
        F:=BaseRing(nf[1]);
        A`HasBaseField:=forall{ E : E in nf[2..#nf] | BaseRing(E) eq F };
        if A`HasBaseField then
            diag:=map< F->A | x:-> A!<nf[i]!x : i in [1..#nf]> >;
            A`BaseField:=<F,diag>;
        end if;
    end if;
    return A`HasBaseField;
end intrinsic;

intrinsic BaseField(A::AlgEt) -> FldNum
{Returns the common base field of the Algebra, if it exists.}
    if not assigned A`BaseField then
        require HasBaseField(A) : "The number fields should all be defined over the same Base ring/field.";
        // if HasBaseField is true, then it is assiged
    end if;
    return Explode(A`BaseField);
end intrinsic;

intrinsic PrimeField(A::AlgEt) -> FldNum
{Returns the prime field of the Algebra}
    if not assigned A`PrimeField then
        nf:=Components(A);
        F:=PrimeField(nf[1]); //this should always be Rationals()
        A`PrimeField:=F;
    end if;
    return A`PrimeField;
end intrinsic;

//------------
// Equality
//------------

intrinsic 'eq'(A1::AlgEt,A2::AlgEt) -> BoolElt
{A1 eq A2}
   return A1 cmpeq A2;
end intrinsic;

/*

// MOVE OR REMOVE THE NEXT
intrinsic HomsToC(A::AlgEt : Precision:=30)->SeqEnum[Map]
{returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30}
    CC:=ComplexField(Precision);
    images:=function(x)
        return &cat[[CC ! z : z in Conjugates(y : Precision:=Precision)] :y in Components(x)];
    end function;
    maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Dimension(A)] >;
    f:=&*[DefiningPolynomial(L[1]) : L in A`Components];
    assert &and [ Abs(Evaluate(f,g(PrimitiveElement(A)))) lt 10^-10 : g in maps]; //the precision here is quite arbitrary...
    return maps;
end intrinsic;

*/

/* TEST

    Attach("~/packages_github/AlgEt/AlgEt.m");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    assert #Basis(A) eq Dimension(A);

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    A;
    F,mFA:=BaseField(A);
    [ mFA(b) : b in Basis(MaximalOrder(F)) ];

    E:=NumberField(p);
    seq:=[E,K];
    A:=EtaleAlgebra(seq);
    HasBaseField(A);
    AbsoluteDimension(A);
    PrimeField(A);

    A:=EtaleAlgebra([E,E]);
    Dimension(A);
    AbsoluteDimension(A);
    BaseField(A);
    PrimeField(A);

*/
