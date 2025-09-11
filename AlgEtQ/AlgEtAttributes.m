/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;


//------------
// Printing
//------------

///hide-all
intrinsic Print(A::AlgEtQ)
{Prints the defining polynomial or the components defining A.}
    if assigned A`DefiningPolynomial then
        f:=DefiningPolynomial(A);
        printf "Etale Algebra over (%o) defined by %o", FieldOfFractions(CoefficientRing(f)),f;
    else
        printf "Etale Algebra product of %o", Components(A);
    end if;
end intrinsic;
///hide-none

//------------
// Access attributes
//------------

///## Attributes

///### Components, equality testing and defining polynomial
/// Two étale algebras are defined to be equal if the ordered sequence of their components are the same.

/// Returns the number fields of which $A$ is a product of, together with embeddings and projections.
intrinsic Components(A::AlgEtQ) -> SeqEnum,SeqEnum,SeqEnum
{Returns the sequence of number fields of which A is a product of, together with embeddings and projections.}
    return Explode(A`Components);
end intrinsic;

intrinsic 'eq'(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt
{Two étale algebras are defined to be equal if they have the same components.}
   c1,e1,p1:=Components(A1);
   c2,e2,p2:=Components(A2);
   return <c1,e1,p1> eq <c2,e2,p2>;
end intrinsic;

/// Returns the defining polynomial of $A$, if the components are distinct number fields.
intrinsic DefiningPolynomial(A::AlgEtQ) -> RngUPolElt
{Returns the defining polynomial of A, if the components are distinct number fields.}
    if not assigned A`DefiningPolynomial then
        polys:=[DefiningPolynomial(L) : L in Components(A)];
        require #polys eq #Seqset(polys) : "The number fields defining A are not distinct.";
        A`DefiningPolynomial:=&*(polys);
    end if;
    return A`DefiningPolynomial;
end intrinsic;

///### Base field and prime field

/// Returns whether the components of $A$ all have the same base field.
intrinsic HasBaseField(A::AlgEtQ) -> BoolElt
{Returns whether the components of A all have the same base field.}
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

/// Returns whether the common base field of the components of $A$ if it exists, together with the diagonal embedding.
intrinsic BaseField(A::AlgEtQ) -> FldNum,Map
{Returns the common base field of the components of the étale algebra, if it exists, together with the diagonal embedding.}
    if not assigned A`BaseField then
        require HasBaseField(A) : "The number fields should all be defined over the same Base ring/field.";
        // if HasBaseField is true, then it is assigned
    end if;
    return Explode(A`BaseField);
end intrinsic;

intrinsic PrimeField(A::AlgEtQ) -> FldNum
{Returns the prime field of the étale algebra.}
    if not assigned A`PrimeField then
        nf:=Components(A);
        F:=PrimeField(nf[1]); //this should always be Rationals()
        A`PrimeField:=F;
    end if;
    return A`PrimeField;
end intrinsic;

///### Dimension

/// Returns the dimension of A over the base field, which in this case is $\mathbb{Q}$.
intrinsic Dimension(A::AlgEtQ)->RngInt
{Returns the dimension of A over the base field, which in this case is the rational field.}
    if not assigned A`Dimension then
        nf:=Components(A);
        require HasBaseField(A) : "The number fields of A should all be defined over the same base ring.";
        A`Dimension:=&+[Degree(E) : E in nf];
    end if;
    return A`Dimension;
end intrinsic;

/// Returns the dimension of $A$ over the prime field.
intrinsic AbsoluteDimension(A::AlgEtQ)->RngInt
{Returns the dimension of A over the prime field, which in this case is the rational field.}
    if not assigned A`AbsoluteDimension then
        A`AbsoluteDimension:=&+[AbsoluteDegree(E) : E in Components(A)];
    end if;
    return A`AbsoluteDimension;
end intrinsic;

/* TESTS

    printf "### Testing Attributes and Equality:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    B:=EtaleAlgebra(f);
    assert not A eq B;
    assert A ne B;
    assert A eq A;
    assert A cmpeq A;
    assert not A cmpeq B;
    SetAssertions(1);
    printf " all good!";

*/