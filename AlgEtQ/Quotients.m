/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

import "Ord.m" : crQZ , crZQ , MatrixAtoQ , inclusion_matrix ;

declare verbose Quotients, 1;

declare attributes AlgEtQIdl : ResidueField,PrimitiveElementResidueField;

///# Quotients
/// Let $I$ and $J$ be orders or fractional ideals such that $I\subseteq J$.
/// Then the quotient $J/I$ is finite.
/// We represent this quotient as a finite abelian group, a finite field, or a finite dimensional vector space over a finite field, depending on the properties of $I$ and $J$.

//------------
// Quotients
//------------

// The maps returned below use the Z-basis zb of I (an ideal or an order) at creation, and Minv, the inverse of the
// matrix of zb (the inclusion matrix of I), so that they stay correct if ZBasis(I) is replaced later, for example by
// ZBasisLLL(I). coords_wrt returns the coordinates of the elements of seq with respect to zb.
coords_wrt:=function(seq, Minv)
    return [ Eltseq(r) : r in Rows(MatrixAtoQ(seq)*Minv) ];
end function;

/// Given a fractional ideal $I$ and the $\mathbb{Z}$-basis of an ideal or order $J$ such that $J \subseteq I$, returns the abelian group $Q=I/J$ together with the quotient map $q:I\to Q$.
intrinsic Quotient(I::AlgEtQIdl, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map
{Given a fractional ideal I and the ZBasis of an ideal or order J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->Q.}
    // if J is not inside I, an error occurs while forming Q. so no need to check in advance
    A:=Algebra(I);
    require forall{ z : z in zbJ | Algebra(z) cmpeq A } : "the algebra is not the same";
    Minv,zbI:=inclusion_matrix(I);
    N := #zbI;
    F := FreeAbelianGroup(N);
    rel := [F ! cc : cc in coords_wrt(zbJ,Minv)];
    //mFI := map<F->A| x:->&+[Eltseq(x)[i]*zbI[i] : i in [1..N]]>;
    mFI := map<F->A| x:->DotProduct(Eltseq(x),zbI)>;
    mIF := map<A->F| x:-> F ! coords_wrt([x],Minv)[1]>;
    Q,qFQ := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to I/J.
    q:=map< A->Q | x:->qFQ(mIF(x)) , y:-> mFI(y@@qFQ) >;
    return Q,q;
end intrinsic;

/// Given fractional ideals $I$ and $J$ such that $J \subseteq I$, returns the abelian group $Q=I/J$ together with the quotient map $q:I \to Q$.
intrinsic Quotient(I::AlgEtQIdl, J::AlgEtQIdl) -> GrpAb, Map
{Given fractional ideals I and J such that J subset I, returns the abelian group Q=I/J together with the quotient map q:I->Q.}
    return Quotient(I,ZBasis(J));
end intrinsic;

/// Given an order $S$ and the $\mathbb{Z}$-basis of an ideal or order $J$ such that $J \subseteq S$, returns the abelian group $Q=S/J$ together with the quotient map $q:S\to Q$.
intrinsic Quotient(S::AlgEtQOrd, zbJ::SeqEnum[AlgEtQElt]) -> GrpAb, Map
{Given an order S and the ZBasis of an ideal or order J such that J subset S, returns the abelian group Q=S/J together with the quotient map q:S->Q.}
    // if J is not inside S, an error occurs while forming Q. so no need to check in advance
    A:=Algebra(S);
    require forall{ z : z in zbJ | Algebra(z) cmpeq A } : "the algebra is not the same";
    Minv,zbS:=inclusion_matrix(S);
    N := #zbS;
    F := FreeAbelianGroup(N);
    rel := [F ! cc : cc in coords_wrt(zbJ,Minv)];
    //mFS := map<F->A| x:->&+[Eltseq(x)[i]*zbS[i] : i in [1..N]]>;
    mFS := map<F->A| x:->DotProduct(Eltseq(x),zbS)>;
    mSF := map<A->F| x:-> F ! coords_wrt([x],Minv)[1]>;
    Q,qFQ := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to S/J.
    q:=map< A->Q | x:->qFQ(mSF(x)) , y:-> mFS(y@@qFQ) >;
    return Q,q;
end intrinsic;

/// Given an order $S$ and integral fractional $S$-ideal $I$, returns the abelian group $S/I$ and the quotient map $q:S \to S/I$ (with preimages). Important: the domain of $q$ is the parent algebra of $S$, since the elements of $S$ are expressed as elements of $A$. We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.
intrinsic ResidueRing(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb , Map
{Given an order S and integral fractional S-ideal I, returns the abelian group S/I and the quotient map q:S -> S/I (with preimages). Important: the domain of q is the Algebra of S, since the elements of S are expressed as elements of A. We stress that the output is a group and does not have a multiplication. This can be obtained by first taking preimages, doing the multiplication, and then applying the projection.}
    require Order(I) eq S : "wrong order";
    return Quotient(S,ZBasis(I));
end intrinsic;

/// Given a prime $P$ of $S$, returns a finite field $F$ isomorphic to $S/P$ and a surjection (with inverse) $S \to F$.
intrinsic ResidueField(P::AlgEtQIdl) -> FldFin, Map
{Given a prime P of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.}
    if not assigned P`ResidueField then
        assert IsPrime(P);
        S := Order(P);
        A := Algebra(S);
        Q,q := ResidueRing(S,P); // q:A->Q, where Q=S/P is an abelian group isomorphic to GF(p)^d
        _,p,d := IsPrimePower(#Q);
        Fp := GF(p);
        // We look for a in S whose image generates S/P over GF(p), that is, such that the images of
        // 1,a,...,a^(d-1) are linearly independent over GF(p). Most elements of S/P have this property.
        // We try the Z-basis of S first: its elements are small, which keeps the preimages under the output map small.
        // If needed, we try all the elements of S/P.
        powers := function(a) // returns 1,a,...,a^d, and their images in S/P = GF(p)^d as the rows of a matrix
            pows := [ One(A) ];
            for i in [1..d] do
                Append(~pows,pows[i]*a);
            end for;
            return Matrix(Fp,[ Eltseq(q(z)) : z in pows ]), pows;
        end function;
        found := false;
        for b in ZBasis(S) do
            M,pows := powers(b);
            if Rank(RowSubmatrix(M,1,d)) eq d then
                found := true;
                break b;
            end if;
        end for;
        if not found then
            for y in Q do
                M,pows := powers(y@@q);
                if Rank(RowSubmatrix(M,1,d)) eq d then
                    found := true;
                    break y;
                end if;
            end for;
            // Here a and its powers can be large: we replace them with small elements with the same images in S/P.
            pows := [ ShortElementModulo(z,P) : z in pows ];
        end if;
        assert found;
        B := RowSubmatrix(M,1,d);
        c := Solution(B,M[d+1]); // a^d = c[1] + c[2]*a + ... + c[d]*a^(d-1) in S/P
        g := Polynomial(Fp,[ -c[i] : i in [1..d] ] cat [ 1 ]); // the minimal polynomial of the image of a
        // We send a to a root r of g in F: this gives an isomorphism S/P -> F.
        F := FiniteField(#Q);
        r := Roots(PolynomialRing(F)!g)[1][1];
        R := Matrix(Fp,[ Eltseq(r^i) : i in [0..d-1] ]); // the coordinates of 1,r,...,r^(d-1) in F over GF(p)
        BR := B^-1*R; // from the coordinates in S/P = GF(p)^d to the coordinates in F
        Rinv := R^-1; // from the coordinates in F to the coordinates with respect to 1,r,...,r^(d-1)
        toF := function(x)
            return F!Eltseq(Vector(Fp,Eltseq(q(x)))*BR);
        end function;
        fromF := function(z)
            // z = e[1] + e[2]*r + ... + e[d]*r^(d-1) is lifted to e[1] + e[2]*a + ... + e[d]*a^(d-1), with 0 <= e[i] < p
            e := Eltseq(Vector(Fp,Eltseq(z))*Rinv);
            return &+[ (Integers()!e[i])*pows[i] : i in [1..d] ];
        end function;
        P`ResidueField:=<F, map< A -> F | x:->toF(x), z:->fromF(z) >>;
    end if;
    return Explode(P`ResidueField);
end intrinsic;

/// Returns an element of the given prime $P$ of $S$ that maps to the primitive element of the residue field $S/P$, that is, a generator of the multiplicative group $(S/P)^*$.
intrinsic PrimitiveElementResidueField(P::AlgEtQIdl)->AlgEtQElt
{Returns an element of the prime P of S that maps to the primitive element of the residue field S/P, that is, a generator of the multiplicative group (S/P)^*.}
    if not assigned P`PrimitiveElementResidueField then
        F,f:=ResidueField(P);
        P`PrimitiveElementResidueField:=PrimitiveElement(F)@@f;
    end if;
    return P`PrimitiveElementResidueField;
end intrinsic;

/// Let $I$, $J$ be orders or fractional ideals, $P$ a fractional $R$-ideals such that:
/// - $P$ is prime of some order $R$, with residue field $K$;
/// - $J \subseteq I$ and $I/J$ is a vector space $V$ over $K$, say of dimension $d$.
/// The intrinsic returns the $K$-vector space $K^d=V$ and the natural surjection $I \to V$ (with preimages).
intrinsic QuotientVS(I::AlgEtQOrd, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map
{Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of some order R, with residue field K;
 - J subset I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).}
    S := Order(P);
    return QuotientVS(S!!OneIdeal(I),S!!OneIdeal(J),P);
end intrinsic;

///ditto
intrinsic QuotientVS(I::AlgEtQOrd, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map
{Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).}
    S := Order(P);
    return QuotientVS(S!!OneIdeal(I),S!!J,P);
end intrinsic;

///ditto
intrinsic QuotientVS(I::AlgEtQIdl, J::AlgEtQOrd, P::AlgEtQIdl) -> ModRng, Map
{Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with preimages).}
    S := Order(P);
    return QuotientVS(S!!I,S!!OneIdeal(J),P);
end intrinsic;

///ditto
intrinsic QuotientVS(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map
{Let I, J, P be fractional R-ideals such that:
 - P is prime of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with preimages).}
    S := Order(P);
    require Order(I) eq S and Order(J) eq S : "the ideals must be over the same order ";
    require J subset I : "The second argument should be a subset of the first.";
    assert2 P*I subset J;
    K,k:=ResidueField(P);
    A := Algebra(S);
    d := Ilog(#K,Integers() ! (Index(J)/Index(I))); // d = dim(I/J) over (S/P)
    V := KModule(K,d);
    //need to find a basis of I/J over R/P.
    Minv,zbI := inclusion_matrix(I); // see coords_wrt
    N := #zbI;
    F := FreeAbelianGroup(N);
    relJ := [F ! cc : cc in coords_wrt(ZBasis(J),Minv)];
    rel:=relJ;
    //mFI := map<F->A| x:->&+[Eltseq(x)[i]*zbI[i] : i in [1..N]]>;
    mFI := map<F->A| x:->DotProduct(Eltseq(x),zbI)>;
    mIF := map<A->F| x:-> F ! coords_wrt([x],Minv)[1]>;
    Q,q := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to I/J.
    bas := [];
    for i in [1..d] do
    //for each iteration of the loop we mod-out from I the S-ideal generated by J and the already found elements of the basis of I/J over S/P
        elt_F := (Q.1@@q);
        elt_I := mFI(elt_F);
        Append(~bas,elt_I);
        rel_i:=[mIF(elt_I*z) : z in ZBasis(S)]; // a Z-basis of elt_I*S (elt_I might be a zero-divisor)
        rel := rel cat rel_i;
        Q, q := quo<F|rel>; //q:F->Q
    end for;
    assert IsTrivial(Q);
    //function mIV using HNF: bas[i]*S+J eq I. exploit this on ZBasis level to find the S-coordinates of ZBasis(I) wrt bas[i]'s
    zbJ:=ZBasis(J);
    zbS:=ZBasis(S);
    gens:=&cat[[ b*z: z in zbS ] : b in bas];
    mat:=MatrixAtoQ(gens cat zbJ);
    den:=Denominator(mat);
    dmat:=crQZ(den*mat);
    H,Tr:=HermiteForm(dmat);//Tr*M eq H
    HI,TrI:=HermiteForm(crQZ(den*MatrixAtoQ(zbI)));
    TrI:=crZQ(TrI);
    Tr1:=crZQ(Matrix(Rows(Tr)[1..#zbI]));
    C:=TrI^-1*Tr1;
    new_coords_zbI:=[];
    for k in [1..#zbI] do
        zbIk:=[];
        for i in [0..d-1] do
            //coord_i:=&+[C[k,i*N+j]*zbS[j] : j in [1..#zbS]];
            coord_i:=DotProduct([C[k,i*N+j] : j in [1..#zbS]],zbS);
            Append(~zbIk,coord_i);
        end for;
        Append(~new_coords_zbI,zbIk);
    end for;
    mIV:=function(x)
        xinI:=coords_wrt([x],Minv)[1];
        //coords_inS:=[ &+[xinI[i]*new_coords_zbI[i][k]: i in [1..#zbI]]  : k in [1..d]];
        coords_inS:=[ DotProduct(xinI,[new_coords_zbI[i][k]: i in [1..#zbI]])  : k in [1..d]];
        coords_inK:=[k(c) : c in coords_inS];
        return &+[coords_inK[i]*V.i : i in [1..d]];
    end function;
    mVI := function(y)
        y:=[Eltseq(y)[j]@@k : j in [1..d]];
        return DotProduct(bas,y);
        //return &+[ bas[j]*(Eltseq(y)[j]@@k) : j in [1..d] ];
    end function;
    return V, map<A->V | x:->mIV(x), y:->mVI(y) >;
end intrinsic;

/* TESTS

    printf "### Testing Quotients:";
    //AttachSpec("../spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    icm:=ICM(E);
    pp:=SingularPrimes(E);
    for P in pp do
        r:=ResidueField(P);
        _:=PrimitiveElementResidueField(P);
        assert #r eq #Quotient(OneIdeal(E),P);
        assert #r eq #ResidueRing(E,P);
        assert #r eq #QuotientVS(OneIdeal(E),P,P);
        printf ".";
    end for;

    for I in icm do
        _:=Quotient(OneIdeal(E),MakeIntegral(I));
        for P in pp do
            PI:=P*I;
            n:=#Quotient(I,PI);
            assert n eq #QuotientVS(I,PI,P);
            printf ".";
        end for;
    end for;

    // ResidueField for residue fields of size 13^4 and 2^4: 13 and 2 are inert in Q(zeta_5).
    // The map O->F is a ring homomorphism, with kernel P and with preimages.
    K:=EtaleAlgebra(x^4+x^3+x^2+x+1);
    O:=MaximalOrder(K);
    for p in [13,2] do
        P:=p*OneIdeal(O);
        F,f:=ResidueField(P);
        assert #F eq p^4;
        assert f(One(K)) eq One(F);
        assert forall{ b : b in ZBasis(P) | f(b) eq Zero(F) };
        for i in [1..10] do
            a:=Random(O : ZeroDivisorsAllowed:=true);
            b:=Random(O : ZeroDivisorsAllowed:=true);
            assert f(a*b) eq f(a)*f(b) and f(a+b) eq f(a)+f(b);
            z:=Random(F);
            assert (z@@f) in O and f(z@@f) eq z;
        end for;
        assert IsPrimitive(f(PrimitiveElementResidueField(P)));
        printf ".";
    end for;

    // the maps of Quotient and QuotientVS stay correct when the stored Z-basis of I is replaced later, as ZBasisLLL(I)
    // does when the basis is not LLL-reduced
    A:=EtaleAlgebra(x^4+16);
    E:=EquationOrder(A);
    I:=Ideal(E,[Random(E,100) : i in [1..3]]);
    P:=PrimesAbove(Ideal(E,[A!3]))[1];
    J:=I*P;
    Q,q:=Quotient(I,J);
    V,v:=QuotientVS(I,J,P);
    zb:=ZBasis(I);
    I`ZBasis:=[ zb[1] ] cat [ zb[i]+(i+1)*zb[1] : i in [2..#zb] ]; // another Z-basis of I
    delete I`inclusion_matrix; // as ZBasisLLL does
    assert forall{ y : y in Q | q(y@@q) eq y } and forall{ z : z in ZBasis(J) | q(z) eq Zero(Q) };
    assert forall{ y : y in V | v(y@@v) eq y } and forall{ z : z in ZBasis(J) | v(z) eq Zero(V) };
    printf ".";

    SetAssertions(1);
    printf " all good!";
*/