/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

declare verbose Quotients, 1;

declare attributes AlgEtIdl : ResidueField,PrimitiveElementResidueField;

/*TODO:

*/


//------------
// Quotients
//------------

intrinsic Quotient(I::AlgEtIdl, zbJ::SeqEnum[AlgEtElt]) -> GrpAb, Map
{ Given an ideal I and the ZBasis of an ideal J such that  J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J. J can also be an order. } 
    // if J is not inside I, an error occurs while forming Q. so no need to check in advance
    A:=Algebra(I);
    zbI:=ZBasis(I);
	N := #zbI;
	F := FreeAbelianGroup(N);
	rel := [F ! cc : cc in AbsoluteCoordinates(zbJ,I)];
	//mFI := map<F->A| x:->&+[Eltseq(x)[i]*zbI[i] : i in [1..N]]>;
	mFI := map<F->A| x:->SumOfProducts(Eltseq(x),zbI)>;
	mIF := map<A->F| x:-> F ! AbsoluteCoordinates([x],I)[1]>;
	Q,qFQ := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to I/J.
    q:=map< A->Q | x:->qFQ(mIF(x)) , y:-> mFI(y@@qFQ) >; 
    return Q,q;
end intrinsic;

intrinsic Quotient(I::AlgEtIdl, J::AlgEtIdl) -> GrpAb, Map
{ given fractional ideals J subset I, returns the abelian group Q=I/J together with the quotient map q:I->J } 
    return Quotient(I,ZBasis(J));
end intrinsic;

intrinsic Quotient(S::AlgEtOrd, zbJ::SeqEnum[AlgEtElt]) -> GrpAb, Map
{ Given an order S and the ZBasis of an ideal J such that  J subset S, returns the abelian group Q=S/J together with the quotient map q:S->J. J can also be an order. } 
    // if J is not inside S, an error occurs while forming Q. so no need to check in advance
    A:=Algebra(S);
    zbS:=ZBasis(S);
	N := #zbS;
	F := FreeAbelianGroup(N);
	rel := [F ! cc : cc in AbsoluteCoordinates(zbJ,S)]; // this absolute coordinates uses inclusion_matrix. fast!
	//mFS := map<F->A| x:->&+[Eltseq(x)[i]*zbS[i] : i in [1..N]]>;
	mFS := map<F->A| x:->SumOfProducts(Eltseq(x),zbS)>;
	mSF := map<A->F| x:-> F ! AbsoluteCoordinates([x],S)[1]>;
	Q,qFQ := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to S/J.
    q:=map< A->Q | x:->qFQ(mSF(x)) , y:-> mFS(y@@qFQ) >; 
    return Q,q;
end intrinsic;

intrinsic ResidueRing(S::AlgEtOrd,I::AlgEtIdl) -> GrpAb , Map
{given an integral ideal I of S, returns the abelian group S/I and the epimorphism pi:S -> S/I (with inverse map). Important: the domain of pi is the Algebra of S, since the elements of S are usually expressed al elements of A. For eg Parent(Random(S)) = Algebra(S)}
    require Order(I) eq S : "wrong order";
    return Quotient(S,ZBasis(I));
//  OLD
//     require IsIntegral(I): "I must be an integral ideal of S";
//     A:=Algebra(S);
//     N:=AbsoluteDimension(A);
//     F:=FreeAbelianGroup(N);
//     S_to_F:=function(x0)
//         assert Parent(x0) eq A;
//         x_inS:=AbsoluteCoordinates([x0],S)[1];
//         return (F ! Eltseq(x_inS)) ;
//     end function;
//     F_to_S:=function(y)
//         y_inA:=&+[ZBasis(S)[i]*Eltseq(y)[i] : i in [1..N]];
//         return y_inA;
//     end function;
//     StoF:=map< A -> F | x :-> S_to_F(x), y :-> F_to_S(y)>;
//     rel:=[F ! x : x in AbsoluteCoordinates(ZBasis(I),S)];
//     Q,q:=quo<F|rel>; //Q=S/I
//     m:=StoF*q; //m is a map from S to Q
//     assert #Q eq Index(S,I);
//     assert2 forall{x : x in ZBasis(I) | m(x) eq Zero(Q)};
//     assert2 forall{x : x in ZBasis(S) | ((m(x))@@m - x) in I};
//     return Q,m;
end intrinsic;

intrinsic ResidueField(P::AlgEtIdl) -> FldFin, Map
{ given P a prime of S, returns a finite field F isomorphic to S/P and a surjection (with inverse) S->F.}
    if not assigned P`ResidueField then
        assert IsPrime(P);
        S := Order(P);
        Q,q := ResidueRing(S,P); //q:S->S/P
        size := #Q;
        F := FiniteField(size);
        min_poly := PolynomialRing(Integers())!DefiningPolynomial(F);
        //the following loop is naive
        for y in Q do
            if q(Evaluate(min_poly,y@@q)) eq Zero(Q) then
                prim_elt_inA := y@@q;
                break y;
            end if;
        end for;
        assert assigned prim_elt_inA;
        //now I need to build  the map
        G,gGF := AdditiveGroup(F); //g:G->F
        hGQ := iso<G->Q | [q(prim_elt_inA^i) : i in [0..Degree(min_poly)-1]]>;
        hQG := Inverse(hGQ);
        map := q*hQG*gGF;
	    P`ResidueField:=<F, map>;
    end if;
    return Explode(P`ResidueField);
end intrinsic;

intrinsic PrimitiveElementResidueField(P::AlgEtIdl)->AlgEtElt
{Returns an element of P that maps to the primitive element of the residue field S/P, that is a multiplicative generator of (S/P)^*.}
    if not assigned P`PrimitiveElementResidueField then
        F,f:=ResidueField(P);
        P`PrimitiveElementResidueField:=PrimitiveElement(F)@@f;
    end if;
    return P`PrimitiveElementResidueField;
end intrinsic;

/* The following version is promising, but still bugged. One assertions fails.
intrinsic QuotientVS(I::AlgEtIdl, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
{Let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
    require J subset I : "Teh second argument should be a subset of the first.";
	S := Order(P);
    assert2 P*(Ideal(S,ZBasis(I))) subset Ideal(S,ZBasis(J));
    assert2 S subset MultiplicatorRing(I);
    assert2 S subset MultiplicatorRing(J);
    K,k:=ResidueField(P);
	A := Algebra(S);
    Q,q:=Quotient(I,J);
	d := Ilog(#K,#Q); // d = dim(I/J) over (S/P)
    // SOMETHING WRONG IN HERE. I Think that there are issues when K is not a prime field
    matrices:=[Matrix(K,[ Eltseq(Q!(q(zS*(Q.j@@q)))) : j in [1..Ngens(Q)]]) : zS in ZBasis(S) ];
    matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
    V:=RModule(matrices);
    assert Dimension(V) eq d;
    assert #Basis(V) eq Ngens(Q);
    bijVQ:=function(y)
        eltseq:=Eltseq(y);
        N:=#Eltseq(eltseq[1]);
        // y is represented in V wrt to generators that are in bijection with Q. 
        // So even if a priori the coeficients of each coordinates of y are in the finite field K, 
        // they are actually all in the prime field of K, and hence can be coerced into integers.
        assert forall{i : i in [1..Ngens(Q)] | IsZero(Eltseq(eltseq[i])[2..N]) };
        out:=&+[Q.i*(Integers()!Eltseq(eltseq[i])[1]) : i in [1..Ngens(Q)]];
        return out; 
    end function;
    bij:=map< Q->V | x :-> &+[V.i*Eltseq(x)[i] : i in [1..Ngens(Q)]],
                     y :-> bijVQ(y) >;
    v:=map< A-> V | x:->bij(q(x)) , y:->(y@@bij)@@q >;
    return V,v;
end intrinsic;
*/

intrinsic QuotientVS(I::AlgEtOrd, J::AlgEtOrd, P::AlgEtIdl) -> ModRng, Map
{Let I, J be orders, P a fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
	S := Order(P);
    return QuotientVS(S!!OneIdeal(I),S!!OneIdeal(J),P);
end intrinsic;

intrinsic QuotientVS(I::AlgEtOrd, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
{Let I be an order, J and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
	S := Order(P);
    return QuotientVS(S!!OneIdeal(I),S!!J,P);
end intrinsic;

intrinsic QuotientVS(I::AlgEtIdl, J::AlgEtOrd, P::AlgEtIdl) -> ModRng, Map
{Let J be an order, I and  P be fractional R-ideals such that:
 - P is prime of of some order R, with residue field K;
 - J in I and I/J is a vector space V over K, say of dimension d.
 The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
	S := Order(P);
    return QuotientVS(S!!I,S!!OneIdeal(J),P);
end intrinsic;

intrinsic QuotientVS(I::AlgEtIdl, J::AlgEtIdl, P::AlgEtIdl) -> ModRng, Map
{
 let I, J, P be fractional R-ideals such that:
 - P is prime of of some order R;
 - J in I and I/J is a vector space over R/P, say of dimension d;
 the function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
	S := Order(P);
    require Order(I) eq S and Order(J) eq S : "the ideals must be over the same order ";
    require J subset I : "Teh second argument should be a subset of the first.";
    assert2 P*I subset J;
    K,k:=ResidueField(P);
	A := Algebra(S);
	d := Ilog(#K,Integers() ! (Index(J)/Index(I))); // d = dim(I/J) over (S/P)
	V := KModule(K,d);
	//need to find a basis of I/J over R/P.
	zbI := ZBasis(I);
	N := #zbI;
	F := FreeAbelianGroup(N);
	relJ := [F ! cc : cc in AbsoluteCoordinates(ZBasis(J),I)];
    rel:=relJ;
	//mFI := map<F->A| x:->&+[Eltseq(x)[i]*zbI[i] : i in [1..N]]>;
	mFI := map<F->A| x:->SumOfProducts(Eltseq(x),zbI)>;
	mIF := map<A->F| x:-> F ! AbsoluteCoordinates([x],I)[1]>;
	Q,q := quo<F|rel>; //q:F->Q. Q is an "abstract" abelian group isomorphic to I/J.
	bas := [];
	for i in [1..d] do
    //for each iteration of the loop we mod-out from I the S-ideal generated by J and the already found elements of the basis of I/J over S/P
		elt_F := (Q.1@@q);
		elt_I := mFI(elt_F);
		Append(~bas,elt_I);
        rel_i:=[mIF(bb) : bb in ZBasis(Ideal(S,elt_I))];
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
            coord_i:=SumOfProducts([C[k,i*N+j] : j in [1..#zbS]],zbS);
            Append(~zbIk,coord_i);
        end for;
        Append(~new_coords_zbI,zbIk);
    end for;
    mIV:=function(x)
        xinI:=AbsoluteCoordinates([x],I)[1];
        //coords_inS:=[ &+[xinI[i]*new_coords_zbI[i][k]: i in [1..#zbI]]  : k in [1..d]];
        coords_inS:=[ SumOfProducts(xinI,[new_coords_zbI[i][k]: i in [1..#zbI]])  : k in [1..d]];
        coords_inK:=[k(c) : c in coords_inS];
        return &+[coords_inK[i]*V.i : i in [1..d]];
    end function;
	mVI := function(y)
        y:=[Eltseq(y)[j]@@k : j in [1..d]];
		return SumOfProducts(bas,y);
		//return &+[ bas[j]*(Eltseq(y)[j]@@k) : j in [1..d] ];
    end function;
    return V, map<A->V | x:->mIV(x), y:->mVI(y) >;
end intrinsic;

/* TEST

*/
