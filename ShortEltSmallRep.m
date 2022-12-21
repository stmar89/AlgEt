/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ShortEltSmallRep, 2;
import "Ord.m" : MatrixQtoA,MatrixAtoQ,MatrixAtoZ;

declare attributes AlgEtIdl : ShortestElement, SmallRepresentative;

//------------
// ShortestElement
//------------

intrinsic ShortestElement(I::AlgEtIdl) ->AlgEtElt
{Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by enumerating short vectors in I, and pick the first one which is a non-zerodivisor.}
    if not assigned I`ShortestElement then
        ZBasisLLL(I);
        L:=Lattice(MatrixAtoQ(ZBasis(I)));
        // b:=Basis(LLL(L)); //we reduce above, so it is cached.
        b:=Basis(L);
        b:=[ Norm(c) : c in b ];
        k:=0;
        stop:=false;
        min:=Min(b);
        repeat
            p:=ShortVectors(L,2^(-k)*min,min*2^k);
            for i in [1..#p] do
                elt:=MatrixQtoA(Algebra(I),Matrix([Eltseq(p[i][1])]))[1];
                if not IsZeroDivisor(elt) then
                    stop:=true;
                    break i;
                end if;
            end for;
            k+:=1;
        until stop;
        I`ShortestElement:= elt;
    end if;
    return I`ShortestElement;
end intrinsic;

//------------
// SmallRepresentative
//------------

intrinsic SmallRepresentative(I::AlgEtIdl) ->AlgEtIdl,AlgEtElt
{Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortestElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.}
    if not assigned I`SmallRepresentative then                                            
        R:=Order(I);
        cRI:=ColonIdeal(R,I);
        a:=ShortestElement(cRI);
        aI:=a*I;
        // the ZBasis of aI might be very big. We make it smaller.
        ZBasisLLL(aI);
        vprintf ShortEltSmallRep,2: "SmallRepresentative:\n
                                I = %o\n,aI = %o\n",PrintSeqAlgEtElt(ZBasis(I)),PrintSeqAlgEtElt(ZBasis(aI));
        I`SmallRepresentative:=<aI,a>;
    end if;
    return Explode(I`SmallRepresentative);
end intrinsic;


/* TEST
    
	AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    //f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    f:=(x^2+5)*(x^2+7)*(x^2+11);
    assert IsSquarefree(f);
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    ff:=Conductor(E);
    _:=ShortestElement(ff);
    oo:=FindOverOrders(E); 
    for S in oo do
        ff:=Conductor(S);
        _:=ShortestElement(ff);
    end for;

	AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-100*x^3-100*x^2-100*x-100;
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true); //this might take a while. timings are very inconsistent
    repeat
        Ii:=Random(P);
    until Ii ne Zero(P);
    I:=p(Ii);
    #Generators(I);

    myHash(I);
    aI:=SmallRepresentative(I);
    myHash(aI);
    time TwoGeneratingSet(I);
    time TwoGeneratingSet(aI);

    // is small representative, canonical? mmm confusing...

	AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-100*x^3-100*x^2-100*x-100;
    //f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    time P,p:=PicardGroup(E);
    I0:=p(Random(P));
    I:=SmallRepresentative(I0);
    for i in [1..10] do
        Index(E,I);
        I:=SmallRepresentative(I);
    end for;

*/
