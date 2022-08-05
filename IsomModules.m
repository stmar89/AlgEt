/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose IsomModules, 2;

/* TODO
    - add the map
    - implement the faster method from Bley-Hofmann-Johnston: https://arxiv.org/abs/2202.03526
*/

//------------
// Isomorphism Testing for modules
//------------

intrinsic IsIsomorphic(I::AlgEtMod,J::AlgEtMod) -> BoolElt
{Given two modules I and J returns wheater they are isomorphic.}
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(I);
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
   
    require EquationOrder(Algebra(S)) subset S : "Implemented only for modules over orders containing the equation order";
    pi:=m(PrimitiveElement(Algebra(S)));
    // based on the following, two Z[pi] modules are isomorphic iff the matrices representing multiplcition by pi are Z-conjugte 
    matI:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I))),Integers());
    matJ:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(J)],ZBasis(J))),Integers());
    test:=AreGLConjugate(matI,matJ);
    return test;
end intrinsic;

/* TEST
   
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
    _<x>:=PolynomialRing(Integers()); 
    g:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    K:=EtaleAlgebra(g);
    nf:=NumberFields(K);
    q:=2;
    pi:=PrimitiveElement(K);
    R:=Order([pi,q/pi]);
    O:=MaximalOrder(K);
    assert IsBass(R);
    assert #FindOverOrders(R) eq 2;
    assert #PicardGroup(R) eq 3;
    assert #PicardGroup(O) eq 1;

    V:=EtaleAlgebra(&cat[nf : i in [1..3]]);
    Vnf:=NumberFields(V);
    m:=NaturalAction(K,V);
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    // the following works for this exmaple, not in general.
    Mff:=Module(R,m,<ff_prod[1] : i in [1..3]>);
    gensMff:=Generators(Mff);
    //time candidates:=IntermediateModules(MO,Mff);
    Q,q:=Quotient(MO,Mff); Q;
    subs:=Subgroups(Q);
    time candidates:=[ Module(R,m,[ g@@q : g in Generators(H`subgroup)] cat gensMff) : H in subs ];
    time candidates:=Seqset(candidates);
    candidates:=Setseq(candidates);

    classes:=[];
    for I in candidates do
        printf ".";
        if not exists{ J : J in classes | IsIsomorphic(I,J) } then
            printf "!";
            Append(~classes,I);
        end if;
    end for;
    assert #classes eq 6; //This is Example 6.1 in "Computing abelian varieties over finite fields isogenous to a power", by Marseglia






*/

