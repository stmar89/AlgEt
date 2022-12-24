/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

//declare verbose ?????, 1;

/*TODO:

*/

//------------
// TwoGeneratingSet
//------------
import "Ord.m" : MatrixAtoQ,MatrixAtoZ;

intrinsic TwoGeneratingSet(I::AlgEtQIdl)
{A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.}
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        // if assigned I`IsIntegral and IsIntegral(I) then
        //     a:=MinimalInteger(I);
        // else
        //     a:=Random(I : ZeroDivisorsAllowed:=false );
        // end if;
        a:=ShortestElement(I);
        S:=Order(I);
        Q,q:=Quotient(I,[a*z : z in ZBasis(S)]);
        if IsTrivial(Q) then
            I`Generators:=[a]; //the ideal is principal
        else
            repeat
                repeat
                    b:=Random(Q);
                until b ne Zero(Q);
                b:=b@@q;
            //until I eq Ideal(S,[a,b]);
            until Q eq sub<Q|[q(b*z):z in ZBasis(S)]>;
            I`Generators:=[a,b];
        end if;
    end if;
end intrinsic;

/* TESTS

    printf "### TwoGenerators:";
	AttachSpec("~/packages_github/AlgEt/spec");
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    P,p:=PicardGroup(E : GRH:=true); //~10 secs

    t0:=Cputime();
    for g in Generators(P) do 
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
    end for;
    Cputime(t0);

    // test if TwoGeneratingSet makes the power faster
    // Conlcusion: yes. By quite a bit!
    f:=x^4-100*x^3-100*x^2-100*x-100;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    repeat
        Ii:=Random(P);
    until Ii ne Zero(P);
    I:=p(Ii);

    delete I`IsInvertible;
    exp:=[ Random(2,30) : i in [1..100]];
    l1:=[ I^i : i in exp ];

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    assert #Generators(I) eq 2;
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;

    I:=SmallRepresentative(I);
    delete I`IsInvertible;
    l1:=[ I^i : i in exp ];

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf " all good!\n"; 

*/
