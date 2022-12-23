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

intrinsic TwoGeneratingSet(I::AlgEtIdl)
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

/* TEST

	AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/MinimalGenerators.m");
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    time P,p:=PicardGroup(E : GRH:=true);

    t0:=Cputime();
    for g in Generators(P) do 
        I:=p(g);
    end for;
    Cputime(t0);
    t0:=Cputime();
    for g in Generators(P) do 
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
    end for;
    Cputime(t0);

    t0:=Cputime();
    SetProfile(true);
    for g in Generators(P) do 
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
    end for;
    SetProfile(false);
    G:=ProfileGraph();
    ProfilePrintByTotalTime(G);
    Cputime(t0);

    t0:=Cputime();
    for g in P do 
        I:=p(g);
    end for;
    Cputime(t0);

    t0:=Cputime();
    for g in P do 
        I:=p(g);
        TwoGeneratingSet(~I);
        assert #Genertors(I) le 2;
    end for;
    Cputime(t0);
*/
