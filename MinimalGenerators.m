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

quo_frac_idl:=function(I,zbJ)
// I->I/J
    N:=AbsoluteDimension(Algebra(I));
    F:=FreeAbelianGroup(N);
    zbI:=ZBasis(I);
    matI:=MatrixAtoQ(zbI);
    matJ:=MatrixAtoQ(zbJ);
    rel:=[F ! Eltseq(x) : x in Rows(matJ*matI^-1)];
    Q,q:=quo<F|rel>;
    qq:=map<Algebra(I) -> Q | x:->q(F!AbsoluteCoordinates([x],zbI)[1]) , y:->&+[Eltseq(y@@q)[i]*zbI[i] : i in [1..N]] >;
    return Q,qq;
end function;

intrinsic TwoGeneratingSet(I::AlgEtIdl)
{A procedure that given an invertible ideal I put in the attibute I`Generators two non-zerodivisors in I that generate I. If I is known to be principal, that is I`Generators consists of one single element, nothing is done.}
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        if assigned I`IsIntegral and IsIntegral(I) then
            a:=Algebra(I)!MinimalInteger(I);
        else
            a:=Random(I : ZeroDivisorsAllowed:=false );
        end if;
        S:=Order(I);
        Q,q:=quo_frac_idl(I,[a*z : z in ZBasis(S)]);
        if #Q eq 0 then
            I`Generators:=[a]; //the ideal is principal
        else
            repeat
                repeat
                    b:=Random(Q);
                until b ne Zero(Q);
                b:=b@@q;
                assert not IsZeroDivisor(b);
            until I eq Ideal(S,[a,b]);
        end if;
        I`Generators:=[a,b];
    end if;
end intrinsic;

// Lemma:
// Let I be an invertible fractional S-ideal generated globally by I=aS+bS.
// For every integer N, we have I^N=a^NS+b^NS.
// Proof: we will show it locally at every prime. Since I is invertible, then it is locally principal.
// Pick a prime P. Then I_P=(a/1)S_P or I_P=(b/1)S_P. wlog assume the first. Then I_P^N=(a^N/1)S_P. QED.
//
// Remark: this trick does not work for product of invertible ideals I and J. One can have I_P1=(a) J_P1=(c) but I_P2=(a), J_P2=(d) so one cannot avoid the mixed products.


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
