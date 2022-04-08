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
{A procedure that given an invertible ideal I creates produces 2 elemetns that generate I. If I is known to be principal nothing is done. }
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        if assigned I`IsIntegral and IsIntegral(I) then
            a:=Algebra(I)!MinimalInteger(I);
        else
            a:=Random(I);
        end if;
        S:=Order(I);
        Q,q:=quo_frac_idl(I,[a*z : z in ZBasis(S)]);
        if #Q eq 0 then
            I`Generators:=[a]; //the ideal is principalA
        else
            repeat
                repeat
                    b:=Random(Q);
                until b ne Zero(Q);
                b:=b@@q;
            until I eq Ideal(S,[a,b]);
        end if;
        I`Generators:=[a,b];
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
