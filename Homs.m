/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtHoms, 1;

declare attributes AlgEt : HomsToC;

//------------
// Homomorphism to the Compelx Numbers
//------------

intrinsic HomsToC(A::AlgEt : Precision:=30)->SeqEnum[Map]
{returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30}
    if not assigned A`HomsToC then
        require Type(PrimeField(A)) eq  FldRat : "The algebra needs to be over Q";
        CC:=ComplexField(Precision);
        images:=function(x)
            return &cat[[CC ! z : z in Conjugates(y : Precision:=Precision)] :y in Components(x)];
        end function;
        maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Dimension(A)] >;
        A`HomsToC:=maps;
    end if;
    return A`HomsToC;
end intrinsic;

/* TESTS

    printf "### Testing Homs:";
    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    homs:=HomsToC(A);
    a:=PrimitiveElement(A);
    assert &and[ Abs(Evaluate(f,h(a))) lt 10^-20 : h in homs ];
    printf " all good!\n"; 

*/
