/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtQHoms, 1;

declare attributes AlgEtQ : HomsToC;

//------------
// Homomorphism to the Complex Numbers
//------------

intrinsic HomsToC(A::AlgEtQ : Prec:=30)->SeqEnum[Map]
{Returns the sequence of homomorphisms from A to the complex field CC. The precision of CC is given by the optional parameter "Prec". Default value is 30}
    if not assigned A`HomsToC or (Prec ne Precision(Codomain(A`HomsToC[1]))) then
        CC:=ComplexField(Prec);
        images:=function(x)
            return &cat[[CC ! z : z in Conjugates(y : Precision:=Prec)] : y in Components(x)];
        end function;
        maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Dimension(A)] >;
        A`HomsToC:=maps;
    end if;
    return A`HomsToC;
end intrinsic;

/* TESTS

    printf "### Testing Homs:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    homs:=HomsToC(A);
    a:=PrimitiveElement(A);
    assert &and[ Abs(Evaluate(f,h(a))) lt 10^-20 : h in homs ];
    old_prec:=Precision(Codomain(homs[1]));
    new_prec:=10*old_prec;
    homs:=HomsToC(A : Prec:=new_prec);
    assert Precision(Codomain(homs[1])) eq new_prec;
    printf " all good!\n"; 

*/
