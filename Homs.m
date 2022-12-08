/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtHoms, 1;

/*TODO:

*/

declare attributes AlgEt : HomsToC;

//------------
// Homomorphism to the Compelx Numbers
//------------

intrinsic HomsToC(A::AlgEt : Precision:=30)->SeqEnum[Map]
{returns Hom(A,\C) as a sequence of maps. The precision of \C is given by the optional parameter "Precision". Default value is 30}
    if not assigned A`HomsToC then
        require Type(PrimeField(A)) eq "FldRat";
        CC:=ComplexField(Precision);
        images:=function(x)
            return &cat[[CC ! z : z in Conjugates(y : Precision:=Precision)] :y in Components(x)];
        end function;
        maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Dimension(A)] >;
        f:=&*[DefiningPolynomial(L[1]) : L in Components(A)];
        assert &and [ Abs(Evaluate(f,g(PrimitiveElement(A)))) lt 10^-10 : g in maps]; //the precision here is quite arbitrary...
        A`HomsToC:=maps;
    end if;
    return A`HomsToC;
end intrinsic;

/* TEST

    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    homs:=HomsToC(A);
    a:=PrimitiveElement(A);
    [ h(a) : h in homs ];

*/
