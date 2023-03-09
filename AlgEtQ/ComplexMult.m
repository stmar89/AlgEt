/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Complex Multiplication Type for AlgEtQ
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ComplexMult, 1;

/////////////////////////////////////////////////////
// CMType for AlgEtQ
/////////////////////////////////////////////////////

declare attributes AlgEtQ : AllCMTypes;
declare type AlgEtQCMType;
declare attributes AlgEtQCMType : Homs, // homs from Q(Frob) of the isogeny class to CC
                                 CMPositiveElement; // a totally imaginary element b in Q(Frob), which is positive for the cm-type. such b is not unique: b and b' define the same cm-type iff b/b' is totally positive


/////////////////////////////////////////////////////
// Creation, print, access and equality testing for AlgEtQCMType
/////////////////////////////////////////////////////

intrinsic CMType(seq::SeqEnum[Map]) -> AlgEtQCMType
{Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType}
    A:=Domain(seq[1]);
    assert &and[Domain(s) eq A : s in seq[2..#seq]];
    g:=Dimension(A) div 2;
    assert #seq eq g;
    F:=PrimitiveElement(A);
    CC_vals:=[ h(F) : h in seq ];
    assert forall{ c : c in CC_vals | not exists{ d : d in CC_vals | Abs( ComplexConjugate(c)-d) lt 10^(2-Precision(d)) }  };
    PHI:=New(AlgEtQCMType);
    PHI`Homs:=seq;
    return PHI;
end intrinsic;

intrinsic CreateCMType(seq::SeqEnum[Map]) -> AlgEtQCMType
{Given a sequence seq of homomorphisms from a CM-algebra to CC, one per conjugate pair, it returns the corresponding CMType.}
    return CMType(seq);
end intrinsic;

intrinsic CMType( b::AlgEtQElt  ) -> AlgEtQCMType
{Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive, that is, Im(phi(b))>0 for every phi in PHI.}
    assert not IsZeroDivisor(b) and (b eq -ComplexConjugate(b));
    PHI:=New(AlgEtQCMType);
    PHI`CMPositiveElement:=b;
    return PHI;
end intrinsic;

intrinsic CreateCMType( b::AlgEtQElt  ) -> AlgEtQCMType
{Given a totally imginary element b, it returns the CMType PHI for which b is PHI-positive.}
    return CMType(b);
end intrinsic;

intrinsic Print( PHI :: AlgEtQCMType)
{Print the AlgEtQCMType.}
    printf "CMType of the Associative Algebra %o determined by the element %o",Domain(Homs(PHI)[1]),CMPositiveElement(PHI);
end intrinsic;

intrinsic CMPositiveElement( PHI::AlgEtQCMType )->AlgEtQElt
{Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).}
    if not assigned PHI`CMPositiveElement then
        A:=Domain(Homs(PHI)[1]);
        F:=PrimitiveElement(A);
        V:=ComplexConjugate(F);
        e:=F-V;
        if forall{ff: ff in Homs(PHI) | Im(ff(e)) gt 0 } then cmposelt:=e; 
        elif forall{ff: ff in Homs(PHI) | Im(ff(-e)) gt 0} then cmposelt:=-e;
        else 
          R:=EquationOrder(A);  
          bb:=5;
          cnter:=0;
          go:=false;
          repeat
            e1:=Random(R : CoeffRange:=bb);
            e:=e1-ComplexConjugate(e1);
            if not IsZeroDivisor(e) and forall{ff: ff in Homs(PHI) | Im(ff(e)) gt 0} then
                go:=true;
            else
                cnter +:=1;
            end if;
            if cnter eq 100 then
                bb:= 2*bb;
                cnter:=1;
            end if;
          until go;
          cmposelt:=e;
        end if;    
        PHI`CMPositiveElement:=cmposelt;
    end if;
    return PHI`CMPositiveElement;
end intrinsic;

intrinsic CMPosElt( PHI::AlgEtQCMType )->AlgEtQElt
{Given a CMType PHI returns a totally imaginary PHI-positive element (which uniquely determines PHI).}
    return CMPositiveElement(PHI);
end intrinsic;

intrinsic Homs( PHI::AlgEtQCMType : prec:=30 )->SeqEnum[Map]
{Given a AlgEtQCMType PHI returns the sequence of maps to the complex field. The vararg prec (default value 30) determines the precision of the codomains of the maps.}
    if not assigned PHI`Homs then
        b:=CMPositiveElement(PHI);
        A:=Parent(b); 
        homs:=HomsToC(A : Precision:=prec);
        phi:=[ ff : ff in homs | Im(ff( b )) gt 0 ];
        assert #phi eq #homs div 2;
        PHI`Homs:=phi;    
    end if;
    return PHI`Homs;
end intrinsic;

intrinsic 'eq'(PHI1 :: AlgEtQCMType, PHI2::AlgEtQCMType : prec:=30)->BoolElt
{Returns whether two cm types are equal. This happens if and only if the quotient of (any) two CMPositiveElements is totally real and totally positive.}
    A:=Domain(Homs(PHI1)[1]);
    assert forall{ phi : phi in Homs(PHI1) cat Homs(PHI2) | Domain(phi) eq A };
    homs:=HomsToC(A : Precision:=prec);
    b1:=CMPositiveElement(PHI1);
    b2:=CMPositiveElement(PHI2);
    if b1 eq b2 then
        return true;
    else
        b:=b1/b2;
        return forall{ h : h in homs | Re(h(b)) gt 0 };
    end if;
end intrinsic;

intrinsic Precision(PHI :: AlgEtQCMType)->RngIntElt 
{Returns the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).}
    if assigned PHI`Homs then
        phi0:=Homs(PHI);
        return Precision(Codomain(phi0[1]));
    else
        return 30; //the default precision
    end if;
end intrinsic;

intrinsic ChangePrecision(PHI0 :: AlgEtQCMType, prec::RngIntElt )->AlgEtQCMType
{Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).}
    require prec gt 0 : "Precision must be a positive integer";
    PHI:=PHI0;
    if assigned PHI`CMPositiveElement then
        b:=CMPositiveElement(PHI);
        A:=Parent(b); 
        homs:=HomsToC(A : Precision:=prec);
        phi:=[ ff : ff in homs | Im(ff( b )) gt 0 ];
        assert #phi eq #homs div 2;
        PHI`Homs:=phi; //this might over-write the attribute
    elif assigned PHI`Homs then
        phi0:=Homs(PHI);
        prec0:=Precision(Codomain(phi0[1]));
        A:=Domain(phi0[1]);
        FA:=PrimitiveElement(A);
        homs:=HomsToC(A : Precision:=prec);
        phi:=[];
        for ff in phi0 do
            guess:=[ gg : gg in homs | Abs(gg(FA)-ff(FA)) lt 10^-(prec0-1)];
            assert #guess eq 1;
            Append(~phi,guess[1]);
        end for;
        assert #phi eq #homs div 2;
        PHI`Homs:=phi; //this might over-write the attribute
    else
        error "there is a problem with the definition of the CM-Type";
    end if;
    return PHI;
end intrinsic;

intrinsic ChangePrecision(~PHI :: AlgEtQCMType, prec::RngIntElt )
{Changes the precision of the given CM-type, that is, the codomain of each homomorphism will be ComplexField(Precision).}
    require prec gt 0 : "Precision must be a positive integer";
    if assigned PHI`CMPositiveElement then
        b:=CMPositiveElement(PHI);
        A:=Parent(b); 
        homs:=HomsToC(A : Precision:=prec);
        phi:=[ ff : ff in homs | Im(ff( b )) gt 0 ];
        assert #phi eq #homs div 2;
        PHI`Homs:=phi; //this might over-write the attribute
    elif assigned PHI`Homs then
        phi0:=Homs(PHI);
        prec0:=Precision(Codomain(phi0[1]));
        A:=Domain(phi0[1]);
        FA:=PrimitiveElement(A);
        homs:=HomsToC(A : Precision:=prec);
        phi:=[];
        for ff in phi0 do
            guess:=[ gg : gg in homs | Abs(gg(FA)-ff(FA)) lt 10^-(prec0-1)];
            assert #guess eq 1;
            Append(~phi,guess[1]);
        end for;
        assert #phi eq #homs div 2;
        PHI`Homs:=phi; //this might over-write the attribute
    else
        error "there is a problem with the definition of the CM-Type";
    end if;
end intrinsic;

/////////////////////////////////////////////////////
// AllCMTypes 
/////////////////////////////////////////////////////

intrinsic AllCMTypes(A::AlgEtQ : Precision := 30 ) -> SeqEnum[AlgEtQCMType]
{Returns all the AlgEtQCMTypes of A. The vararg Precision determined the precision of the codomain of the maps defining the CMTypes.}
    if not assigned A`AllCMTypes then
        // this uses the fact that if A HasComplexConjugation then HomsToC come in conjugate pairs
        cc:=CartesianProduct(Partition([ h: h in HomsToC(A : Precision:=Precision )],2));
        cc:=[ [ci : ci in c] : c in cc ]; //from tuple to seq
        A`AllCMTypes:=[ CMType(c) : c in cc ];
    end if;
    return A`AllCMTypes;
end intrinsic;

/* TESTS

    printf "### Testing CM-types:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    polys:=[
    x^4+x^2+529,
    x^4+11*x^3+73*x^2+319*x+841,
    x^4-4*x^3+8*x^2-116*x+841,
    x^4+4*x^3+8*x^2+116*x+841,
    x^4-17*x^2+841,
    x^4-x^3+26*x^2-31*x+961,
    x^4-6*x^3+43*x^2-186*x+961,
    x^4-4*x^3+8*x^2-124*x+961,
    x^4+2*x^3+52*x^2+62*x+961,
    x^4+x^3+26*x^2+31*x+961
    ];

    for f in polys do
        A:=EtaleAlgebra(f);
        all:=AllCMTypes(A);
        _:=[ CMPositiveElement(PHI) : PHI in all ];
        for i,j in [1..#all] do
            assert (i eq j) eq (all[i] eq all[j]);
        end for;
        for i,j in [1..#all] do
            assert (i eq j) eq (all[i] eq ChangePrecision(all[j],2*Precision(all[j])));
        end for; 
        printf ".";
    end for;
    printf " all good!\n";

*/
