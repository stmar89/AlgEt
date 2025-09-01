/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare attributes AlgEtQOrd:   ICM,
                                ICM_bar;

///# Ideal class monoid
/// Let $R$ be an order in an étale algebra $A$ over $\mathbb{Q}$. Ideal multiplication induces the structure of commutative monoid on the set of ideal classes of $R$, which we then call `ideal class monoid` of $R$.
/// We denote it by $\mathcal{I}(R)$.
/// The unit element of $\mathcal{I}(R)$ is the class of any principal fractional $R$-ideal.
///   
/// In MAGMA there are two ways of computing the ideal class monoid, which are described in the next subsections.

///## Representatives of the ideal class monoid
/// The first method to compute the ideal class monoid returns a sequence of representatives of the ideal classes.
/// We have a partitioning $\mathcal{I}(R) = \bigsqcup_S \mathcal{I}_S(R)$ where the disjoint union is taken over the overorders $S$ of $R$ and $\mathcal{I}_S(R)$ is the subset of $\mathcal{I}(R)$ consisting of ideal classes with multiplicator ring $S$.
/// The computation is then performed by first computing $\mathcal{W}(R)$ and then observing that for each overorder $S$ of $R$, the Picard group $\mathrm{Pic}(S)$ acts freely on $\mathcal{I}_S(R)$ with quotient space $\mathcal{W}_S(R)$.

/// Given an order $S$, returns a sequence of fractional $S$-ideals representing the ideal classes. The parameter `GRH`  determines whether the computation of the Picard groups of the overorders of $S$ is done using the GRH bound.
intrinsic IdealClassMonoid(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Given an order S, returns a sequence of fractional S-ideals representing the ideal classes. The parameter GRH  determines whether the computation of the Picard groups of the overorders of S is done using the GRH bound.}
    return ICM(S : GRH:=GRH);
end intrinsic;

///ditto
intrinsic ICM(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Given an order S, returns a sequence of fractional S-ideals representing the ideal classes. The parameter GRH  determines whether the computation of the Picard groups of the overorders of S is done using the GRH bound.}
    if not assigned S`ICM then
        _:=WKICM(S); //to populate all attributes T`WKICM_bar
        seqOO:=OverOrders(S);
        seqICM:=[];
        for T in seqOO do
            ICM_barT := [(S!!I) : I in ICM_bar(T: GRH:=GRH )];
            seqICM:=seqICM cat ICM_barT;
        end for;
        assert forall{I: I in seqICM | Order(I) eq S};
        for I in seqICM do
            ZBasisLLL(I);
        end for;
        S`ICM:=seqICM;
    end if;
    return S`ICM;
end intrinsic;

/// Given an order $S$, returns a sequence containing all representatives of the ideal classes of fractional $S$-ideals with multiplicator ring $S$.
intrinsic IdealClassesWithPrescribedMultiplicatorRing(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Given an order S, returns a sequence containing all representatives of the ideal classes of fractional S-ideals with multiplicator ring S.}
    return ICM_bar(S : GRH:=GRH);
end intrinsic;

///ditto
intrinsic ICM_bar(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Given an order S, returns a sequence containing all representatives of the ideal classes of fractional S-ideals with multiplicator ring S.}
    if not assigned S`ICM_bar then
        seqWKS_bar:=WKICM_bar(S);
        GS,gS:=PicardGroup(S : GRH:=GRH );
        repS:=[gS(x) : x in GS];
        ICM_barS := &cat[[I*J : I in seqWKS_bar] : J in repS];
        for iI in [1..#ICM_barS] do
            I:=ICM_barS[iI];
            I`MultiplicatorRing:=S;
            ZBasisLLL(I);
        end for;
        assert2 forall{J : J in ICM_barS | MultiplicatorRing(J) eq S};
        assert forall{J : J in ICM_barS | Order(J) eq S};
        S`ICM_bar:=ICM_barS;
    end if;
    return S`ICM_bar;
end intrinsic;

/* TESTS

    printf "### Testing ICM:";
	SetAssertions(2);
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    polys:=[
        x^3-100*x^2-100*x-100,
        x^4+291*x^3-988*x^2-1000*x-1000,
        x^3+31*x^2+43*x+77
        ];
    for f in polys do
        K:=EtaleAlgebra(f);
        E:=EquationOrder(K);
        _:=ICM(E);
        printf ".";
    end for;
    SetAssertions(1);
    printf " all good!";

*/