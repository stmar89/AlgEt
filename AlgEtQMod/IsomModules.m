/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the GNU Lesser General Public License (L-GPL)
//      http://www.gnu.org/licenses/
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3.0 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
// 
// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose IsomModules, 2;

import "../AlgEtQ/Ord.m" : crQZ,crZQ,Columns;
import "PowerBass.m" : is_pure_power_internal;

/* TODO

*/

//------------
// Isomorphism Testing for modules
//------------

intrinsic IsIsomorphic(I::AlgEtQMod,J::AlgEtQMod : Method:="Magma",UseSpecializedMethod:=true) -> BoolElt
{Given two modules I and J returns wheater they are isomorphic.
The vararg UseSpecializedMethod (default true) triggers the use of specialized code in the following two special cases: when I and J are fractional ideals, or when they are modules in V over a Bass order S in K with V = K^s for some s.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
To use the latter option, install julia, add the package Hecke and build the package.
Detailed instructions are provided in the example file in the GitHub repository.
Method should be of the form "julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO, according to the output of the julia command "Hecke.build()").}
//TODO return also the isomorphism. Below there are some comments about it. Add it above as well.

    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(J);
    S:=Order(I);
    require S eq Order(J) : "The modules must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
    
    if UseSpecializedMethod then
        if V eq Algebra(S) and forall{b:b in AbsoluteBasis(V)|m(b) eq b} then 
            // fractional ideals case
            dsr1:=DirectSumRepresentation(I)[1];
            dsr2:=DirectSumRepresentation(J)[1];
            Iid:=dsr1[1];
            Jid:=dsr2[1];
            test,x:=IsIsomorphic(Iid,Jid);
            if not test then
                return false,_;
            else
                assert Iid eq x*Jid;
                m1:=dsr1[2](One(V));
                m2:=dsr2[2](One(V));
                isom:=Hom(V,V,[(b/m1)*(1/x)*m2 : b in AbsoluteBasis(V)]);
                assert Module(S,m,[isom(z):z in ZBasis(I)]) eq J;
                assert Module(S,m,[z@@isom:z in ZBasis(J)]) eq I;
                return true,isom;
            end if;
        end if;
        if IsBass(S) and is_pure_power_internal(m) then
            // power of Bass case
            return IsIsomorphicOverBassOrder(I,J);
        end if;
    end if;

    //require EquationOrder(Algebra(S)) subset S : "Implemented only for modules over orders containing the equation order";
    pi:=PrimitiveElement(Algebra(S));
    if not pi in S then
        pi:=pi*Exponent(Quotient(pi*S + OneIdeal(S),OneIdeal(S)));
    end if;
    pi:=m(pi);
    // pi is now a primitive element of K such that Z[pi] c S.
    
    // based on the following, two Z[pi] modules are isomorphic iff the matrices representing multiplcition by pi are Z-conjugate 
    matI:=crQZ(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I))));
    matJ:=crQZ(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(J)],ZBasis(J))));
    if Method eq "Magma" then
        test,T:=AreGLConjugate(matI,matJ);
        if test then
            // matI*T eq T*matJ
            /*
            // TODO add isom. it should be something like what is below
            zz:=ZBasis(I);
            assert forall{r:r in Rows(T)|SumOfProducts(Eltseq(r),ZBasis(J)) in I};
            mat:=crZQ(matI)^-1*T^-1*matJ;

            imgs:=[ SumOfProducts(Eltseq(r),AbsoluteBasis(V)) : r in Columns(mat) ];
            isom:=Hom(V,V,imgs);
            assert Module(S,m,[isom(z):z in ZBasis(I)]) eq J;
            return true,isom;
            */
            return true,_;
        else
            return false,_;
        end if;
    else //Method eq "julia ...."
        ID:=&cat[ Sprint(Random(0,9)) : i in [1..20] ];
        file:="tmp_" cat ID cat ".txt";
        str:=Sprintf("[\n%o,\n%o\n]\n", Eltseq(matI),Eltseq(matJ));
        fprintf file, "%o",str;
        julia_str:=Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,"");
        assert julia_str[1..7] eq "Vector[";
        indices:=eval(julia_str[8..Index(julia_str,"]")]);
        // file is erase inside the julia script immediately after loading
        if #indices eq 2 then
        // the matrices are not conjugate
            return false,_;
        else
            assert #indices eq 1;
            /*
            //TODO
            // add isom. it should be similar to what we have above
            return true,isom;
            */
            return true,_;
        end if;
    end if;
end intrinsic;

//------------
// Isomorphism Classes for modules
//------------

intrinsic IsomorphismClasses(R::AlgEtQOrd,m::Map : Method:="Magma", UseSpecializedMethod:=true) -> SeqEnum[AlgEtQMod]
{Given an order R in some AlgEtQ K, where K acts on some V, by m:K->V, returns representatives of the isomorphism classes of the S-module lattices in V.
The vararg UseSpecializedMethod (default true) triggers the use of specialized code in the following two special cases: when K=V by using IdealClassMonoid, or when V = K^s for some s and R is Bass, by using IsomorphismClassesOverBassOrder.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
To use the latter option, install julia, add the package Hecke and build the package.
Detailed instructions are provided in the example file in the GitHub repository.
Method should be of the form "julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO, according to the output of the julia command "Hecke.build()").}
    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";

    if UseSpecializedMethod then
        if Codomain(m) eq Domain(m) and forall{b:b in AbsoluteBasis(Domain(m))|m(b) eq b} then 
            // early exit: using ICM
            return [ ModuleFromDirectSum(R,m,[<I,m>]) : I in ICM(R) ];
        end if;
        if IsBass(R) and is_pure_power_internal(m) then
            // early exit: using power-of-Bass code
            return IsomorphismClassesOverBassOrder(R,m);
        end if;
    end if;

    V:=Codomain(m);
    K:=Domain(m);
    pi:=PrimitiveElement(K);
    if not pi in R then
        pi:=pi*Exponent(Quotient(pi*R + OneIdeal(R),OneIdeal(R)));
    end if;
    pi:=m(pi);
    O:=MaximalOrder(K);
    ff:=Conductor(R);
    Vnf:=Components(V);
    Knfpoly:=[ DefiningPolynomial(L) : L in Components(K) ];
    Vnfpoly:=[ DefiningPolynomial(L) : L in Vnf ];
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    mat:=Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(MO)],ZBasis(MO)));
    min_poly:=MinimalPolynomial(mat);
    char_poly:=CharacteristicPolynomial(mat);
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    ind:=[ #Vnfpoly+1 - Index(Reverse(Vnfpoly),fK)  : fK in Knfpoly ]; // last occurence of each number field from the dec of K 
                                                                       // in the decomposition of V
    tup:=<ff_prod[Index(Knfpoly,Vnfpoly[i])] : i in [1..#Vnf]>;
    assert forall{i : i in [1..#tup] | Vnf[i] eq NumberField(Order(tup[i]))};
    Mff:=Module(R,m,tup);
    vprint IsomModules,1 : "candidates:";
    vtime IsomModules,1 : candidates:=IntermediateModulesWithTrivialExtension(MO,Mff,O);
    // candidates contains all representatives of the isomorphism classes of modules whose O-extension is the trivial Steinitz class.
    vprint IsomModules,1 : #candidates;
    classes:=[]; // the final output

    // now we sieve out a minimal set of representatives of the isomorphism classes from candidates.
    if Method eq "Magma" then
        classes_O:=[];
        for M in candidates do
            if not exists{ N : N in classes_O | IsIsomorphic(M,N : Method:=Method) } then
                Append(~classes_O,M);
            end if;
        end for;
    else //Method eq "julia ...."
        ID:=&cat[ Sprint(Random(0,9)) : i in [1..20] ];
        file:="tmp_" cat ID cat ".txt";
        str:="[\n";
        for i->I in candidates do
            mat:=Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I)));
            assert MinimalPolynomial(mat) eq min_poly;
            assert CharacteristicPolynomial(mat) eq char_poly;
            str cat:=Sprintf("%o,\n", Eltseq(mat));
        end for;
        str:=Prune(Prune(str)) cat "\n]\n";
        fprintf file, "%o",str;
        julia_str:=Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,"");
        assert julia_str[1..7] eq "Vector[";
        indices_of_classes_O:=eval(julia_str[8..Index(julia_str,"]")]);
        vprintf IsomModules,1 : "indices_of_classes_k = %o\n",indices_of_classes_O;
        classes_O:=[candidates[i] : i in indices_of_classes_O];
    end if;
    assert2 forall{i:i,j in [1..#classes_O]|(i eq j) eq IsIsomorphic(classes_O[i],classes_O[j])};

    PO,pO:=PicardGroup(O);
    for g in PO do
        if g eq Zero(PO) then
            classes_k:=classes_O; 
        else
            _,Ik:=CoprimeRepresentative(pO(g),ff); // Ik+ff=O, Ik cap ff = Ik*ff
            assert not IsPrincipal(Ik);
            test,Ik_prod:=IsProductOfIdeals(Ik);
            assert test;
            MffIk:=Module(R,m,<i notin ind select ff_prod[Index(Knfpoly,Vnfpoly[i])] 
                                        else ff_prod[Index(Knfpoly,Vnfpoly[i])]*Ik_prod[Index(Knfpoly,Vnfpoly[i])]
                                        : i in [1..#Vnf]>); // MffIk = f1^(s1-1)+f1I1 + ... + fn^(sn-1)+fnIn
            zbMffIk:=ZBasis(MffIk);
            e:=ChineseRemainderTheorem(ff,Ik,One(K),Zero(K)); // e in Ik, e-1 in ff.
            assert not IsZeroDivisor(e);
            e:=Components(e);
            e:=V!< i in ind select e[Index(Knfpoly,Vnfpoly[i])] else Vnf[i] ! 1 : i in [1..#Vnf] >;
            ik:=map<V->V | x:->x*e >;
            // ik induces the isomorphism between (O1^s1 + ... + On^sn)/(f1^s1 + ... + fn^sn),
            // and (O1^(s1-1)+I1 + ... + On^(sn-1)+In)/(f1^(s1-1)+f1I1 + ... + fn^(sn-1)+fnIn),
            // where Ik = I1 + .... + In.
            classes_k:=[ Module(R,m, [ik(z) : z in ZBasis(M)] cat zbMffIk ) : M in classes_O];
            assert2 forall{i:i,j in [1..#classes_k]|(i eq j) eq IsIsomorphic(classes_k[i],classes_k[j])};
        end if;
        Append(~classes,classes_k);
    end for;
    classes:=&cat(classes);
    assert2 forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j])};
    return classes;
end intrinsic;

/* TESTS
    // see all_tests_AlgEtQMod.m

*/
