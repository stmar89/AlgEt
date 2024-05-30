/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, stefano.marseglia89@gmail.com
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

/* TODO

*/

//------------
// Isomorphism Testing for modules
//------------

intrinsic IsIsomorphic(I::AlgEtQMod,J::AlgEtQMod : Method:="Magma") -> BoolElt
{Given two modules I and J returns wheater they are isomorphic.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
To use the latter option, install julia, add the package Hecke and build the package.
Detailed instructions are provided in the example file in the GitHub repository.
Method should be of the form "julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO, according to the output of the julia command "Hecke.build()").}
    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(I);
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
   
    //require EquationOrder(Algebra(S)) subset S : "Implemented only for modules over orders containing the equation order";
    pi:=PrimitiveElement(Algebra(S));
    if not pi in S then
        pi:=pi*Exponent(Quotient(pi*S + OneIdeal(S),OneIdeal(S)));
    end if;
    pi:=m(pi);
    // pi is now a primitive element of K such that Z[pi] c S.
    
    // based on the following, two Z[pi] modules are isomorphic iff the matrices representing multiplcition by pi are Z-conjugte 
    matI:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I))),Integers());
    matJ:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(J)],ZBasis(J))),Integers());
    if Method eq "Magma" then
        test:=AreGLConjugate(matI,matJ);
    else //Method eq "julia ...."
        ID:=&cat[ Sprint(Random(0,9)) : i in [1..20] ];
        file:="tmp_" cat ID cat ".txt";
        str:=Sprintf("[\n%o,\n%o\n]\n", Eltseq(matI),Eltseq(matJ));
        fprintf file, "%o",str;
        indices:=eval(Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,""));
        // file is erase inside the julia script immediately after loading
        if #indices eq 2 then
        // the matrices are not conjugate
            test:=false;
        else
            assert #indices eq 1;
            test:=true;
        end if;
    end if;
    return test;
end intrinsic;

//------------
// Isomorphism Classes for modules
//------------

intrinsic IsomorphismClasses(R::AlgEtQOrd,m::Map : Method:="Magma") -> SeqEnum[AlgEtQMod]
{Given an order R in some AlgEtQ K, where K acts on some V, by m:K->V, returns representatives of hte isomorphism classes of the S-module lattices in V.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
To use the latter option, install julia, add the package Hecke and build the package.
Detailed instructions are provided in the example file in the GitHub repository.
Method should be of the form "julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO, according to the output of the julia command "Hecke.build()").}
    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";
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
    Mff:=Module(R,m,<ff_prod[Index(Knfpoly,Vnfpoly[i])] : i in [1..#Vnf]>);
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
        indices_of_classes_O:=eval(Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,""));
        vprintf IsomModules,1 : "indices_of_classes_k = %o\n",indices_of_classes_O;
        classes_O:=[candidates[i] : i in indices_of_classes_O];
    end if;


    PO,pO:=PicardGroup(O);
    for g in PO do
        if g eq Zero(PO) then
            classes_k:=classes_O; 
        else
            _,Ik:=CoprimeRepresentative(pO(g),ff); // Ik+ff=O, Ik cap ff = Ik*ff
            test,Ik_prod:=IsProductOfIdeals(Ik);
            assert test;
            MffIk:=Module(R,m,<i in ind select 
                                            ff_prod[Index(Knfpoly,Vnfpoly[i])] 
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
        end if;
        Append(~classes,classes_k);
    end for;
    classes:=&cat(classes);
    return classes;
end intrinsic;

/* TESTS
    // see all_tests_AlgEtQMod.m

*/
