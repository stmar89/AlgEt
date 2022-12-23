/* vim: set syntax=magma :*/

//freeze;

declare verbose OverOrders,3;

/////////////////////////////////////////////////////
// OverOrders for Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
// and Edgar Costa, MIT
/////////////////////////////////////////////////////

// most of the code here is an implementation of  Hofman, Sircana "On the computation of overorders"

/*TODO
 - prime per prime
 - can I use Check:=false in creating the orders?
*/


import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

declare attributes AlgEtQOrd : MinimalOverOrders,
                              OverOrders;

intrinsic MinimalOverOrders(R::AlgEtQOrd : singular_primes := [], orders := {@ @}) -> SetIndx[AlgEtQOrd]
{Returns the minimal over orders of R given the singular primes of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
if not assigned R`MinimalOverOrders then
    min_oo := { };
    if not IsMaximal(R) then
      zbR := ZBasis(R);
      if singular_primes ne [] then
        pp := [(R!!P) meet (OneIdeal(R)) : P in singular_primes];
        pp := Setseq(Seqset(pp)); //remove duplicates
        pp := [P : P in pp | Index(P, P*P) ne Index(R,P)]; //only the sing ones
        assert2 SequenceToSet(pp) eq SequenceToSet(PrimesAbove(Conductor(R)));
      else
        pp := PrimesAbove(Conductor(R));
      end if;
      for P in pp do
        pot_min_oo := {@ @}; // will contain all potential minimal over-orders.
        pot_min_oo_2 := {@ @}; // will contain all potential minimal over-orders of dimension ge 2
        F, f := ResidueField(P);
        T := MultiplicatorRing(P);
        V,mTV := QuotientVS(T, R, P);
        assert2 forall{ v : v in Basis(V) | mTV((v@@mTV)^2) in V };
        assert2 forall{ v : v in Basis(V) | mTV((v@@mTV)) eq v };
        assert2 forall{ t : t in ZBasis(T) | t-((mTV(t))@@mTV) in R };
        d := Dimension(V);
        //see Proposition 5.3 of Tommy's paper
        if d eq 1 then
          Include(~min_oo, T);
        else
          q:=#F; //need to use q! not the characterstic of F!
          qpow:=hom<V->V | [mTV((v@@mTV)^q) : v in Basis(V)]>;
          eigen_vals:=[e[1] : e in Setseq(Eigenvalues(Matrix(qpow)))];
          eigen_spaces:=[Kernel(hom<V->V | [qpow(v)-e*v : v in Basis(V)]>)
                               : e in eigen_vals]; //in this way there are naturally embedded in V
          subs_1:=[ W: W in &cat[Submodules(E) : E in eigen_spaces] | Dimension(W) eq 1];
          for W in subs_1 do //dim eq 1
            wT:=W.1@@mTV;
            if q eq 2 or mTV(wT^2) in W then
            //doing this square directly on the finite field level would probably be faster.
            //for p eq 2 being a subspace of the eigenspace garantuees that it is mult closed
                S:=Order([wT] cat zbR);
                Include(~pot_min_oo,S);
                Include(~min_oo,S);//necessarly minimal
            end if;
          end for;
          dims := PrimesUpTo(d+1); //the plus one is to prevent issues when d=2.
          subs_2 := Submodules(V : CodimensionLimit := d-2); //we exclude dim 0 and 1
          subs_2 := [W : W in subs_2 | Dimension(W)+1 in dims];
          //only subs of dim a prime number
          //the +1 comes from using (P:P)/R instead of (P:P)/P, see Remark 5.4
          for W in subs_2 do //dim at least 2
            S := Order([(w@@mTV) : w in Basis(W)] cat zbR);
            Include(~pot_min_oo,S);
            Include(~pot_min_oo_2,S);
          end for;
          //we remove non-minimals
          for S in pot_min_oo_2 do
            if not exists {T : T in pot_min_oo | S ne T and T subset S} then
              Include(~min_oo, S);
            end if;
          end for;
        end if;
      end for;
    end if;
    R`MinimalOverOrders := {@ @};
    for S in min_oo do
      i := Index(orders, S);
      if i eq 0 then
        ZBasisLLL(S);
        Include(~R`MinimalOverOrders, S);
      else
        T:=orders[i];
        ZBasisLLL(T);
        Include(~R`MinimalOverOrders, T);
      end if;
    end for;
  end if;
  return R`MinimalOverOrders;
end intrinsic;


intrinsic FindOverOrders_Minimal(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]
{Given an order R returns all the over orders by a recursive search of the minimal overordes. Based on "On the computations of overorders" by TommyHofmann and Carlo Sircana.}
    A := Algebra(R);
    singular_primes := PrimesAbove(MaximalOrder(A)!!Conductor(R));
    queue := {@ R @};
    output:={@ R @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[ MinimalOverOrders(elt : singular_primes := singular_primes, orders := output) : elt in queue ];
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    output:={@ T : T in output | Index(MaximalOrder(A),T) ne 1 @} join {@ MaximalOrder(A) @};
    return output;
  /*OLD
  A := Algebra(R);
  singular_primes := PrimesAbove(MaximalOrder(A)!!Conductor(R));
  queue := {@ R @};
  done := {@  @};
  output := {@ @};
  while #queue gt 0 do
    output join:=  queue;
    done join:= queue;
    for elt in queue do
      output join:= MinimalOverOrders(elt : singular_primes := singular_primes, orders := output);
    end for;
    queue := output diff done;
  end while;
  // remove and add the maximal order to avoid creating it twice.
  output:={@ T : T in output | Index(MaximalOrder(A),T) ne 1 @} join {@ MaximalOrder(A) @};
  return output;
  */
end intrinsic;


intrinsic FindOverOrders(E::AlgEtQOrd:  populateoo_in_oo := false) -> SetIndx[AlgEtQOrd]
{Returns all the overorders of E. The boolean VarArg populateoo_in_oo determines whether to populate the attribute OverOrders of each overorder of E.}
  if not assigned E`OverOrders then
      E`OverOrders := FindOverOrders_Minimal(E);
  end if;

  // there might be a better way to do this
  // like looping over MaximalUnderOrders
  if populateoo_in_oo then
    for i in [1..#E`OverOrders] do
      S := E`OverOrders[i];
      if not assigned S`OverOrders then
        S`OverOrders := {@ T : T in E`OverOrders | S subset T @};
      end if;
    end for;
  end if;
  return E`OverOrders;
end intrinsic;

/*
// for this one, I need Driscriminant and AritmeticRadical

intrinsic pMaximalOrder(O::AlgEtQOrd, p::RngIntElt) -> AlgEtQOrd
{given O, retuns the maximal p over order}
  if (Abs(Integers() ! Discriminant(O)) mod p^2) ne 0 then
    return O;
  end if;

  OO := O;
  // Theorem 6.1.3 Cohen
  while true do
    I := ArithmeticRadical(OO, BaseRing(OO)*p);
    OO := MultiplicatorRing(I);
    if OO eq Order(I) then
      return OO;
    end if;
  end while;
end intrinsic;
*/

/* TESTS

    printf "### OverOrders:";
	  AttachSpec("~/packages_github/AlgEt/spec");

    SetVerbose("OverOrders",1);
    SetAssertions(1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    oo:=FindOverOrders(E);

    printf " all good!\n"; 
 
*/
