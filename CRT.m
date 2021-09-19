/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ?????, 1;

/*TODO:

*/

// THIS IS A TEMPLATE

//------------
// What do we do here?
//------------

intrinsic ???(???::???) -> ???
{???}
  return ???;
end intrinsic;



//  intrinsic ChineseRemainderTheorem(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
//  {given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J}
//      require IsCoprime(I,J) : "the ideals must be coprime";
//      S:=Order(I);
//      require a in S and b in S:"the elements must lie in order of definition of the ideals";
//      require S eq Order(J): "the ideals must be of the same order";
//      K:=Algebra(S);
//      n:=Degree(K);
//      //I need to modify the ZBasis(S) in a way that One(K) is the first element of Zbasis_S
//      Zbasis_S:=ZBasis(S);
// //      pos:=Position(ZBasis(S),One(K));
// //      if pos ne 1 then
// //          if pos eq 0 then //One(K) not in Zbasis_S
// //              coord:=Coordinates([One(K)],Zbasis_S)[1];
// //              pos:=Position(coord,1);
// //              if pos eq 0 then
// //                  pos:=Position(coord,-1);
// //              end if;
// //  //test
// //  if pos eq 0 then coord; end if;
// //  //assert pos ne 0;
// //  //replacing Zbasis_S[pos] with One(K) since they generate the same Z-Span
// //  end if; 
// //          temp:=Zbasis_S[1];
// //          Zbasis_S[1]:=One(K);
// //          Zbasis_S[pos]:=temp;
// //      end if;
//      M:=Matrix(Zbasis_S);
//      Minv:=M^-1;
//      A:=ChangeRing(Matrix(ZBasis(I))*Minv,Integers());
//      B:=ChangeRing(Matrix(ZBasis(J))*Minv,Integers());
//      I_min:=MinimalInteger(I);
//      J_min:=MinimalInteger(J);
//      g,c1,d1:=XGCD(I_min,J_min);
//      if g ne 1 then
//          C:=VerticalJoin(A,B);
//          H,U:=HermiteForm(C); //U*C = H;
//          z:=ZeroMatrix(Integers(),n,n);
//          s:=ScalarMatrix(n,1);
//          assert2 H eq VerticalJoin(s,z);
//          P:=VerticalJoin(HorizontalJoin(z,s),HorizontalJoin(s,z));
//          U1:=Transpose(U)*P; //I need the (n+1)st column of U1
//          Z:=Transpose(U1)[n+1];
//          X:=Matrix(Integers(),1,n,[Z[i] : i in [1..n]]);
//          Y:=X*A;
//          c:=&+[Y[1,i]*Zbasis_S[i] : i in [1..n]];
//          assert2 c in I;
//          d:=One(K)-c;
//          assert2 d in J;
//      else
//          //g:=c1*I_min+d1*J_min
//          c:=c1*I_min;
//          d:=d1*J_min;
//      end if;
//      e:=a*d+b*c;
//      assert e-a in I;
//      assert e-b in J;
//      return e;
//  end intrinsic;


/* TEST

*/
