## List of instrinsics in AlgEtQMtrx/Mtrx.m:

> <pre><b>Print</b>(x::AlgEtQMtrx)</pre>
<em>Print the element.</em>

> <pre><b>Universe</b>(x::AlgEtQMtrx) -> AlgEtQ</pre>
<em>Returns the étale algebra to which the entries of the matrix belongs to.</em>

> <pre><b>NumberOfRows</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
<em>Returns the number of rows of the matrix.</em>

> <pre><b>NumberOfColumns</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
<em>Returns the number of colums of the matrix.</em>

> <pre><b>Components</b>(x::AlgEtQMtrx) -> Tup</pre>
<em>Returns the tuple of components of x. If x is defined over K=K1x...xKn with each Ki a number field, the ith compoenent is a matrix over Ki.</em>

> <pre><b>Entries</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtElt]]</pre>
<em>Returns a sequence of sequences containing the entries of the matrix.</em>

> <pre><b>Matrix</b>(entries::SeqEnum[SeqEnum[AlgEtQElt]]) -> AlgEtQMtrx</pre>
<em>Given a sequence of s sequences of r elements of an étale algebra, it returns the corresponding sxr matrix.</em>

> <pre><b>Matrix</b>(s::RngIntElt, r::RngIntElt, Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx</pre>
<em>Given a sequence of sxr elements of an étale algebra, it returns the corresponding sxr matrix.</em>

> <pre><b>DiagonalMatrix</b>(Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx</pre>
<em>Given a sequence of elements of an étale algebra, it returns the corresponding diagonal matrix.</em>

> <pre><b>'[,]'</b>(x::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQElt</pre>
<em>Returns x[i,j].</em>

> <pre><b>Rows</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]</pre>
<em>Returns the sequence of rows.</em>

> <pre><b>Columns</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]</pre>
<em>Returns the sequence of columns.</em>

> <pre><b>Transpose</b>(x::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>Returns the transpose matrix.</em>

> <pre><b>IsSquareMatrix</b>(x::AlgEtQMtrx) -> BoolElt</pre>
<em>Returns whether the matrix is a square matrix.</em>

> <pre><b>ScalarMatrix</b>(n::RngIntElt,s::AlgEtQElt) -> AlgEtQMtrx</pre>
<em>Returns the scalar matrix s\*I where I is the nxn identity matrix.</em>

> <pre><b>IdentityMatrix</b>(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx</pre>
<em>Returns the nxn identity matrix over the étale algebra A.</em>

> <pre><b>ZeroMatrix</b>(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx</pre>
<em>Returns the nxn identity matrix over the étale algebra A.</em>

> <pre><b>RandomMatrix</b>(A::AlgEtQ, s::RngIntElt, r::RngIntElt : bd:=3) -> AlgEtQMtrx</pre>
<em>Returns a random sxr identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.</em>

> <pre><b>RandomMatrix</b>(A::AlgEtQ, n::RngIntElt : bd:=3) -> AlgEtQMtrx</pre>
<em>Returns a random square nxn identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.</em>

> <pre><b>IsInvertible</b>(x::AlgEtQMtrx) -> BoolElt,AlgEtQMtrx</pre>
<em>Returns whether the matrix is invertible and, if so, it returns also the inverse.</em>

> <pre><b>Rank</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
<em>The rank of x, which is defined as the minimum of the rank of his components.</em>

> <pre><b>Trace</b>(x::AlgEtQMtrx) -> AlgEtQElt</pre>
<em>Determinant.</em>

> <pre><b>Determinant</b>(x::AlgEtQMtrx) -> AlgEtQElt</pre>
<em>Determinant.</em>

> <pre><b>'eq'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> BoolElt</pre>
<em>Is x1=x2 ?</em>

> <pre><b>'+'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1+x2.</em>

> <pre><b>'-'</b>(x1::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>-x1.</em>

> <pre><b>'-'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1-x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'+'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1+x2.</em>

> <pre><b>'-'</b>(x1::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>-x1.</em>

> <pre><b>'-'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1-x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::.,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQMtrx,x2::.) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::RngIntElt,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::FldRatElt,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQMtrx,x2::RngIntElt) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'*'</b>(x1::AlgEtQMtrx,x2::FldRatElt) -> AlgEtQMtrx</pre>
<em>x1\*x2.</em>

> <pre><b>'^'</b>(x::AlgEtQMtrx,n::RngIntElt) -> AlgEtQMtrx</pre>
<em>x^n.</em>

> <pre><b>'&+'</b>(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx</pre>
<em>Given a sequence of AlgEtQMtrx returns the sum of the matrices.</em>

> <pre><b>'&*'</b>(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx</pre>
<em>Given a sequence of AlgEtQMtrx returns the product of the matrices.</em>

> <pre><b>InsertBlock</b>(x::AlgEtQMtrx,y::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQMtrx</pre>
<em>Return the matrix obtained by inserting matrix Y in X at position [i,j].</em>

> <pre><b>Solution</b>(x::AlgEtQMtrx,y::AlgEtQMtrx)->AlgEtQMtrx</pre>
<em>Returns a solution V to the system V.x=y.</em>

> <pre><b>???</b>(x::AlgEtQMtrx) -> ???</pre>
<em>????.</em>


