-- (OGr(3,11), S): OGr(3,11) is embedded in OGr(3,12) as the 3-planes orthogonal to e6-f6, as in
-- [SK, Proposition 5.38]; the generic stabilizer is sl5 = sl(E5), E5 = <e1,...,e5>, killing e6 and f6.
-- The stabilizer of P in sl5 has dimension 10, so the SL5-orbit of [P] has dimension 24 - 10 = 14 = dim X.

-- Clifford action on Delta = wedge^. C^6: e_i = wedge, f_i = contraction
n = 6;
subsetsList = subsets(toList(1..n));
N = #subsetsList;
indexOf = new HashTable from apply(N, i -> (subsetsList#i, i));
wedgeOp = (i,S) -> (
    if member(i,S) then null
    else ( p := #select(S, s -> s < i); sgn := (-1)^p; Snew := sort append(S,i); (sgn, Snew) )
);
contractOp = (i,S) -> (
    if not member(i,S) then null
    else ( pos := position(S, s -> s == i); sgn := (-1)^pos; Snew := delete(i,S); (sgn, Snew) )
);
buildOpMatrix = (opFunc, i) -> (
    M := mutableMatrix(QQ, N, N);
    for j from 0 to N-1 do (
        S := subsetsList#j; r := opFunc(i,S);
        if r =!= null then ( (sgn,Snew) := r; k := indexOf#Snew; M_(k,j) = sgn; );
    );
    matrix M
);
Emat = apply(toList(1..n), i -> buildOpMatrix(wedgeOp,i));
Fmat = apply(toList(1..n), i -> buildOpMatrix(contractOp,i));
cl = v -> sum(n, i -> v_(i,0) * Emat#i) + sum(n, i -> v_(n+i,0) * Fmat#i);
omega = Pm -> product apply(numColumns Pm, k -> cl(Pm_{k}));   -- omega_P acting on spinors

delta = mutableMatrix(QQ, N, 1);
delta_(indexOf#{}, 0) = 1;
delta_(indexOf#{1,2,3,4,5,6}, 0) = 1;
delta = matrix delta;   -- 1 + e123456

J = matrix apply(12, i -> apply(12, j -> if (i<6 and j==i+6) or (j<6 and i==j+6) then 1 else 0));

-- gl6 = gl(E), E = <e1,...,e6>, inside so12: X = diag(A, -A^T)
R = QQ[(symbol a)_(1,1)..(symbol a)_(6,6)];
A = matrix table(6, 6, (i,j) -> a_(i+1,j+1));
X = mutableMatrix(R, 12, 12);
for i from 0 to 5 do for j from 0 to 5 do ( X_(i,j) = A_(i,j); X_(6+j,6+i) = -A_(i,j) );
X = matrix X;
-- equations of the stabilizer of the subspace spanned by the columns of Pm
stabIdeal = Pm -> ( Pp := sub(gens kernel transpose Pm, R); ideal flatten entries (transpose(Pp) * X * sub(Pm, R)) );

-- P = < e1-f1+e6+f6, e2+f5, e3+f4 >, columns in e1..e6, f1..f6
P = transpose matrix {
    {1,0,0,0,0,1, -1,0,0,0,0,1},
    {0,1,0,0,0,0, 0,0,0,0,1,0},
    {0,0,1,0,0,0, 0,0,0,1,0,0}};
w = transpose matrix {{0,0,0,0,0,1, 0,0,0,0,0,-1}};   -- e6-f6
assert(rank P == 3);
assert(transpose(P) * J * P == 0);   -- P isotropic
assert(transpose(w) * J * P == 0);   -- P orthogonal to e6-f6
assert(omega P * delta == 0);        -- omega_P kills delta

-- sl5: a_(i,6) = a_(6,j) = 0 and trace zero
I = trim(stabIdeal P + ideal trace A + ideal apply(6, i -> a_(i+1,6)) + ideal apply(6, j -> a_(6,j+1)));
<< "dim stab_{sl5}(P) = " << dim(R/I) << endl;   -- 10, so the orbit of [P] has dimension 14 = dim X
<< "generic element of the stabilizer (compare with the matrix in the proof):" << endl << submatrix(sub(A, R/I), {0..4}, {0..4}) << endl;
