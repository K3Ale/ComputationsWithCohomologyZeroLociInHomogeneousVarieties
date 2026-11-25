-- (OGr(3,10), S_+ + S_-): the stabilizer of P in sl4 (the copy of sl4 in so10 described in the
-- proof of Proposition 7.2.26) has dimension 4, so the orbit of [P] has dimension 15 - 4 = 11 = dim X
-- and is dense. The computation is done over QQ(z), with z = zeta as in delta_- = e5 + zeta e12345.

F = frac(QQ[symbol z]);

-- Clifford action on Delta = wedge^. C^5: e_i = wedge, f_i = contraction
n = 5;
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
    M := mutableMatrix(F, N, N);
    for j from 0 to N-1 do (
        S := subsetsList#j; r := opFunc(i,S);
        if r =!= null then ( (sgn,Snew) := r; k := indexOf#Snew; M_(k,j) = sgn; );
    );
    matrix M
);
Emat = apply(toList(1..n), i -> buildOpMatrix(wedgeOp,i));
Fmat = apply(toList(1..n), i -> buildOpMatrix(contractOp,i));
cl = v -> sum(n, i -> v_(i,0) * Emat#i) + sum(n, i -> v_(n+i,0) * Fmat#i);

deltaPlus = mutableMatrix(F, N, 1);
deltaPlus_(indexOf#{}, 0) = 1;
deltaPlus_(indexOf#{1,2,3,4}, 0) = 1;
deltaPlus = matrix deltaPlus;     -- 1 + e1234
deltaMinus = mutableMatrix(F, N, 1);
deltaMinus_(indexOf#{5}, 0) = 1;
deltaMinus_(indexOf#{1,2,3,4,5}, 0) = z;
deltaMinus = matrix deltaMinus;   -- e5 + zeta e12345

-- P = < zeta e1 + zeta f4 + (1-zeta) f5, e2 - f3, zeta e4 + zeta e5 - f1 >, columns in e1..e5, f1..f5
v1 = transpose matrix{{z,0,0,0,0, 0,0,0,z,1-z}};
v2 = transpose matrix{{0,1,0,0,0, 0,0,-1,0,0}};
v3 = transpose matrix{{0,0,0,z,z, -1,0,0,0,0}};
Pmat = v1|v2|v3;
J = matrix apply(10, i -> apply(10, j -> if (i<5 and j==i+5) or (j<5 and i==j+5) then 1_F else 0_F));
assert(rank Pmat == 3);
assert(transpose(Pmat) * J * Pmat == 0);                 -- P isotropic
assert(cl(v1) * cl(v2) * cl(v3) * deltaPlus == 0);        -- omega_P kills delta_+
assert(cl(v1) * cl(v2) * cl(v3) * deltaMinus == 0);       -- omega_P kills delta_-

-- sl4 = sl(E4) inside so10, E4 = <e1,...,e4>: X = diag(A, 0, -A^T, 0), trace A = 0
R = F[(symbol a)_(1,1)..(symbol a)_(4,4)];
A = matrix table(4, 4, (i,j) -> a_(i+1,j+1));
X = mutableMatrix(R, 10, 10);
for i from 0 to 3 do for j from 0 to 3 do ( X_(i,j) = A_(i,j); X_(5+j,5+i) = -A_(i,j) );
X = matrix X;
assert(transpose(X) * sub(J,R) + sub(J,R) * X == 0);   -- X in so10

-- X(P) in P  <=>  Pperp^T * X * P = 0
Pperp = sub(gens kernel transpose Pmat, R);
I = trim(ideal flatten entries (transpose(Pperp) * X * sub(Pmat, R)) + ideal trace A);
<< "dim stab_{sl4}(P) = " << dim(R/I) << endl;   -- 4, so the orbit of [P] has dimension 11 = dim X

-- the matrices written in the proof (basis e1, e2, e3, e4):
--   | 0  0  -zeta d  0 |
--   | b -a   c       d |
--   | 0  0   a       0 |
--   | 0  0   b       0 |
Itex = ideal(a_(1,1), a_(1,2), a_(1,3) + z*a_(2,4), a_(1,4),
             a_(2,1) - a_(4,3), a_(2,2) + a_(3,3),
             a_(3,1), a_(3,2), a_(3,4),
             a_(4,1), a_(4,2), a_(4,4));
<< "the stabilizer is the space of matrices written in the proof: " << (I == trim Itex) << endl;   -- true
