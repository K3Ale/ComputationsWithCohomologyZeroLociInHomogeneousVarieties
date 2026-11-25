-- (OGr(3,10), S_+): the stabilizer of P in s = stab_{so10}([delta]) has codimension 13,
-- so the orbit of [P] is dense.

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
    M := mutableMatrix(QQ, N, N);
    for j from 0 to N-1 do (
        S := subsetsList#j; r := opFunc(i,S);
        if r =!= null then ( (sgn,Snew) := r; k := indexOf#Snew; M_(k,j) = sgn; );
    );
    matrix M
);
Emat = apply(toList(1..n), i -> buildOpMatrix(wedgeOp,i));
Fmat = apply(toList(1..n), i -> buildOpMatrix(contractOp,i));
Id32 = id_(QQ^N);
cl = v -> sum(n, i -> v_(i,0) * Emat#i) + sum(n, i -> v_(n+i,0) * Fmat#i);

-- basis of so10 (basis e1..e5, f1..f5) and its action on spinors via u^w -> (1/2)(uw - wu)
pairsA = flatten apply(toList(1..n), i -> apply(toList(1..n), j -> (i,j)));
pairsBC = flatten apply(toList(1..n), i -> apply(toList(i+1..n), j -> (i,j)));
basisList = apply(pairsA, p -> ("A", p)) | apply(pairsBC, p -> ("B", p)) | apply(pairsBC, p -> ("C", p));
vecRepMatrix = bt -> (
    (typ, p) := bt; (i,j) := p;
    M := mutableMatrix(QQ,10,10);
    if typ == "A" then ( M_(i-1, j-1) = 1; M_(5+j-1, 5+i-1) = -1; )
    else if typ == "B" then ( M_(i-1, 5+j-1) = 1; M_(j-1, 5+i-1) = -1; )
    else if typ == "C" then ( M_(5+i-1, j-1) = 1; M_(5+j-1, i-1) = -1; );
    matrix M
);
spinAction = bt -> (
    (typ, p) := bt; (i,j) := p;
    if typ == "A" then ( R := Emat_(i-1)*Fmat_(j-1); if i==j then R = R - (1/2)*Id32; R )
    else if typ == "B" then Emat_(i-1)*Emat_(j-1)
    else Fmat_(i-1)*Fmat_(j-1)
);

deltaPlus = mutableMatrix(QQ, N, 1);
deltaPlus_(indexOf#{}, 0) = 1;
deltaPlus_(indexOf#{1,2,3,4}, 0) = 1;
deltaPlus = matrix deltaPlus;   -- 1 + e1234

-- P = < e1+e2+e5+f3-f4, f1-f2, e4+f5 >
v1 = matrix{{1},{1},{0},{0},{1},{0},{0},{1},{-1},{0}};
v2 = matrix{{0},{0},{0},{0},{0},{1},{-1},{0},{0},{0}};
v3 = matrix{{0},{0},{0},{1},{0},{0},{0},{0},{0},{1}};
Pmat = v1|v2|v3;
J = matrix apply(10, i -> apply(10, j -> if (i<5 and j==i+5) or (j<5 and i==j+5) then 1 else 0));
assert(rank Pmat == 3);
assert(transpose(Pmat) * J * Pmat == 0);               -- P isotropic
assert(cl(v1) * cl(v2) * cl(v3) * deltaPlus == 0);      -- omega_P kills delta

-- s = stab_{so10}([delta]):  X.delta = t delta  (extra unknown t)
Mspin = matrix {apply(basisList, bt -> (spinAction bt) * deltaPlus)} | (-deltaPlus);
dimS = rank source gens kernel Mspin;

-- X(P) in P:  NullGens^T * X * v_k = 0
NullGens = gens kernel transpose Pmat;
Mvec = matrix {apply(basisList, bt -> (
    Xk := vecRepMatrix bt;
    (transpose(NullGens) * (Xk * v1)) || (transpose(NullGens) * (Xk * v2)) || (transpose(NullGens) * (Xk * v3))
    ))} | map(QQ^(3 * numColumns NullGens), QQ^1, 0);
dimStab = rank source gens kernel (Mspin || Mvec);

<< "dim s = " << dimS << ", dim stab_s(P) = " << dimStab << ", codim = " << dimS - dimStab << endl;   -- 30, 17, 13
