-- (OGr(4,14), S_+): the stabilizer of Q in g2 x g2 is a copy of sl_2,
-- so the G2 x G2 action has cohomogeneity at most one.

K = QQ;

-- g2 x g2 = stabilizer of delta in so14 (basis e1..e7, f1..f7), 28 parameters
R = K[symbol x11,symbol x12,symbol x13,symbol x21,symbol x22,symbol x23,symbol x31,symbol x32,symbol x33,
       symbol y11,symbol y12,symbol y13,symbol y21,symbol y22,symbol y23,symbol y31,symbol y32,symbol y33,
       symbol a,symbol b,symbol c,symbol d,symbol e,symbol f,symbol u,symbol v,symbol w,symbol l,symbol m,symbol n];
A1 = matrix{
  {x11, x12, x13,   0,   0,   0, a},
  {x21, x22, x23,   0,   0,   0, b},
  {x31, x32, x33,   0,   0,   0, c},
  {  0,   0,   0, y11, y12, y13, d},
  {  0,   0,   0, y21, y22, y23, e},
  {  0,   0,   0, y31, y32, y33, f},
  {  u,   v,   w,   l,   m,   n, 0}
};
A2 = matrix{
  {  0,  w, -v,   0,   0,   0, -a},
  { -w,  0,  u,   0,   0,   0, -b},
  {  v, -u,  0,   0,   0,   0, -c},
  {  0,  0,  0,   0,  -n,   m,  d},
  {  0,  0,  0,   n,   0,  -l,  e},
  {  0,  0,  0,  -m,   l,   0,  f},
  {  a,  b,  c,  -d,  -e,  -f,  0}
};
A3 = matrix{
  {  0, -c,  b,   0,   0,   0,  u},
  {  c,  0, -a,   0,   0,   0,  v},
  { -b,  a,  0,   0,   0,   0,  w},
  {  0,  0,  0,   0,   f,  -e, -l},
  {  0,  0,  0,  -f,   0,   d, -m},
  {  0,  0,  0,   e,  -d,   0, -n},
  { -u, -v, -w,   l,   m,   n,  0}
};
A = (A1 | A2) || (A3 | -transpose(A1));

-- Clifford action on Delta = wedge^. C^7: e_k = wedge, f_k = contraction
evenSubs = sort select(subsets(toList(1..7)), s -> even(#s));
oddSubs  = sort select(subsets(toList(1..7)), s -> odd(#s));
evenIdx = hashTable apply(#evenSubs, i -> evenSubs#i => i);
oddIdx  = hashTable apply(#oddSubs,  i -> oddSubs#i => i);
cliffMat = (k, src, tgtIdx, isWedge) -> (
    M := mutableMatrix(K, #src, #src);
    for j from 0 to #src-1 do (
        S := src#j;
        if isWedge and not member(k,S) then
            M_(tgtIdx#(sort append(S,k)), j) = (-1)^(#select(S, s -> s < k));
        if not isWedge and member(k,S) then
            M_(tgtIdx#(delete(k,S)), j) = (-1)^(position(S, s -> s == k));
        );
    matrix M
    );
GammaPlus  = join(apply(7, i -> cliffMat(i+1, evenSubs, oddIdx, true)),
                  apply(7, i -> cliffMat(i+1, evenSubs, oddIdx, false)));   -- Delta_+ -> Delta_-
GammaMinus = join(apply(7, i -> cliffMat(i+1, oddSubs, evenIdx, true)),
                  apply(7, i -> cliffMat(i+1, oddSubs, evenIdx, false)));   -- Delta_- -> Delta_+
actPlus  = v -> sum apply(14, i -> (v#i) * (GammaPlus#i));
actMinus = v -> sum apply(14, i -> (v#i) * (GammaMinus#i));

delta = mutableMatrix(K, #evenSubs, 1);
scan({{},{1,2,3,7},{4,5,6,7},{1,2,3,4,5,6}}, S -> delta_(evenIdx#S, 0) = 1);
delta = matrix delta;   -- 1 + e1237 + e4567 + e123456

-- omega_Q . delta = v1.(v2.(v3.(v4.delta)))
checkDeltaProduct = Q0 -> (
    vs := apply(4, i -> flatten entries Q0_{i});
    (actMinus vs#0) * (actPlus vs#1) * (actMinus vs#2) * (actPlus vs#3) * delta
    );

JK = (map(K^7,K^7,0) | id_(K^7)) || (id_(K^7) | map(K^7,K^7,0));

-- Q = < p1, p2, p3, e7 >, columns in e1..e7, f1..f7
--   p1 = e1+e4+f1-f4, p2 = e2+e5+f2-f5, p3 = e3+e6+f3-f6
Q = matrix{
  {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {0,0,0,1},
  {1,0,0,0},
  {0,1,0,0},
  {0,0,1,0},
  {-1,0,0,0},
  {0,-1,0,0},
  {0,0,-1,0},
  {0,0,0,0}
};
assert(rank Q == 4);
assert(transpose(Q) * JK * Q == 0);       -- Q isotropic
assert(checkDeltaProduct Q == 0);         -- Q is a point of Y

-- X(Q) in Q  <=>  Qperp^T * A * Q = 0  (trace conditions x33, y33 of the A1 blocks)
Qperp = sub(gens kernel transpose(Q), R);
conds = transpose(Qperp) * A * sub(Q, R);
I = trim((ideal flatten entries conds) + ideal(x11+x22+x33, y11+y22+y33));
<< "dim stab_{g2 x g2}(Q) = " << dim(R/I) << endl;   -- 3, so the orbit of Q has codimension one
