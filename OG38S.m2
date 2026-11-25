-- (OGr(3,8), S_+): the stabilizer of P in spin7 has dimension 13, so the orbit of [P]
-- (dimension 21 - 13 = 8) is dense in X.
-- Basis of V8: (e1, e2, e3, e4, f1, f2, f3, f4) as in [SK, (5.46)].

-- spin7 inside so8, as in [SK, (5.46)] (entry (4, f3) corrected from -b_{2,3} to -b_{3,4})
R = QQ[(symbol a)_{1,1}..(symbol a)_{4,3},(symbol a)_{1,4},(symbol a)_{2,4},(symbol a)_{3,4},(symbol b)_{1,2},(symbol b)_{1,3},(symbol b)_{1,4},(symbol b)_{2,3},(symbol b)_{2,4},(symbol b)_{3,4}];
A = matrix {
      { a_{1,1}, a_{1,2},  a_{1,3},  a_{1,4},                   0,        b_{1,2},  b_{1,3},  b_{1,4} },
      { a_{2,1}, a_{2,2},  a_{2,3},  a_{2,4},                  -b_{1,2},  0,        b_{2,3},  b_{2,4} },
      { a_{3,1}, a_{3,2},  a_{3,3},  a_{3,4},                  -b_{1,3}, -b_{2,3},  0,        b_{3,4} },
      { a_{4,1}, a_{4,2},  a_{4,3}, -a_{1,1}-a_{2,2}-a_{3,3},  -b_{1,4}, -b_{2,4}, -b_{3,4},  0       },
      { 0,       b_{3,4}, -b_{2,4},  b_{2,3},                  -a_{1,1}, -a_{2,1}, -a_{3,1}, -a_{4,1} },
      {-b_{3,4}, 0,        b_{1,4}, -b_{1,3},                  -a_{1,2}, -a_{2,2}, -a_{3,2}, -a_{4,2} },
      { b_{2,4},-b_{1,4},  0,        b_{1,2},                  -a_{1,3}, -a_{2,3}, -a_{3,3}, -a_{4,3} },
      {-b_{2,3}, b_{1,3}, -b_{1,2},  0,                        -a_{1,4}, -a_{2,4}, -a_{3,4}, a_{1,1}+a_{2,2}+a_{3,3} }};

-- P = < e1+f2, e2-f1, e3-f4 >
P = transpose matrix {
    {1,0,0,0, 0,1,0,0},
    {0,1,0,0,-1,0,0,0},
    {0,0,1,0, 0,0,0,-1}};

-- X(P) in P  <=>  Pperp^T * A * P = 0
Pperp = sub(gens kernel transpose P, R);
I = trim ideal flatten entries (transpose(Pperp) * A * sub(P, R));
<< "conditions on spin7 (compare with the 8 conditions in the proof):" << endl << gens gb I << endl;
<< "dim stab_{spin7}(P) = " << dim(R/I) << endl;   -- 13, so the orbit of [P] has dimension 8 = dim X
