-- Lemma: if X_delta in OGr(2,9) is singular, then delta is pure.
-- Checks: (1) the ten quadrics I of the text define OGr(4,9) (= pure spinors, via the Mukai pairing);
--         (2) J = radical(4-minors of M + (d0,d3,d4,d34)) is the ideal written in the text;
--         (3) I is contained in J.
--

-- Delta = wedge^. <e1..e4>, coordinates d_I
S4 = sort(subsets {1,2,3,4}, x -> (#x, x));
varName = I -> if #I == 0 then "d0" else "d" | concatenate(apply(I, toString));
R = QQ[apply(S4, varName)];
d = I -> R_(varName I);
idx = new HashTable from apply(#S4, k -> (S4#k, k));

-- Clifford action of V9 = E + F + <v9>: e_j wedge, f_j contraction, v9 . e_I = (-1)^|I| e_I
signIns = (j,I) -> if odd number(I, i -> i < j) then -1 else 1;
creat = j -> ( Mx := mutableMatrix(R,16,16); for k to 15 do ( I := S4#k;
    if not member(j,I) then Mx_(idx#(sort append(I,j)), k) = signIns(j,I) ); matrix Mx );
annih = j -> ( Mx := mutableMatrix(R,16,16); for k to 15 do ( I := S4#k;
    if member(j,I) then ( J := select(I, i -> i != j); Mx_(idx#J, k) = signIns(j,J) ) ); matrix Mx );
parop = diagonalMatrix(R, apply(S4, I -> (-1)^(#I)));
gammas = join(apply(1..4, creat), apply(1..4, annih), {parop});

-- Mukai pairing B(e_I, e_J) = sign of e_{rev I} ^ e_J in e_1234
permSign = L -> ( s := 1; L = new MutableList from L;
    for i to #L-1 do for j to #L-2-i do if L#j > L#(j+1) then (t := L#j; L#j = L#(j+1); L#(j+1) = t; s = -s); s );
B = matrix table(16, 16, (a,b) -> ( I := S4#a; J := S4#b;
    if #I + #J == 4 and #(set I * set J) == 0 then permSign(join(reverse I, J)) * 1_R else 0_R ));
delta = transpose matrix {apply(S4, d)};
Ipure = ideal prepend((transpose delta * B * delta)_(0,0), apply(gammas, g -> (transpose delta * B * g * delta)_(0,0)));

-- the ten quadrics of the text
Itext = ideal(
    d{3,4}*d{1,2,4} - d{2,4}*d{1,3,4} + d{1,4}*d{2,3,4} - d{4}*d{1,2,3,4},
    d{3,4}*d{1,2,3} - d{2,3}*d{1,3,4} + d{1,3}*d{2,3,4} - d{3}*d{1,2,3,4},
    d{2,4}*d{1,2,3} - d{2,3}*d{1,2,4} + d{1,2}*d{2,3,4} - d{2}*d{1,2,3,4},
    d{1,4}*d{1,2,3} - d{1,3}*d{1,2,4} + d{1,2}*d{1,3,4} - d{1}*d{1,2,3,4},
    d{4}*d{1,2,3} - d{3}*d{1,2,4} + d{2}*d{1,3,4} - d{1}*d{2,3,4},
    d{1,4}*d{2,3} - d{1,3}*d{2,4} + d{1,2}*d{3,4} - d{}*d{1,2,3,4},
    d{4}*d{2,3} - d{3}*d{2,4} + d{2}*d{3,4} - d{}*d{2,3,4},
    d{4}*d{1,3} - d{3}*d{1,4} + d{1}*d{3,4} - d{}*d{1,3,4},
    d{4}*d{1,2} - d{2}*d{1,4} + d{1}*d{2,4} - d{}*d{1,2,4},
    d{3}*d{1,2} - d{2}*d{1,3} + d{1}*d{2,3} - d{}*d{1,2,3});
<< "(1) I_text = ideal of pure spinors: " << (Itext == Ipure)
   << "; dim = " << dim Itext << ", degree = " << degree Itext << ", prime: " << isPrime Itext << endl;   -- true; 11, 12, true

-- smoothness matrix of X_delta at [P], P = <e1,e2> (with d0 = d3 = d4 = d34 = 0)
M = matrix {
  { d{1,2},     d{1},      d{2},      0,      0,      0,      0,     d{1,3},   d{1,4},   d{2,3},   d{2,4}   },
  { d{1,2,3},  -d{1,3},   -d{2,3},    d{1},   d{2},   0,      0,     0,       -d{1,3,4}, 0,       -d{2,3,4} },
  { d{1,2,4},  -d{1,4},   -d{2,4},    0,      0,     -d{1},  -d{2},  d{1,3,4}, 0,        d{2,3,4}, 0        },
  { d{1,2,3,4}, d{1,3,4},  d{2,3,4},  d{1,4}, d{2,4}, d{1,3}, d{2,3}, 0,       0,        0,        0        }};
J = radical(minors(4, M) + ideal(d{}, d{3}, d{4}, d{3,4}));

-- J as written in the text
Jtext = ideal(d{}, d{3}, d{4}, d{3,4},
    d{2,4}*d{1,3,4} - d{1,4}*d{2,3,4}, d{2,3}*d{1,3,4} - d{1,3}*d{2,3,4}, d{2}*d{1,3,4} - d{1}*d{2,3,4},
    d{2,4}*d{1,2,3} - d{2,3}*d{1,2,4} + d{1,2}*d{2,3,4} - d{2}*d{1,2,3,4},
    d{1,4}*d{1,2,3} - d{1,3}*d{1,2,4} + d{1,2}*d{1,3,4} - d{1}*d{1,2,3,4},
    d{1,4}*d{2,3} - d{1,3}*d{2,4}, d{2}*d{1,4} - d{1}*d{2,4}, d{2}*d{1,3} - d{1}*d{2,3});
<< "(2) J = J_text: " << (J == Jtext) << endl;                     -- true
<< "(3) I subset J: " << isSubset(Itext, J) << endl;               -- true: singular => pure
