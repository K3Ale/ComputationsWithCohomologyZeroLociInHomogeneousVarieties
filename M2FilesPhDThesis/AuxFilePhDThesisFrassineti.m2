-- Auxiliary file to support the computations needed in the PhD thesis

-- These computations are based on the M2 Package "CohomologyZeroLociInHomogeneousVarieties"
-- written by A. Frassineti

-- This package is available at https://github.com/K3Ale/package-HomogeneousVarieties
-- and it can be downloaded in the following way:


-- Open an M2 environment in the folder containing this file and launch the following commands

run "git clone https://github.com/K3Ale/package-HomogeneousVarieties.git"
path = prepend(currentDirectory() | "package-HomogeneousVarieties/", path);
needsPackage "CohomologyZeroLociInHomogeneousVarieties";

-- Now you can run the desired computations, which follow the order of appearance in the thesis
-- Note that for homogeneous varieties of high dimension, these computations can take some minutes

-- Convention on spinor bundles (as in Chapter 1): on OGr(k,2n), Sp = S+ = E_{omega_n} and Sm = S- = E_{omega_{n-1}}
-- OGr{n,2n} denotes the component OGr(n,2n)_+ = D_n/P_n

--------------------------- Beginning of the computations ---------------------------------------------------------------

-- elapsedTime shows how long the computation took on my computer
-- I have included that for the computations which may take several minutes

------------------------------------- Chapter 2 -------------------------------------------------------------------------
-- Example 2.2.4
X = Gr{3,10};
F = det dual first bundles X;
Y = embeddedVariety F;
dim Y, chern(1,Y)
elapsedTime displayHN Y -- 263.649s elapsed

------------------------------------- Chapter 4 -------------------------------------------------------------------------
-- Table 4.1 and Proposition 4.3.4
-- be aware that this may take a while to finish 
X = Gr{6,10};
(U,Q) = bundles X;
F = exteriorPower(3,dual U);
Z = embeddedVariety F;
elapsedTime K = koszul(X,F); -- this may take a while
-- list of the partitions needed in Proposition 4.3.4
listOfPartition = {{1, 0, 0, 0},{1, 1, 0, 0},
{2, 0, 0, 0},{2, 1, 0, 0},{2, 1, 1, 0},{2, 2, 0, 0},
{3, 0, 0, 0},{3, 1, 0, 0},{3, 1, 1, 0},{3, 2, 1, 0},
{4, 0, 0, 0},{4, 1, 0, 0},{4, 1, 1, 0},{4, 2, 1, 0},
{4, 2, 2, 0}};
-- conversion to the weight notation
listOfWeights = apply(listOfPartition, p -> {p#0-p#1,p#1-p#2,p#2-p#3,p#3});
listOfBundles = apply(listOfWeights, v -> homogeneousVectorBundle({{0,0,0,0,0} | reverse v},{1},X));

-- Cohomology of End(S)|_Z imposing Serre duality SUMMAND BY SUMMAND:
-- S*dual S = sum m_k E_k and, since omega_Z is trivial (Z hyperkaehler),
-- H^i(Z,E_k) = H^(d-i)(Z,dual E_k)^* for every irreducible E_k.
-- Each irreducible summand gets its own block of unknowns z_(k,j).
endCohomology = (S, d) -> (
    E := S * dual S;
    W := E#"weights";
    m := E#"multiplicities";
    n := #W;
    H := apply(W, w -> take(cohomologyRestriction(X, F, homogeneousVectorBundle({w},{1},X)), d+1));
    -- common ring with one block of unknowns per irreducible summand
    z := symbol z;
    R := ZZ[flatten apply(n, k -> apply(numgens ring (H#k)#0, j -> z_(k,j)))];
    off := 0;
    H = for k from 0 to n-1 list (
        A := ring (H#k)#0;
        phi := map(R, A, apply(numgens A, j -> R_(off + j)));
        off = off + numgens A;
        apply(H#k, e -> phi e)
        );
    -- position of the dual summand
    dualIndex := k -> position(W, w -> w == ((dual homogeneousVectorBundle({W#k},{1},X))#"weights")#0);
    I := ideal flatten apply(n, k -> apply(d+1, i -> H#k#i - H#(dualIndex k)#(d-i)));
    apply(d+1, i -> (sum(n, k -> m#k * H#k#i)) % I)
);

for i from 0 to #listOfBundles-1 do (
    S = listOfBundles#i;
    << listOfPartition#i << endl;
    << endCohomology(S, 4) << endl << endl;
);


-------------------------------- Chapter 6 -------------------------------------------------------------------------

-- Proposition 6.3.2
X = OGr{2,9};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
cohomologyRestriction(X, S, U') -- 8

-- Lemma 6.4.1
X = homogeneousVariety(rootSystemG2,set {1});
O1 = homogeneousVectorBundle({{0,1}},{1},X);
--
F = O1*O1;
g2 = embeddedVariety(F);
dim g2, chern(1,g2)
displayHN g2
eulerCharacteristicTangent (1,g2) -- chi(T_Y) = -62
volumeFano g2 -- 36
cohomologyRestriction(X,F,O1) -- h^0(-K) = 14


-- Lemma 6.4.2
X = SGr{2,8};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,1,0}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0,0}},{1},X);
--
F = S + U' + O1; 
sb5 = embeddedVariety(F);
dim sb5, chern(1,sb5)
displayHN sb5 
eulerCharacteristicTangent (1,sb5) -- chi(T_Y) = -28
volumeFano sb5 -- 57
cohomologyRestriction(X,F,O1) -- h^0(-K) = 18

-- Lemma 6.4.3
X = SGr{3,8};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,0,1,0}},{1},X);
--
F = U' + exteriorPower(2,U') + S;
sc1 = embeddedVariety(F);
dim sc1, chern(1,sc1)
displayHN sc1 
eulerCharacteristicTangent (1,sc1) -- chi(T_Y) = -25
volumeFano sc1 -- 66
cohomologyRestriction(X,F,O1) -- h^0(-K) = 20

-- Lemma 6.4.4
X = OGr{2,7};
O1 = homogeneousVectorBundle({{0,1,0}},{1},X);
F = O1 + O1 + O1;
b8 = embeddedVariety(F);
dim b8, chern(1,b8)
displayHN b8 
eulerCharacteristicTangent (1,b8) -- chi(T_Y) = -33
volumeFano b8 -- 56 
cohomologyRestriction(X,F,O1) -- h^0(-K) = 18

-- Lemma 6.4.5
X = OGr{3,7};
O1 = homogeneousVectorBundle({{0,0,1}},{1},X);
F = O1*O1 + O1*O1*O1;
oy2 = embeddedVariety(F);
dim oy2, chern(1,oy2)
displayHN oy2 
eulerCharacteristicTangent (1,oy2) -- chi(T_Y) = -108
volumeFano oy2 -- 12
cohomologyRestriction(X,F,O1) -- h^0(-K) = 8

-- Lemma 6.4.6
X = homogeneousVariety(rootSystemG2,set{2});
O1 = homogeneousVectorBundle({{1,0}},{1},X);
--
F = O1*O1*O1*O1;
g1 = embeddedVariety F;
dim g1, chern(1,g1)
displayHN g1
eulerCharacteristicTangent (1,g1) -- chi(T_Y) = -160
volumeFano g1 -- 8
cohomologyRestriction(X,F,O1) -- h^0(-K) = 7

-- Lemma 6.4.7
X = OGr{6,12};
U' = homogeneousVectorBundle({{1,0,0,0,0,0}},{1},X);
O1 = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
F = O1+O1+O1+O1+O1 + (dual U')*O1;
oz5 = embeddedVariety(F);
dim oz5, chern(1,oz5)
displayHN oz5 
eulerCharacteristicTangent (1,oz5) -- chi(T_Y) = -40
volumeFano oz5 -- 42
cohomologyRestriction(X,F,O1) -- h^0(-K) = 15

-----------------------------------------------------------------------------

-- Proposition 6.5.1
X = OGr{3,9};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,0,1,0}},{1},X);
cotangent = homogeneousCotangentBundle X;
F = S+S+S+S;
F1 = embeddedVariety(F);
dim F1, chern(1,F1)
displayHN F1 
eulerCharacteristicTangent (1,F1) -- chi(T_Y) = -12 
volumeFano F1 -- 116
cohomologyRestriction(X,F,O1) -- h^0(-K) = 30

-- condition (a)
cohomologyRestriction(X,F, dual cotangent) -- H^1(Y,T_X|Y) = 0
-- condition (b')
M = matrix(F#"koszul"/(k -> first cohomology (k*F))) -- M_{p,q} = H^p(\wedge^q E^\vee \otimes E)
-- each element in the diagonal except for the first is zero

-----------------------------------------------------------------------------

-- Proposition 6.5.2
X = OGr{5,10};
U' = homogeneousVectorBundle({{1,0,0,0,0}},{1},X);
O1 = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
cotangent = homogeneousCotangentBundle X;
F = O1*O1 + O1+O1+O1+O1+O1;
F2 = embeddedVariety(F);
dim F2, chern(1,F2)
displayHN F2 
eulerCharacteristicTangent (1,F2) -- chi(T_Y) = -65
volumeFano F2 -- 24
cohomologyRestriction(X,F,O1) -- h^0(-K) = 11

-- condition (a)
cohomologyRestriction(X,F, dual cotangent) -- H^1(Y,T_X|Y) = 0
-- condition (b')
M = matrix(F#"koszul"/(k -> first cohomology (k*F))) -- M_{p,q} = H^p(\wedge^q E^\vee \otimes E)
-- each element in the diagonal except for the first is zero

-- Proof of h^0(T_F2) = 0 
-- condition (v1)
(cohomologyRestriction(X,F,exteriorPower(3,cotangent)*O1))#0
(cohomologyRestriction(X,F,exteriorPower(2,cotangent)))#1
(cohomologyRestriction(X,F,exteriorPower(1,cotangent)*dual O1))#2
(cohomologyRestriction(X,F,exteriorPower(0,cotangent)*dual(O1*O1)))#3
-- condition (v2)
(cohomologyRestriction(X,F,exteriorPower(2,cotangent)*dual(O1)))#1
(cohomologyRestriction(X,F,exteriorPower(1,cotangent)*dual(O1*O1)))#2
(cohomologyRestriction(X,F,exteriorPower(0,cotangent)*dual(O1*O1*O1)))#3
-- condition (v3)
(cohomologyRestriction(X,F,exteriorPower(1,cotangent)*dual(O1*O1*O1)))#2
(cohomologyRestriction(X,F,exteriorPower(0,cotangent)*dual(O1*O1*O1*O1)))#3
-- condition (v4)
(cohomologyRestriction(X,F,dual(O1*O1*O1*O1*O1)))#3

--------------------------------------------------------------

-- Proposition 6.5.3
X = OGr{2,7};
U' = homogeneousVectorBundle({{1,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0}},{1},X);
cotangent = homogeneousCotangentBundle X;
F = symmetricPower(2,S);
F3 = embeddedVariety(F);
dim F3, chern(1,F3)
displayHN F3
eulerCharacteristicTangent (1,F3) -- chi(T_Y) = -14
volumeFano F3 -- 72
cohomologyRestriction(X,F,O1) -- h^0(-K) = 21

-- condition (a) (Failure of local completeness)
cohomologyRestriction(X,F, dual cotangent) -- H^1(Y,T_X|Y) = 1
cohomologyRestriction(X,F,F) -- H^1(Y,F|Y) = 0


-- Model in OGr(2,9)
X = OGr{2,9};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0,0}},{1},X);
cotangent = homogeneousCotangentBundle X;
F = symmetricPower(2,U')+S;
F3 = embeddedVariety(F);
dim F3, chern(1,F3)
displayHN F3
eulerCharacteristicTangent (1,F3) -- chi(T_Y) = -14
volumeFano F3 -- 72
cohomologyRestriction(X,F,O1) -- h^0(-K) = 21

-- condition (a) (Failure of local completeness)
cohomologyRestriction(X,F, dual cotangent) -- H^1(Y,T_X|Y) = 1
cohomologyRestriction(X,F,F) -- H^1(Y,F|Y) = 0


-- Family in SGr(2,8)
X = SGr{2,8};
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
E = homogeneousVectorBundle({{0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0,0}},{1},X);
cotangent = homogeneousCotangentBundle X;
F = U'+E;
sb4 = embeddedVariety(F);
dim sb4, chern(1,sb4)
displayHN sb4
eulerCharacteristicTangent (1,sb4) -- chi(T_Y) = -14
volumeFano sb4 -- 72
cohomologyRestriction(X,F,O1) -- h^0(-K) = 21


-------------------------------- Picard rank greater than one
-- Section 6.6

X = OGr{2,11};
U' = homogeneousVectorBundle({{1,0,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0,0,0}},{1},X);
F = S + symmetricPower(2,U');
ob3 = embeddedVariety(F);
dim ob3, chern(1,ob3)
displayHN ob3 
eulerCharacteristicTangent (1,ob3) -- chi(T_Y) = -26
volumeFano ob3 -- 80
cohomologyRestriction(X,F,O1) -- h^0(-K) = 23

X = OGr{2,11};
U' = homogeneousVectorBundle({{1,0,0,0,0}},{1},X);
S = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0,0,0}},{1},X);
F = S + O1 + O1 + O1;
ob4 = embeddedVariety(F);
dim ob4, chern(1,ob4)
displayHN ob4 
eulerCharacteristicTangent (1,ob4) -- chi(T_Y) = -36
volumeFano ob4 -- 70
cohomologyRestriction(X,F,O1) -- h^0(-K) = 21

-- (oy1): model [P^1 x P^1 x P^1 x P^1 x P^1, O(1,1,1,1,1)], see Section 6.6
X = Gr{1,2} ** Gr{1,2} ** Gr{1,2} ** Gr{1,2} ** Gr{1,2};
O1 = homogeneousVectorBundle({{1,1,1,1,1}},{1},X);
F = O1;
oy1 = embeddedVariety(F);
dim oy1, chern(1,oy1)
displayHN oy1 
eulerCharacteristicTangent (1,oy1) -- chi(T_Y) = -16
volumeFano oy1 -- 120
cohomologyRestriction(X,F,O1) -- h^0(-K) = 31

-- (ow11): model [Gr(2,7), (Sym^2 U^vee)^2], written as [OGr(2,7), Sym^2 U^vee], see Section 6.6
X = OGr{2,7};
U' = homogeneousVectorBundle({{1,0,0}},{1},X);
O1 = homogeneousVectorBundle({{0,1,0}},{1},X);
F = symmetricPower(2,U');
ow11 = embeddedVariety(F);
dim ow11, chern(1,ow11)
displayHN ow11 
eulerCharacteristicTangent (1,ow11) -- chi(T_Y) = -4
volumeFano ow11 -- 80
cohomologyRestriction(X,F,O1) -- h^0(-K) = 21

-- (ob7): model [Fl(1,3,4), O(2,2)], by Proposition 7.2.22 and triality, see Section 6.6
X = Fl{1,3,4};
O11 = homogeneousVectorBundle({{1,0,1}},{1},X);
F = O11 * O11;
ob7 = embeddedVariety(F);
dim ob7, chern(1,ob7)
displayHN ob7 
eulerCharacteristicTangent (1,ob7) -- chi(T_Y) = -68
volumeFano ob7 -- 40
cohomologyRestriction(X,F,O11) -- h^0(-K) = 15

X = OGr{3,8};
O1 = homogeneousVectorBundle({{0,0,1,1}},{1},X);
L = homogeneousVectorBundle({{0,0,1,0}},{1},X);
L' = homogeneousVectorBundle({{0,0,0,1}},{1},X);
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
F = O1 + L + L + L' + L';
oe6 = embeddedVariety(F);
dim oe6, chern(1,oe6)
displayHN oe6 
eulerCharacteristicTangent (1,oe6) -- chi(T_Y) = -23
volumeFano oe6 -- 100
cohomologyRestriction(X,F,O1) -- h^0(-K) = 27

X = OGr{3,8};
O1 = homogeneousVectorBundle({{0,0,1,1}},{1},X);
L = homogeneousVectorBundle({{0,0,1,0}},{1},X);
L' = homogeneousVectorBundle({{0,0,0,1}},{1},X);
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
F = L*L + L + L' + L' + L';
oe9 = embeddedVariety(F);
dim oe9, chern(1,oe9)
displayHN oe9
eulerCharacteristicTangent (1,oe9) -- chi(T_Y) = -20
volumeFano oe9 -- 100
cohomologyRestriction(X,F,O1) -- h^0(-K) = 27

X = OGr{3,8};
O1 = homogeneousVectorBundle({{0,0,1,1}},{1},X);
L = homogeneousVectorBundle({{0,0,1,0}},{1},X);
L' = homogeneousVectorBundle({{0,0,0,1}},{1},X);
U' = homogeneousVectorBundle({{1,0,0,0}},{1},X);
F = exteriorPower(2,U') + L + L';
oe11 = embeddedVariety(F);
dim oe11, chern(1,oe11)
displayHN oe11 
eulerCharacteristicTangent (1,oe11) -- chi(T_Y) = -12
volumeFano oe11 -- 132
cohomologyRestriction(X,F,O1) -- h^0(-K) = 33



----------------------------------- Chapter 7 -------------------------------------------------------------------------

-- Generation of table 7.2 and Betti numbers
-- index is obtained from the divisibility of chern(1,Y)
-- and by checking that the value of the volume does not allow any further divisibility

-- (OGr(3,8), Sp)
X = OGr{3,8};
Sp = homogeneousVectorBundle({{0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y -- in this case, we can conclude that the index is 1 from the P^3-bundle structure on Q^5
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,8), Sp + Sm)
X = OGr{3,8};
Sp = homogeneousVectorBundle({{0,0,0,1}},{1},X);
Sm = homogeneousVectorBundle({{0,0,1,0}},{1},X);
Y = embeddedVariety (Sp+Sm);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,9), S)
X = OGr{3,9};
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
Y = embeddedVariety S;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(2,10), Sp)
X = OGr{2,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp+Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety (Sp+Sp);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp+Sm)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Sm = homogeneousVectorBundle({{0,0,0,1,0}},{1},X);
Y = embeddedVariety (Sp+Sm);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,10)_+, Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,10)_+, Sp+Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety (Sp+Sp);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,11), S)
X = OGr{3,11};
S = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Y = embeddedVariety S;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,12), Sp)
X = OGr{3,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(4,12), Sp)
X = OGr{4,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(6,12)_+, Sp)
X = OGr{6,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,14), Sp)
X = OGr{3,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(4,14), Sp)
X = OGr{4,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,14), Sp)
X = OGr{5,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(7,14)_+, Sp)
X = OGr{7,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(7,14)_+, Sm)
X = OGr{7,14};
Sm = homogeneousVectorBundle({{0,0,0,0,0,1,0}},{1},X);
Y = embeddedVariety Sm;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y


---------------------------
-- verification of local rigidity
-- Theorem 7.1.5 (see also Propositions 7.2.32 and 7.2.37)

-- (OGr(3,8), Sp)
X = OGr{3,8};
Sp = homogeneousVectorBundle({{0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,8), Sp + Sm)
X = OGr{3,8};
Sp = homogeneousVectorBundle({{0,0,0,1}},{1},X);
Sm = homogeneousVectorBundle({{0,0,1,0}},{1},X);
F = Sp+Sm;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3') is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,9), S)
X = OGr{3,9};
S = homogeneousVectorBundle({{0,0,0,1}},{1},X);
F = S;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(2,10), Sp)
X = OGr{2,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,10), Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,10), Sp+Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = Sp+Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3') is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,10), Sp+Sm)
X = OGr{3,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
Sm = homogeneousVectorBundle({{0,0,0,1,0}},{1},X);
F = Sp+Sm;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3') is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(5,10)_+, Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(5,10)_+, Sp+Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = Sp+Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3') is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,11), S) 
X = OGr{3,11};
S = homogeneousVectorBundle({{0,0,0,0,1}},{1},X);
F = S;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- in this case, conditions (1) and (2) are not automatically verified, see the proof of  Proposition 7.2.32 
-- with that additional argument, we obtain the expected:
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,12), Sp) 
X = OGr{3,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- in this case, conditions (1) and (2) are not automatically verified, see the proof of  Proposition 7.2.37 
-- with that additional argument, we obtain the expected:
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(4,12), Sp)
X = OGr{4,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(6,12)_+, Sp)
X = OGr{6,12};
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,14), Sp)
X = OGr{3,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(4,14), Sp)
X = OGr{4,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(5,14), Sp)
X = OGr{5,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(7,14)_+, Sp)
X = OGr{7,14};
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(7,14)_+, Sm)
X = OGr{7,14};
Sm = homogeneousVectorBundle({{0,0,0,0,0,1,0}},{1},X);
F = Sm;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that conditions (1) and (2) are verified, in fact: 
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F


-- additional checks and computations needed

-- Proposition 7.2.1
-- X = (OGr(2,7),S) is connected
S = homogeneousVectorBundle({{0,0,1}},{1},OGr{2,7});
cohomologyRestriction(OGr{2,7},S,structureSheaf OGr{2,7}) -- cohomologies of O_X

-- Proposition 7.2.5
load (currentDirectory() | "OG38S.m2")

-- Remark 7.2.10
Z = OGr{3,8};
O1 = homogeneousVectorBundle({{0,0,1,1}},{1},Z);
L = homogeneousVectorBundle({{0,0,1,0}},{1},Z);
L' = homogeneousVectorBundle({{0,0,0,1}},{1},Z);
U' = homogeneousVectorBundle({{1,0,0,0}},{1},Z);
F = L + L + L + L' + L' + L';
X = embeddedVariety(F);
dim X, chern(1,X)
displayHN X -- Picard rank 2
volumeFano X -- 28
cohomologyRestriction(Z,F,O1) -- h^0(-K) = 17
-- Since X is a Fano threefold with Picard rank 2, 
-- volume 28 and h^0(-K) = 17, it must be MM 2-21


-- Lemma 7.2.14
load (currentDirectory() | "OG29pure.m2")

-- Proposition 7.2.16
load (currentDirectory() | "OG39S.m2")

-- Proposition 7.2.19
load (currentDirectory() | "OG210S.m2")

-- Proposition 7.2.21
-- X = (OGr(2,10), S+ + S+) is connected
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},OGr{2,10});
cohomologyRestriction(OGr{2,10},Sp+Sp,structureSheaf OGr{2,10}) -- cohomologies of O_X

-- Proposition 7.2.22
-- X = (OGr(2,10), S+ + S-) is connected
Sp = homogeneousVectorBundle({{0,0,0,0,1}},{1},OGr{2,10});
Sm = homogeneousVectorBundle({{0,0,0,1,0}},{1},OGr{2,10});
cohomologyRestriction(OGr{2,10},Sp+Sm,structureSheaf OGr{2,10}) -- cohomologies of O_X

-- Proposition 7.2.23
load (currentDirectory() | "OG310S.m2")

-- Proposition 7.2.24
load (currentDirectory() | "OG310SS.m2")

-- Proposition 7.2.26
load (currentDirectory() | "OG310SpSm.m2")

-- Proposition 7.2.27
-- X = (OGr(5,10)_-, S+) is connected
-- by the outer automorphism, (OGr(5,10)_-, S+) = (OGr(5,10)_+, S-)
Sm = homogeneousVectorBundle({{0,0,0,1,0}},{1},OGr{5,10});
cohomologyRestriction(OGr{5,10},Sm,structureSheaf OGr{5,10}) -- cohomologies of O_X

-- Proposition 7.2.28
-- X = (OGr(2,11), S) is connected
S = homogeneousVectorBundle({{0,0,0,0,1}},{1},OGr{2,11});
cohomologyRestriction(OGr{2,11},S,structureSheaf OGr{2,11}) -- cohomologies of O_X

-- Proposition 7.2.30
load (currentDirectory() | "OG311S.m2")

-- Proposition 7.2.33
-- X = (OGr(2,12), S+) is connected
Sp = homogeneousVectorBundle({{0,0,0,0,0,1}},{1},OGr{2,12});
cohomologyRestriction(OGr{2,12},Sp,structureSheaf OGr{2,12}) -- cohomologies of O_X

-- Proposition 7.2.35
load (currentDirectory() | "OG312S.m2")

-- Proposition 7.2.38
load (currentDirectory() | "OG412S.m2")

-- Proposition 7.2.41
-- X = (OGr(2,14), S+) has two connected components
Sp = homogeneousVectorBundle({{0,0,0,0,0,0,1}},{1},OGr{2,14});
cohomologyRestriction(OGr{2,14},Sp,structureSheaf OGr{2,14}) -- cohomologies of O_X

-- Proposition 7.2.44
load (currentDirectory() | "OG314S.m2")

-- Proposition 7.2.46
load (currentDirectory() | "OG414S.m2")

-- Proposition 7.2.47 
load (currentDirectory() | "OG514S.m2")

-- Proposition 7.3.2
X = embeddedVariety homogeneousVectorBundle({{0,0,0,0,0,1}},{1},OGr{3,13});
eulerCharacteristicTangent(1,X) -- chi(T_X) = 15
load (currentDirectory() | "OG313S.m2")

-- Proposition 7.3.3
X = embeddedVariety homogeneousVectorBundle({{0,0,0,0,0,1}},{1},OGr{4,13});
eulerCharacteristicTangent(1,X) -- chi(T_X) = 15
