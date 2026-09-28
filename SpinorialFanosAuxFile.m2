-- Ausiliary file to support the computations needed in
-- Spinorial Fano manifolds - A. Frassineti, L. Manivel

-- These computations are based on the M2 Package "CohomologyZeroLociInHomogeneousVarieties"
-- written by A. Frassineti

-- This package is available at https://github.com/K3Ale/package-HomogeneousVarieties
-- and it can be downloaded in the following way:
-- Open a M2 environment and launch the following commands

run "git clone https://github.com/K3Ale/package-HomogeneousVarieties.git"
path = prepend("/full/path/to/package-HomogeneousVarieties/", path);
needsPackage "CohomologyZeroLociInHomogeneousVarieties";

-- IMPORTANT:
-- This package relies on Python and particularly on numPy, which should already be available on your machine

-- Now you can run the desired computations, which follow the order of appearance in the paper
-- Note that for high dimensions of the homogeneous variety, these computations can take some minutes

--------------------------- Beginning of the computations ---------------------------------------------------------------
-- Generation of table 2 and Betti numbers
-- index is obtained by index is obtained from the divisibility of chern(1,Y)
-- and by checking that the value of the volume does not allow any further divisibility

-- (OGr(3,8), Sp)
X = OGr{3,8};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y -- in this case, we can conclude that the index is 1 from the P^3-bundle structure on Q^5
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,8), Sp + Sm)
X = OGr{3,8};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,1,0})},{1},X);
Y = embeddedVariety (Sp+Sm);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,9), S)
X = OGr{3,9};
S = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
Y = embeddedVariety S;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(2,10), Sp)
X = OGr{2,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp+Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety (Sp+Sp);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,10), Sp+Sm)
X = OGr{3,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1,0})},{1},X);
Y = embeddedVariety (Sp+Sm);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,10)_+, Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,10)_+, Sp+Sp)
X = OGr{5,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety (Sp+Sp);
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,11), S)
X = OGr{3,11};
S = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Y = embeddedVariety S;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,12), Sp)
X = OGr{3,12};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(4,12), Sp)
X = OGr{4,12};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(6,12)_+, Sp)
X = OGr{6,12};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(3,14), Sp)
X = OGr{3,14};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(4,14), Sp)
X = OGr{4,14};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(5,14), Sp)
X = OGr{5,14};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(7,14)_+, Sp)
X = OGr{7,14};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
Y = embeddedVariety Sp;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y

-- (OGr(7,14)_+, Sm)
X = OGr{7,14};
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1,0})},{1},X);
Y = embeddedVariety Sm;
dim Y, chern(1,Y)
factor volumeFano Y
hodgeNumbers(Y,Verbose => true);
displayHN Y


---------------------------
-- verification of local rigidity

-- (OGr(3,8), Sp)
X = OGr{3,8};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,1,0})},{1},X);
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
S = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1})},{1},X);
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

-- (OGr(3,10), Sp)
X = OGr{3,10};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,1,0})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
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
S = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,1})},{1},X);
F = S;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- in this case, conditions (1) and (2) are not automatically verified, see the proof of  Proposition 24 
-- with that additional argument, we obtain the expected:
cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(3,12), Sp) 
X = OGr{3,12};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
F = Sp;
Y = embeddedVariety F;
K = koszul(X,F);
K/(k -> (cohomology(k*homogeneousTangentBundle X))#0) -- shows the cohomologies of each term of the twisted Koszul complex
-- in this case, conditions (1) and (2) are not automatically verified, see the proof of  Proposition 27 
-- with that additional argument, we obtain the expected:cohomologyRestriction(X,F,homogeneousTangentBundle X)
cohomology homogeneousTangentBundle X
K/(k -> cohomology(k*F))-- shows the cohomologies of each term of the twisted Koszul complex
-- this shows that condition (3) is verified, in fact: 
cohomologyRestriction(X,F,F)
cohomology F

-- (OGr(4,12), Sp)
X = OGr{4,12};
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
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
Sp = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,0,1})},{1},X);
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
Sm = homogeneousVectorBundle({weight(X#"rootSystem",{0,0,0,0,0,1,0})},{1},X);
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



