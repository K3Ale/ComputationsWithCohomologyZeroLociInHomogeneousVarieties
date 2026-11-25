# Macaulay2 files for the PhD thesis of Alessandro Frassineti

This folder contains the *Macaulay2* computations used in the PhD thesis

> A. Frassineti, *On the geometry of zero loci of homogeneous vector bundles: Fano and hyper-Kähler manifolds*, PhD thesis, Università degli Studi di Genova (2026).

## Contents

- `AuxFilePhDThesisFrassineti.m2`: the main ancillary file. It collects all the computations
  mentioned in the thesis, in their order of appearance, labelled by the corresponding
  statement (Chapters 2, 4, 6 and 7).
- `OG*.m2`: standalone scripts for the stabilizer computations of Chapter 7 (Lemma 7.2.14 for
  `OG29pure.m2`, and quasi-homogeneity and cohomogeneity of spinorial Fano manifolds for the
  other files). The main file loads each of them at the point where it is needed:

| File             | Statement in the thesis | Variety                     |
|------------------|-------------------------|-----------------------------|
| `OG38S.m2`       | Proposition 7.2.5       | (OGr(3,8), S+)              |
| `OG29pure.m2`    | Lemma 7.2.14            | singular X_δ ⊆ OGr(2,9)     |
| `OG39S.m2`       | Proposition 7.2.16      | (OGr(3,9), S)               |
| `OG210S.m2`      | Proposition 7.2.19      | OGr(2,10) and (OGr(2,10), S+) |
| `OG310S.m2`      | Proposition 7.2.23      | (OGr(3,10), S+)             |
| `OG310SS.m2`     | Proposition 7.2.24      | (OGr(3,10), S+ ⊕ S+)        |
| `OG310SpSm.m2`   | Proposition 7.2.26      | (OGr(3,10), S+ ⊕ S−)        |
| `OG311S.m2`      | Proposition 7.2.30      | (OGr(3,11), S)              |
| `OG312S.m2`      | Proposition 7.2.35      | (OGr(3,12), S+)             |
| `OG412S.m2`      | Proposition 7.2.38      | (OGr(4,12), S+)             |
| `OG314S.m2`      | Proposition 7.2.44      | (OGr(3,14), S+)             |
| `OG414S.m2`      | Proposition 7.2.46      | (OGr(4,14), S+)             |
| `OG514S.m2`      | Proposition 7.2.47      | (OGr(5,14), S+)             |
| `OG313S.m2`      | Proposition 7.3.2       | [OGr(3,13), S]              |

## How to run

Most computations in the main file use the package
[`CohomologyZeroLociInHomogeneousVarieties`](https://github.com/K3Ale/package-HomogeneousVarieties).
Start *Macaulay2* in this folder and run

```
run "git clone https://github.com/K3Ale/package-HomogeneousVarieties.git"
path = prepend(currentDirectory() | "package-HomogeneousVarieties/", path);
needsPackage "CohomologyZeroLociInHomogeneousVarieties";
```

Then evaluate the computations in `AuxFilePhDThesisFrassineti.m2` one block at a time. Some of
them, for homogeneous varieties of high dimension, take several minutes, and the main file
gives the running time for these. The `OG*.m2` files must stay in the same folder as the main
file, which loads them through `currentDirectory()`.

The `OG*.m2` scripts do not need the package. Each of them can also be run on its own, e.g.

```
M2 --script OG314S.m2
```

and prints the result stated in the thesis. All `OG*.m2` scripts together run in about 30 seconds.

The computations were tested with *Macaulay2* version 1.26.06.

## Conventions

Spinor bundles follow Chapter 1 of the thesis: on OGr(k,2n), `Sp` = S+ = E_{ω_n} and
`Sm` = S− = E_{ω_{n−1}}. `OGr{n,2n}` denotes the component OGr(n,2n)_+ = D_n/P_n.
