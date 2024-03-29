# cubic-fourfolds

This repository houses the code used to work with cubic fourfolds over F2. The vast majority
of the software is written in the Magma computer algebra language.

## Dependancies

- Macaulay2 (For `FunctionalEquationSign`)
- Sage (For certain scripts in `article-scripts`)
- The `magma-parallel-cookbook` is needed to run the parallelized versions of the scripts in
the `surveys/` directory.
- C++ for point_counting functionality. Specifically, `g++` must be available.

## Installation

1. Clone this repository.
2. Use `AttachSpec(".../cubic-fourfolds/CubicLib/CubicLib.spec");`

# Repository contents

### 1. CubicLib

This is where the source code of the library is located.


### 2. ABC

This is a series of questions Jack asked about lines in cubic fourfolds, and the code
set up to answer them.

### 3. article_scripts

This directory contains the scripts that are important for the paper. One must have available
the data generated by the scripts in the *orbit_calculation, surveys* folders. 
A full breakdown is:

- **artin_tate.m :** This script computes the set of observed Artin-Tate special values.


- **feasibility_check.m :** This script compares for various $(n, d, q)$ whether applying
union-find or our methods are computationally feasible. The script is meant to validate
a pair of tables appearing in the article.


- **newton-polygons-survey.m :** Computes the list of Newton Polygons of the zeta functions
observed in the database.


- **weil_compare.m :** This script runs various queries about the set of zeta functions
from cubic fourfolds in relation to the zeta functions of K3-type discussed in 
[Kedlaya-Sutherland].

- **weil_search.sage :** Used to generate the list of Weil polynomials of degrees up to 22
satisfying the conditions outlined in [Kedlaya-Sutherland].

**To generate the data from scratch, these scripts should be run third.**

### 4. discovery

These are scripts for submitting batch jobs on the Discovery cluster at Dartmouth College.
We advise ignoring this directory.

### 5. legacy

This directory contains depreciated code. Please ignore.

### 6. linear\_space\_calculations

This directory contains the two scripts for determining the set of lines (resp. planes)
in P5(F2) passing through the chosen orbit representatives of the cubic fourfolds.

### 7. macaulay2

Contains the Macaulay2 code used to determine the sign of the functional equation of the
zeta functions. The function `FunctionalEquationSign` within `CubicLib` accomplishes this
task for a single example. However, to avoid overhead costs in generating the database,
we use a Macaulay2 script.

### 8. orbit_calculation

Contains scripts used to determine the list of orbit representatives. Two scripts exist for
this task for the sake of timing comparisons, distinguished by the particular choice of 
filtration used for the task.

**To generate the data from scratch, these scripts should be run first.**

### 9. surveys

Contains scripts used to compute various features of each orbit representatives
(Orbit sizes, point counts). A single threaded and parallelize version of the scripts are
available. 

**To generate the data from scratch, these scripts should be run second.**

# Documentation

Specific documentation strings can be found for most intrinsics in CubicLib. In general, the
most relevant functions are as follows:

```
// Reading data into a session.
LoadCubicOrbitData;
ReadOrbitSizeData;
ReadStabilizerSizeData;
ReadZetaFunctions;

// Group theory
CountOrbits;
IsFeasible;
IsFeasibleUnionFind;

// Working with cubics
IsEquivalentCubics;
PointCounts;

// Working with Polynomials/Zeta functions
FunctionalEquationSign;
Charpoly;
IsWeilPolynomial;
TranscendentalFactor;
TranscendentalRank;
IrrationalFactor;
TateTwist;
CubicLPolynomialToPointCounts;
CubicWeilPolynomialToLPolynomial; // (Alias CWTL)
CubicWeilPolynomialToPointCounts;

ArtinTateValue;
IsOrdinary;

// Lines and Planes
LinesThrough;

// Standard forms
ConicFibrationCoefficients;
ConifFibrationForm;
```
