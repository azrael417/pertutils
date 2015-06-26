# pertutils
This repository contains some optimized implementation of functions from perturbative QCD (such as the QCD beta function) or effective field theories
(e.g. finite volume correction formulas, Luscher Zeta function).
It is compiled in a dynamic library and can used by other programs. I am going to provide a doxygen documentation soon.
An important dependency is my library libmathutils, which I have to keep in a private repository due to licensing issues. I am trying to refactor it in such a way that it does not lose any functionality but does not contain routines from Numerical Recipes any more.
