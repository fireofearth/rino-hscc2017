# RINO

This is a refactored version of the Robust INner and Outer Approximated Reachability (RINO) software with clean code, testing rig, and added examples. This particular version only has ODE (and not DDE code).

# Additional Notes

This is a library to compute guaranteed inner and outer approximations of reachable sets for uncertain continous-time dynamical systems, with (possibly time-varying) perturbations and control inputs.

It relies on Taylor model based reachability analysis to compute outer envelopes of all possible trajectories of an uncertain system, as implemented in other reachability tools (but with the specificity to rely on affine arithmetic for the evaluation of the Taylor models). Additionally, it uses a generalized mean-value theorem to deduce inner tubes, that contain only states guaranteed to be reached. Finally, it also studies robust versions of these tubes, when there can be control inputs and perturbations to the system.

## Dependencies

You need g++, LAPACK and BLAS installed.

Install the FILIB++ Interval Library, available from http://www2.math.uni-wuppertal.de/~xsc/software/filib.html (we used Version 3.0.2), and set variable $FILIBHOME

Get and unzip the FADBAD++ automatic diffentiation package, available from http://www.fadbad.com/fadbad.html (we used FADBAD++ 2.1), and set variable $FADBADHOME

A modified of the third party package for Affine Arithmetic aaflib-0.1 (http://aaflib.sourceforge.net) has been included in the current directory, because some modifications were needed (additional functions and modifications of trigonometric functions). 
Future plans include separating more cleanly the initial version and our modifications...

## Installing

Go to directory aaflib-0.1 within the current package and compile by "make static". 

Returning to the main directory, you can now compile by "make" and obtain the "main" executable. 
The installation has been mostly tested on MacOS, but should also work on Ubuntu. 

## Running existing examples

Hopefully in the future, a parser will allow to read the models of the systems to analyze from input files. For now, they are defined in ode_def.h/ode_def.cpp, and given some fixed ids.
Running an existing example is then performed at command line, by 

```
./main system_id
```

The corresponding systems (both for ODEs and DDES) are defined in ode_def.h (system and constant parameters) and ode_def.cpp (dimensions of the system, initial conditions, uncertain control inputs and perturbations, whether they are constant or time-varying, and control inputs or perturbations, and finally the integration settings - order of Taylor models, initial final time, time step etc). 
More documentation on how to use these (and better input mechanisms) should come.

## Authors and References

This package, written by [Sylvie Putot](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/), implements the ideas presented in:
-  [HSCC 2019] Inner and Outer Reachability for the Analysis of Control Systems, by Eric Goubault and Sylvie Putot, Proceedings of the 22th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2019, Montreal [ [DOI](https://dl.acm.org/citation.cfm?id=3311794) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc19.pdf) ]
-  [CAV 2018] Inner and Outer Approximating Flowpipes for Delay Differential Equations, by Eric Goubault, Sylvie Putot and Lorenz Sahlmann, Proceedings of 30th International Conference on Computer Aided Verification, CAV 2018, Springer LNCS volume 10982 [ [DOI](https://www.springer.com/us/book/9783319961415) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/cav18.pdf) ]  
-  [HSCC 2017] Forward inner-approximated reachability of non-linear continuous systems, by Eric Goubault and Sylvie Putot, Proceedings of the 20th ACM International Conference on Hybrid Systems: Computation and Control, HSCC 2017 [ [DOI](https://dl.acm.org/citation.cfm?id=3049811) | [pdf](http://www.lix.polytechnique.fr/Labo/Sylvie.Putot/Publications/hscc17.pdf) ]
