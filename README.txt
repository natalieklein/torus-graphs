Matlab code for fitting phase coupling networks using Torus Graphs.

Modifications to work with Octave:
- Install statistics and mapping packages (see octave_setup.m for install commands for these and dependencies)
- In calls to chi2cdf, use 1 - chi2cdf(X, N) as there is no 'upper' argument in Octave statistics package version

A tutorial is included along with this code in Torus_Graph_Tutorial.pdf. We recommend using the tutorial to familiarize yourself with the code.

The code and tutorial supplement our publication:

Klein, N., Orellana, J., Brincat, S., Miller, E. K. and Kass, R. E. (2019). “Torus Graphs for Multivariate Phase Coupling Analysis”, The Annals of Applied Statistics. 

==========
COPYRIGHT
==========

@ 2019 Natalie Klein    neklein@stat.cmu.edu
       Josue Orellana   josue@cmu.edu

N.K. and J.O. contributed equally to this work.
       
This code may not be distributed in any form without prior consent of the authors. 
Use and distribution of the included LFP data is also subject to approval from Brincat and Miller.

While the code has been carefully checked for bugs and applied to many datasets, we do not provide warranties of any kind.  
The authors take no responsibility for support or maintenance of the code, although feedback is appreciated.  
