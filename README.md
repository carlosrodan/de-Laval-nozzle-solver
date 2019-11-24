# This code simulates a convergent-divergent nozzle in Matlab and is based on the book "Computational fluid dynamics. The basics with applications" by Anderson J.D. (Chapter 7). 

--------------------------------------------------------------
*The main scripts to run are NozzleSolverIsent.m for the fully isentropic nozzle and NozzleSolverShock.m for a nozzle with a normal shock in the divergent section.
--------------------------------------------------------------

-- fluxTerms.m calculates the fluxes F1,F2,F3,J2 from U1,U2,U3.

-- mach.m is used for fsolve() to get the Mach number at a certain section A/A*

-- isentFlow.m gives the flow properties for a isentropic section of the nozzle specified by A/A*



*I did this as a project at Linköpings Universitet in the course TMMV56. 

