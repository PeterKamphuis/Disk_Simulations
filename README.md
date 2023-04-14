# Disk_Simulations
A fortran program to create artificial fits images of a exponentional disk.


This branch contains a Cplus version of the code in an attempt to speed it up.


Current Status:

It is unclear how to compile it

c++ -o emissionmodel_c EmissionModel.C -L/usr/lib/ -lcfitsio

does not throw an error
It runs, not particularly faster than the fortran version. It looks to be even slightly slower. The output is completely wrong however.
