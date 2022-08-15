Each folder contains the code files corresponding to grid network experiments for every subproblem formulation in the section of numerical experiments. The coding files in "LinearModel_BR" are recommended to get to know the code pieces by comments; Other directory are associated with different subproblem formulation as described in the manuscript.

The results are generated based on data in "grid network instances" folder. Make sure to put the instances in the proper dierctory. 

The default program is set to generate the instance file name and open them to identify the network structure. But, you can comment out some pieces of code in the main cpp file (C2plusANDcplex.cpp) to randomly generate a new network structure. In this case, you must remove the main "for" loop over the default instances.

The results are automatically recorded in ".csv" format in the root directory with this order:
1.master problem optimal value
2.subproblem optimal value
3.Gap
4.number of iteration
5.number of MIP iteration
6.Total time
7.Total MIP time
x.location of disruption epicenters

For more information, please contact me: meansar@okstate.edu
