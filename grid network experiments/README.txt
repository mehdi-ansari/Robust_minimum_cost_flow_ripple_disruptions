Each folder contains the code files corresponding to grid network experiments for every subproblem formulation in the numerical experiments section:
  1. LinearModel_BR: Linear Cost Function with Binary Reformulation subproblem
  2. LinearModel_DF: Linear Cost Function with Decomposed Formulation subproblem
  3. LinearModel_EF: Linear Cost Function with Extended Formulation subproblem
  4. LinearModel_PA: Linear Cost Function with Polynomial Algorithm subproblem
  5. MaxModel_BR: Maximum Cost Function with Binary Reformulation subproblem
  6. MaxModel_EF: Maximum Cost Function with Extended Formulation subproblem

The coding files in "LinearModel_BR" are recommended to get to know the code pieces by comments; Other directory are associated with different subproblem formulation as described in the manuscript.

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
