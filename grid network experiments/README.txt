Each folder contains the code files of every subproblem formulation in the numerical experiments of paper. "LinearModel_BR" recommended to get to know properly the pieces of codes by comments; The main difference between experiments is about the solving approach of subproblem.

Make sure to put the instances of "grid network instances" folder in the root dierctory. 

The default program is set to generate the instance file name and open them to identify the network structure. But, you can uncomment some pieces of code in the main cpp file (C2plusANDcplex.cpp) to randomly generate a new network structure. In this case, you must remove the main "for" loop over the default instances.

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