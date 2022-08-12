#include "parameters.h"
#include "InitialCosts.h"


vector <vector <double>> parameters(int n0, int m0, int Ly, int H0, 
	vector <double> x0, vector<vector <double>> y0, int k0) {

	vector <vector <double>> Output_parameters(n0 + 6);
	//PARAMETERS MasterProblem----------------------------------------------------------
	int countB = 0;
	double** B;
	B = new double* [n0]();
	for (size_t i = 0; i < n0; i++)
		B[i] = new double[m0]();

	vector <vector <double>> CallC0;
	CallC0 = InitialCosts(n0, m0, Ly, H0, x0, y0);


	for (size_t i = 0; i < n0; i++)
		for (size_t j = 0; j < n0; j++)
			if (CallC0[i][j] > 0)
			{
				B[i][countB] = 1;
				B[j][countB] = -1;
				countB++;
			}

	for (size_t i = 0; i < n0; i++)
		for (size_t j = 0; j < m0; j++)
			Output_parameters[i].push_back(B[i][j]);

	//Paarameters Subproblem ------------------------------------------------------------
	//Damage level "d" on regions in n0-th row:
	for (size_t i = 0; i < k0; i++) 
	{
		Output_parameters[n0].push_back(40.00 / pow(1.5, i + 1));
	}
	Output_parameters[n0].push_back(0); //Modification

	//Damage area "q" in (n0+1)-th row:
	//finding max and min of Y :| 
	double Ymax = y0[0][0];
	double Ymin = y0[0][0];
	if (y0[Ly + 1][0] > y0[0][0])
		Ymax = y0[Ly + 1][0];
	else
		Ymin = y0[Ly + 1][0];
	for (size_t i = 1; i < Ly + 1; i++) {
		if (Ymax < y0[i][0])
			Ymax = y0[i][0];
		if (Ymin > y0[i][H0 - 1])
			Ymin = y0[i][H0 - 1];
	}

	double Qarea = (x0[n0 - 1] - x0[0]) + (Ymax - Ymin);
	
	Output_parameters[n0 + 1].push_back(-0.001);
	for (size_t i = 1; i < k0 + 1; i++)
		Output_parameters[n0 + 1].push_back((i)*Qarea / k0);

	//Some useful parameters: 	
	Output_parameters[n0 + 2].push_back(Ymax);
	Output_parameters[n0 + 2].push_back(Ymin);
	Output_parameters[n0 + 2].push_back(Qarea);

	//center of arc on x-coordiante:
	for (size_t i = 0; i < n0; i++)
		for (size_t j = 0; j < n0; j++)
			if (CallC0[i][j] > 0)
				Output_parameters[n0 + 3].push_back((x0[i] + x0[j]) / 2);

	//center of arc on y-coordiante:
	for (size_t i = 0; i < H0; i++)
		Output_parameters[n0 + 4].push_back((y0[0][0] + y0[1][i]) / 2);
	
	for (size_t i = 1; i < Ly; i++)
		for (size_t j = 0; j < H0; j++)
		{
			if (j >= 1)
				Output_parameters[n0 + 4].push_back((y0[i][j] + y0[i][j - 1]) / 2);
			if (j <= H0 - 2)
				Output_parameters[n0 + 4].push_back((y0[i][j] + y0[i][j + 1]) / 2);
			if (j >= 1)
				Output_parameters[n0 + 4].push_back((y0[i][j] + y0[i + 1][j - 1]) / 2);
			Output_parameters[n0 + 4].push_back((y0[i][j] + y0[i + 1][j]) / 2);
			if (j <= H0 - 2)
				Output_parameters[n0 + 4].push_back((y0[i][j] + y0[i + 1][j + 1]) / 2);
		}

	for (size_t i = 0; i < H0; i++)
	{
		if (i >= 1)
			Output_parameters[n0 + 4].push_back((y0[Ly][i] + y0[Ly][i - 1]) / 2);
		if (i <= H0 - 2)
			Output_parameters[n0 + 4].push_back((y0[Ly][i] + y0[Ly][i + 1]) / 2);
		Output_parameters[n0 + 4].push_back((y0[Ly][i] + y0[Ly + 1][0]) / 2);
	}

	for (size_t i = 0; i < m0; i++)
	{
		vector <double> MAX_DIST(4);
		MAX_DIST[0] = abs(Output_parameters[n0 + 3][i] - x0[0]) + abs(Output_parameters[n0 + 4][i] - Ymin);
		MAX_DIST[1] = abs(Output_parameters[n0 + 3][i] - x0[n0 - 1]) + abs(Output_parameters[n0 + 4][i] - Ymin);
		MAX_DIST[2] = abs(Output_parameters[n0 + 3][i] - x0[0]) + abs(Output_parameters[n0 + 4][i] - Ymax);
		MAX_DIST[3] = abs(Output_parameters[n0 + 3][i] - x0[n0 - 1]) + abs(Output_parameters[n0 + 4][i] - Ymax);
		double max00 = MAX_DIST[0];
		for (size_t i = 1; i < MAX_DIST.size(); i++)
			if (max00 < MAX_DIST[i])
				max00 = MAX_DIST[i];
		
		Output_parameters[n0 + 5].push_back(max00);
	}



	return Output_parameters;
}
