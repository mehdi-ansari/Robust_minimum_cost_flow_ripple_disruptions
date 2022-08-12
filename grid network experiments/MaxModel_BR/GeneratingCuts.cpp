#include "GeneratingCuts.h"

vector<double> MYaddingC(int m0, vector<int> k0, int NumberOfAttacks0, vector<double> vals0,
	vector<double> Yoptimal0, vector<double> p1a0, vector<double> p2a0, vector<double> Alp10,
	vector<double> Alp20, vector<vector<double>> d0, vector<vector<double>> q0, vector<double> C0a0) {
	vector<double> addingC0;
	int iter = 0;
	for (size_t i = 0; i < m0; i++)
	{
		if (Yoptimal0[i] > 0)
		{
			addingC0.push_back(vals0[iter]);
			iter++;
		}
		else
		{
			double CC = 0;
			for (size_t t = 0; t < NumberOfAttacks0; t++)
				for (size_t j = 0; j < k0[t]; j++)
					if (fabs(p1a0[i] - Alp10[t]) + fabs(p2a0[i] - Alp20[t]) > q0[t][j] && fabs(p1a0[i] - Alp10[t]) + fabs(p2a0[i] - Alp20[t]) <= q0[t][j + 1]) {
						if (d0[t][j] > CC) {
							CC = d0[t][j];
						}
					}

			addingC0.push_back(C0a0[i] + CC);
		}
	}

	return addingC0;
}