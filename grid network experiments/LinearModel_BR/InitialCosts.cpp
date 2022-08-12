#include "InitialCosts.h"

vector <vector <double>> InitialCosts(int n0, int m0, int Ly0, int H0,
	vector <double> x0, vector<vector <double>> y0) {
	double** C0;
	C0 = new double* [n0]();
	for (size_t i = 0; i < n0; i++)
	{
		C0[i] = new double[n0]();
	}

	for (size_t i = 1; i < Ly0; i++) // Distance between Ly0s
	{
		for (size_t j = 0; j < H0; j++)
		{
			C0[(i - 1) * H0 + j + 1][i * H0 + j + 1] =
				sqrt(pow(x0[(i - 1) * H0 + j + 1] - x0[i * H0 + j + 1], 2.0)
					+ pow(y0[i][j] - y0[i + 1][j], 2.0));

			if (j >= 1)
			{
				C0[(i - 1) * H0 + j + 1][(i - 1) * H0 + j] =
					sqrt(pow(x0[(i - 1) * H0 + j + 1] - x0[(i - 1) * H0 + j], 2.0)
						+ pow(y0[i][j] - y0[i][j - 1], 2.0));

				C0[(i - 1) * H0 + j + 1][i * H0 + j] =
					sqrt(pow(x0[(i - 1) * H0 + j + 1] - x0[i * H0 + j], 2.0)
						+ pow(y0[i][j] - y0[i + 1][j - 1], 2.0));
			}
			if (j <= H0 - 2)
			{
				C0[(i - 1) * H0 + j + 1][(i - 1) * H0 + j + 2] =
					sqrt(pow(x0[(i - 1) * H0 + j + 1] - x0[(i - 1) * H0 + j + 2], 2.0)
						+ pow(y0[i][j] - y0[i][j + 1], 2.0));

				C0[(i - 1) * H0 + j + 1][i * H0 + j + 2] =
					sqrt(pow(x0[(i - 1) * H0 + j + 1] - x0[i * H0 + j + 2], 2.0)
						+ pow(y0[i][j] - y0[i + 1][j + 1], 2.0));
			}
		}
	}

	for (size_t i = 0; i < H0; i++) // Distance from source to the first Ly0
	{
		C0[0][i + 1] =
			sqrt(pow(x0[0] - x0[i + 1], 2.0)
				+ pow(y0[0][0] - y0[1][i], 2.0));
	}

	for (size_t i = 0; i < H0; i++) // distance from last Ly0 to sink
	{
		C0[(Ly0 - 1) * H0 + i + 1][n0 - 1] =
			sqrt(pow(x0[(Ly0 - 1) * H0 + i + 1] - x0[n0 - 1], 2.0)
				+ pow(y0[Ly0][i] - y0[Ly0 + 1][0], 2.0));


		if (i >= 1)
		{
			C0[(Ly0 - 1) * H0 + i + 1][(Ly0 - 1) * H0 + i] =
				sqrt(pow(x0[(Ly0 - 1) * H0 + i + 1] - x0[(Ly0 - 1) * H0 + i], 2.0)
					+ pow(y0[Ly0][i] - y0[Ly0][i - 1], 2.0));
		}
		if (i <= H0 - 2)
		{
			C0[(Ly0 - 1) * H0 + i + 1][(Ly0 - 1) * H0 + i + 2] =
				sqrt(pow(x0[(Ly0 - 1) * H0 + i + 1] - x0[(Ly0 - 1) * H0 + i + 2], 2.0)
					+ pow(y0[Ly0][i] - y0[Ly0][i + 1], 2.0));
		}
	}

	// Defining a Vector for C by0 C0 ----------------------------------------------	
	vector <vector <double>> myC0(n0 + 1);

	for (size_t i = 0; i < n0; i++)
		for (size_t j = 0; j < n0; j++)
		{
			myC0[i].push_back(C0[i][j]);
			if (C0[i][j] > 0)
				myC0[n0].push_back(C0[i][j]);
		}

	return myC0;
}