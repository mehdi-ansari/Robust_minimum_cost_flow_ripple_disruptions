#include "GeneratingNetwork.h"

vector<vector <double> > GeneratingNetwork(int n0, int Ly, int H0, int Ux, int Uy)
{
	vector < vector <double> > networkXY(Ly + 3);

	// Generating X-Coordinate values in the range [0,Ux] -------------------------
	double xInt;
	
	for (size_t i = 0; i < n0; i++)
	{
		xInt = rand() % 1000;
		networkXY[0].push_back((xInt / 1000) * Ux);
	}
	sort(networkXY[0].begin(), networkXY[0].end());
	// ----------------------------------------------------------------------------

	// Generating Y-Coordinate values in the range [0,Uy]--------------------------
	double yInt;
	yInt = rand() % 1000;
	networkXY[1].push_back((yInt / 1000) * Uy);
	
	for (size_t i = 2; i < Ly + 2; i++)
	{
		for (size_t j = 0; j < H0; j++)
		{
			yInt = rand() % 1000;
			networkXY[i].push_back((yInt / 1000) * Uy);
		}
	}

	yInt = rand() % 1000;
	networkXY[Ly + 2].push_back((yInt / 1000) * Uy);
	// Ordering Y -----------------------------------------------------------------
	for (size_t l = 2; l < Ly + 2; l++)
	{
		sort(networkXY[l].begin(), networkXY[l].end(), greater <double>());
	}
	// ----------------------------------------------------------------------------

	return networkXY;
}