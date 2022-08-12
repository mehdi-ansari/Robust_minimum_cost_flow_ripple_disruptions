#include "PrintNetwork.h"

int PrintNetwork(vector <double> x0, vector <vector <double> > y0, int Ly, int H0) {
	cout << "Source" << endl;
	cout << "(" << x0[0] << " , " << y0[0][0] << ") ";
	cout << endl;
	cout << endl;
	for (size_t i = 1; i < Ly + 1; i++)
	{
		cout << "Layer " << i << endl;
		for (size_t j = 0; j < H0; j++)
		{
			cout << "(" << x0[(i - 1) * H0 + j + 1] << " , " << y0[i][j] << ") ";
		}
		cout << endl;
		cout << endl;
	}

	cout << "Sink" << endl;
	cout << "(" << x0[x0.size() - 1] << " , " << y0[Ly + 1][0] << ") " << endl;

	return 0;
}