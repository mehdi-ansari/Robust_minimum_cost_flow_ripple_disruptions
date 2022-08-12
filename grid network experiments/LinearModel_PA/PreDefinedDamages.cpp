#include "PreDefinedDamages.h"

vector <vector <double>> PreDef_Dmg(vector <double> p1a, vector <double> p2a, vector <double> q, vector <double> d) {
	double x_cand1, x_cand2, x_cand3, x_cand4;

	vector <vector <double>> Loc_cand;

	for (size_t i = 0; i < p1a.size(); i++)
	{
		for (size_t j = 0; j < p1a.size(); j++)
		{
			for (size_t k = 1; k < q.size(); k++)
			{
				for (size_t l = 1; l < q.size(); l++)
				{
					x_cand1 = (p2a[j] + p1a[j] - p2a[i] + p1a[i] + q[k] + q[l]) / 2;
					if ((x_cand1 >= p1a[i] && x_cand1 <= p1a[i] + q[k]) && 
						(x_cand1 >= p1a[j] && x_cand1 <= p1a[j] + q[l]))
					{
						vector <double> coord_cand(2);
						coord_cand[0] = x_cand1;
						coord_cand[1] = x_cand1 + p2a[i] - p1a[i] - q[k];
						Loc_cand.push_back(coord_cand);
					}

					x_cand2 = (p2a[j] + p1a[j] - p2a[i] + p1a[i] + q[k] - q[l]) / 2;
					if ((x_cand2 >= p1a[i] && x_cand2 <= p1a[i] + q[k]) &&
						(x_cand2 >= p1a[j] - q[l] && x_cand2 <= p1a[j]))
					{
						vector <double> coord_cand(2);
						coord_cand[0] = x_cand2;
						coord_cand[1] = x_cand2 + p2a[i] - p1a[i] - q[k];
						Loc_cand.push_back(coord_cand);
					}

					x_cand3 = (p2a[j] + p1a[j] - p2a[i] + p1a[i] - q[k] + q[l]) / 2;
					if ((x_cand3 >= p1a[i] - q[k] && x_cand3 <= p1a[i]) &&
						(x_cand3 >= p1a[j] && x_cand3 <= p1a[j] + q[l]))
					{
						vector <double> coord_cand(2);
						coord_cand[0] = x_cand3;
						coord_cand[1] = x_cand3 + p2a[i] - p1a[i] + q[k];
						Loc_cand.push_back(coord_cand);
					}

					x_cand4 = (p2a[j] + p1a[j] - p2a[i] + p1a[i] - q[k] - q[l]) / 2;
					if ((x_cand4 >= p1a[i] - q[k] && x_cand4 <= p1a[i]) &&
						(x_cand4 >= p1a[j] - q[l] && x_cand4 <= p1a[j]))
					{
						vector <double> coord_cand(2);
						coord_cand[0] = x_cand4;
						coord_cand[1] = x_cand4 + p2a[i] - p1a[i] + q[k];
						Loc_cand.push_back(coord_cand);
					}
				}
			}
		}
	}

	vector <vector <double>> Pre_damage(Loc_cand.size());

	for (size_t t = 0; t < Loc_cand.size(); t++)
	{
		for (size_t i = 0; i < p1a.size(); i++)
		{
			for (size_t j = 0; j < d.size(); j++)
			{
				if (abs(p1a[i] - Loc_cand[t][0]) + abs(p2a[i] - Loc_cand[t][1]) > q[j] &&
					abs(p1a[i] - Loc_cand[t][0]) + abs(p2a[i] - Loc_cand[t][1]) <= q[j + 1])
				{
					Pre_damage[t].push_back(d[j]);
				}
			}
		}
		Pre_damage[t].push_back(Loc_cand[t][0]);
		Pre_damage[t].push_back(Loc_cand[t][1]);
	}

	return Pre_damage;
}