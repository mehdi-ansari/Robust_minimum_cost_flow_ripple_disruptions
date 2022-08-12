#include <iostream>
#include <string>
#include <cmath>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <vector>
#include <math.h>       /* fmod */
#include <fstream>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <chrono>


#include <ilcplex/ilocplex.h>

using namespace std;

ILOSTLBEGIN

int main() {
	ofstream MyExcelFile;
	MyExcelFile.open("RESULTS_MaxModel_max.ver.csv");
	vector<string> Cities = { "Paris_Edgelist.txt", "Houston_Edgelist.txt", "SanFrancisco_Edgelist.txt", "Chicago_Edgelist.txt", "LosAngeles_Edgelist.txt", "Tokyo_Edgelist.txt", "Madrid_Edgelist.txt", "London_Edgelist.txt",
								"Bogota_Edgelist.txt", "Tehran_Edgelist.txt", "Shangai_Edgelist.txt", "Delhi_Edgelist.txt" };

	for (size_t city = 0; city < Cities.size(); city++)
	{
		//BEGINNING TIME--------------------------------------------
		auto begin = chrono::steady_clock::now();

		int NumberOfAttacks = 2;

		vector <int> k(NumberOfAttacks);	 //# of ripples 
		k[0] = 6;
		k[1] = 4;

		//Damage 40/(1.5)^j:
		vector <vector <double>> d(NumberOfAttacks);
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			for (size_t j = 0; j < k[t]; j++)
			{
				d[t].push_back((40.00 / pow(1.5, (j + 1))));
			}
		}

		//Reading datasets in the directory
		string dataset_name = Cities[city];

		string line;

		vector <vector <double>> Edge_list;

		ifstream Edgelist_file(dataset_name);
		if (Edgelist_file.is_open())
		{
			while (getline(Edgelist_file, line))
			{
				vector <double> XPath;
				istringstream in(line);
				copy(istream_iterator<double>(in), istream_iterator<double>(), back_inserter(XPath));
				if (XPath.size() == 6)
				{
					Edge_list.push_back(XPath);
				}
				XPath.clear();
			}
		}
		else
		{
			cout << dataset_name << " is not open!!" << endl;
		}

		Edgelist_file.close();

		/* Edge_list Guide:
		0: x-coordinate
		1: y-coordinate
		2: tail node
		3: head node
		4: edge number
		5: length
		*/

		vector <double> x_node;
		vector <double> y_node;

		vector <int> node_recorder;

		x_node.push_back(Edge_list[0][0]);		//x-coordinate node location
		y_node.push_back(Edge_list[0][1]);		//y-coordinate node location

		node_recorder.push_back(Edge_list[0][2]);

		double ux = x_node[0];
		double uy = y_node[0];
		double lx = x_node[0];
		double ly = y_node[0];

		double y_ux, x_uy, y_lx, x_ly;

		int m = Edge_list.size();		//Number of edges
		for (size_t i = 0; i < Edge_list.size(); i++)
		{
			// x and y:
			if (i != 0)
			{
				if (Edge_list[i][2] != Edge_list[i - 1][2])
				{
					x_node.push_back(Edge_list[i][0]);
					y_node.push_back(Edge_list[i][1]);

					node_recorder.push_back(Edge_list[i][2]);		//tail node
				}
			}

			//Determining the region of network; upper bounds and lower bounds
			if (ux <= x_node[x_node.size() - 1])
			{
				ux = x_node[x_node.size() - 1];
				y_ux = y_node[y_node.size() - 1];
				//node_ux = node_recorder[node_recorder.size() - 1];
			}
			if (uy <= y_node[y_node.size() - 1])
			{
				uy = y_node[y_node.size() - 1];
				x_uy = x_node[x_node.size() - 1];
				//node_uy = node_recorder[node_recorder.size() - 1];
			}

			//lx and ly:
			if (lx >= x_node[x_node.size() - 1])
			{
				lx = x_node[x_node.size() - 1];
				y_lx = y_node[y_node.size() - 1];
				//node_lx = node_recorder[node_recorder.size() - 1];
			}
			if (ly >= y_node[y_node.size() - 1])
			{
				ly = y_node[y_node.size() - 1];
				x_ly = x_node[x_node.size() - 1];
				//node_ly = node_recorder[node_recorder.size() - 1];
			}
		}

		int n = Edge_list[Edge_list.size() - 1][2];		//number of nodes

		cout << "Number of nodes: " << n << ", Number of arcs: " << m << endl;

		//middle of arcs:
		vector <double> p1a(m);
		vector <double> p2a(m);

		//Cost of arcs
		vector <double> C0a(m);

		//Adjacency_list
		vector <vector <int>> Adj_list(n);
		vector <vector <int>> arc_Adj_list(n);

		for (size_t i = 0; i < Edge_list.size(); i++)
		{
			int keyindex;
			vector<int>::iterator it;
			it = find(node_recorder.begin(), node_recorder.end(), Edge_list[i][3]);
			if (it != node_recorder.end())
				keyindex = distance(node_recorder.begin(), it);
			else
			{
				cout << "Element not found in myvector\n";
				int erty;
				cin >> erty;
			}

			Adj_list[Edge_list[i][2] - 1].push_back(Edge_list[i][3] - 1);
			arc_Adj_list[Edge_list[i][2] - 1].push_back(i);

			p1a[i] = (Edge_list[i][0] + x_node[keyindex]) / 2;
			p2a[i] = (Edge_list[i][1] + y_node[keyindex]) / 2;

			C0a[i] = Edge_list[i][5];
		}

		//_________Breadth-First search________________++++++++++++++++++++++++++++:
		int max_degree = Adj_list[0].size();
		int source = 0;

		//Selecting a source with max degree:
		for (size_t i = 1; i < Adj_list.size(); i++)
		{
			if (max_degree < Adj_list[i].size())
			{
				max_degree = Adj_list[i].size();
				source = i;
			}
		}

		//Algorithm:
		vector <int> BFS_dist(n);
		for (size_t i = 0; i < n; i++)
		{
			BFS_dist[i] = -1;
		}

		vector <int> BFS_temp_list;

		//Initialize
		BFS_dist[source] = 0;
		BFS_temp_list.push_back(source);
		//LOOP:
		while (BFS_temp_list.size() != 0)
		{
			int first_v = BFS_temp_list[0];
			for (size_t i = 0; i < Adj_list[first_v].size(); i++)
			{
				int first_v_neighb = Adj_list[first_v][i];
				if (BFS_dist[first_v_neighb] == -1)
				{
					BFS_temp_list.push_back(first_v_neighb);
					BFS_dist[first_v_neighb] = BFS_dist[first_v] + 1;
				}
			}

			BFS_temp_list.erase(BFS_temp_list.begin());
		}

		int furthest_dist = 0;
		int sink;

		for (size_t i = 0; i < n; i++)
		{
			if (furthest_dist < BFS_dist[i])
			{
				furthest_dist = BFS_dist[i];
				sink = i;
			}
		}
		cout << "Path: " << source << "  " << sink << "  " << furthest_dist << endl;
		//END_________Breadth-First search________________++++++++++++++++++++++++++++

		int indx_source;
		vector<int>::iterator it_source;
		it_source = find(node_recorder.begin(), node_recorder.end(), source + 1);
		if (it_source != node_recorder.end())
			indx_source = distance(node_recorder.begin(), it_source);
		else
		{
			cout << "Element not found in myvector\n";
			int erty;
			cin >> erty;
		}

		int indx_sink;
		vector<int>::iterator it_sink;
		it_sink = find(node_recorder.begin(), node_recorder.end(), sink + 1);
		if (it_sink != node_recorder.end())
			indx_sink = distance(node_recorder.begin(), it_sink);
		else
		{
			cout << "Element not found in myvector\n";
			int erty;
			cin >> erty;
		}

		double Qarea0 = abs(x_node[indx_source] - x_node[indx_sink]) + abs(y_node[indx_source] - y_node[indx_sink]);

		//Largest ripple:
		double Qarea = (ux - lx) + (uy - ly);

		vector <vector <double>> q(NumberOfAttacks);
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			q[t].push_back(-0.001);
			for (size_t j = 1; j < k[t] + 1; j++)
			{
				q[t].push_back((j)*Qarea0 / (k[t]));
			}
		}

		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			q[t][k[t]] = Qarea;
			d[t][k[t] - 1] = 0;
		}

		Edge_list.clear();

		double Total_time = 0;

		//DEfining var and const. for MASTER Problem
		IloEnv env1;
		IloModel model1(env1);

		// VARIABLES ---------------------------------------------------------------
		IloNumVarArray Yvar(env1, m);
		for (size_t i = 0; i < m; i++)
			Yvar[i] = IloNumVar(env1, 0, 1);

		IloNumVar w(env1);
		model1.add(Yvar);		//Never remove it! why?! I do not know!!

		// CONSTRAINTS -------------------------------------------------------------
		IloExprArray constraint(env1, n);
		for (size_t i = 0; i < n; i++)
		{
			constraint[i] = IloExpr(env1);
			constraint[i] += 0;
		}

		for (size_t i = 0; i < arc_Adj_list.size(); i++)
		{
			for (size_t j = 0; j < arc_Adj_list[i].size(); j++)
			{
				constraint[i] += Yvar[arc_Adj_list[i][j]];
				constraint[Adj_list[i][j]] += -Yvar[arc_Adj_list[i][j]];
			}
		}

		model1.add(constraint[source] == 1);
		model1.add(constraint[sink] == -1);
		for (size_t i = 0; i < n; i++)
		{
			if ((i != source) && (i != sink))
			{
				model1.add(constraint[i] == 0);
			}
		}

		constraint.end();

		IloExpr Objective(env1);

		model1.add(IloMinimize(env1, w));
		IloCplex cplex1(model1);

		// Initial Solutions -----------------------------------------------------------
		double Vforcheck = 0;
		for (size_t i = 0; i < m; i++) //better cutting palne - look at the obj of MIP
			Vforcheck += C0a[i] + NumberOfAttacks * d[0][0];
		double Wforcheck = 0;

		// First Cut --------------------------------------------------	
		vector <vector <double>> C;
		C.push_back(C0a);

		vector <double> Alp1(NumberOfAttacks);
		vector <double> Alp2(NumberOfAttacks);

		int MIP_iteration = 0;
		int iteration = 0;
		bool flag = false;

		while (flag == false)
		{
			// MAIN PROBLEM --- MAIN PROBLEM --- MAIN PROBLEM --- MAIN PROBLEM --- 
			//OBJECTIVE FUNCTION --------------------------------------------------------
			for (size_t i = 0; i < m; i++)
			{
				Objective += C[C.size() - 1][i] * Yvar[i];
			}
			model1.add(Objective <= w);
			//SOLVE ----------------------------------------------------------------------
			cplex1.solve();

			Objective.clear();

			Wforcheck = cplex1.getObjValue();
			cout << "The objective value (MAIN) = " << Wforcheck << "  " << cplex1.getStatus() << endl;
			cout << endl;

			//----- Getting Y for next MIP
			//----- Recognizing Indicies of active arcs
			vector<double> Yoptimal(m);
			vector<int> activeArcs;

			for (size_t i = 0; i < m; i++)
			{
				Yoptimal[i] = cplex1.getValue(Yvar[i]);

				if (Yoptimal[i] > 0)
				{
					activeArcs.push_back(i);
					//cout << "y" << i + 1 << "  " << Yoptimal[i] << endl;
				}
			}
			cout << "Active arcs size: " << activeArcs.size() << endl;


			// Heuristics -----------------------------------------------------------------------------
			int MIP_needed = 1;

			if (MIP_needed == 1)
			{
				vector <double> C_added(activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++)
				{
					C_added[i] = C0a[activeArcs[i]];
				}

				vector <double> vals(activeArcs.size());
				// SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP ---
				IloEnv env;
				IloModel model(env);

				//Varibales -------------------------------------------------------------------
				typedef IloArray<IloNumVarArray> FloatMatrix;
				typedef IloArray<FloatMatrix> Float3Matrix;
				Float3Matrix Lambda(env, activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++) {
					Lambda[i] = FloatMatrix(env, NumberOfAttacks);
					for (size_t j = 0; j < NumberOfAttacks; j++)
						Lambda[i][j] = IloNumVarArray(env, k[j], 0, 1, ILOINT);
				}

				IloNumVarArray Ca(env, activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++)
					Ca[i] = IloNumVar(env);

				/*IloNumArray CAsMIPStart(env, activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++)
					CAsMIPStart[i] = C_H[activeArcs[i]];*/

					//constraints ------------------------------------------------------------------	
				for (size_t i = 0; i < activeArcs.size(); i++)
				{
					IloExpr constraint4a(env); // Constraint 4a ---------------------------
					constraint4a += C0a[activeArcs[i]];
					IloExpr constraint4b(env);
					for (size_t t = 0; t < NumberOfAttacks; t++) {
						for (size_t j = 0; j < k[t]; j++) {
							constraint4a += d[t][j] * Lambda[i][t][j];
							constraint4b += Lambda[i][t][j];
						}
					}
					model.add(constraint4a == Ca[i]);
					constraint4a.end();

					model.add(constraint4b == 1);
					constraint4b.end();
				}


				for (size_t t = 0; t < NumberOfAttacks; t++)
				{
					for (size_t i = 0; i < activeArcs.size() - 1; i++)
					{
						for (size_t j = i + 1; j < activeArcs.size(); j++)
						{
							for (size_t r1 = 0; r1 < k[t]; r1++)
							{
								for (size_t r2 = r1; r2 < k[t]; r2++)
								{
									if (abs(p1a[activeArcs[i]] - p1a[activeArcs[j]]) + abs(p2a[activeArcs[i]] - p2a[activeArcs[j]]) > q[t][r1 + 1] + q[t][r2 + 1])
									{
										model.add(Lambda[i][t][r1] + Lambda[j][t][r2] <= 1);
										if (r1 != r2)
										{
											model.add(Lambda[i][t][r2] + Lambda[j][t][r1] <= 1);
										}
									}
								}
							}
						}
					}
				}

				//Objective -------------------------------------------------------------------
				IloExpr obj(env);
				for (size_t i = 0; i < activeArcs.size(); i++)
					obj += Ca[i] * Yoptimal[activeArcs[i]];

				model.add(IloMaximize(env, obj));

				//MIP_begin = clock(); // TIME of beginning to solve MIP
				auto MIP_begin = chrono::steady_clock::now();

				//SOLVE ----------------------------------------------------------------------
				IloCplex cplex(model);
				cplex.setParam(IloCplex::Param::Threads, 4);
				cplex.setParam(IloCplex::TiLim, 600); //TIMElimit of solving
				cplex.solve();
				MIP_iteration++;

				chrono::duration <double> elapsed_secs_MIP = chrono::steady_clock::now() - MIP_begin;
				Total_time += elapsed_secs_MIP.count();

				//COut -----------------------------------------------------------------------
				cout << "The objective value (MIP) = " << cplex.getObjValue() << "  " << cplex.getStatus() << endl;

				Vforcheck = cplex.getObjValue();

				for (size_t i = 0; i < activeArcs.size(); i++)
				{
					vals[i] = cplex.getValue(Ca[i]);
					C_added[i] = vals[i];
				}

				vector <vector <vector <int>>> opt_lambda(activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++)
				{
					vector <vector <int>> t_opt_lambda(NumberOfAttacks);
					for (size_t t = 0; t < NumberOfAttacks; t++)
					{
						for (size_t j = 0; j < k[t]; j++)
						{
							t_opt_lambda[t].push_back(cplex.getValue(Lambda[i][t][j]));
						}
					}
					opt_lambda[i] = t_opt_lambda;
					t_opt_lambda.clear();
				}

				env.end();

				//Identify epicenters ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				for (size_t t = 0; t < NumberOfAttacks; t++)
				{
					double m1 = -Qarea;
					double M1 = Qarea + ux + uy;
					double m2 = -Qarea - uy;
					double M2 = Qarea + ux;

					for (size_t i = 0; i < activeArcs.size(); i++)
					{
						for (size_t j = 0; j < opt_lambda[i].size(); j++)
						{
							if (opt_lambda[i][t][j] == 1)
							{
								// Upper and lower bounds of optimal epicenters l^1-norm-----------------------------------------------------
								if (m1 < -q[t][j + 1] + p1a[activeArcs[i]] + p2a[activeArcs[i]])
									m1 = -q[t][j + 1] + p1a[activeArcs[i]] + p2a[activeArcs[i]];
								if (M1 > q[t][j + 1] + p1a[activeArcs[i]] + p2a[activeArcs[i]])
									M1 = q[t][j + 1] + p1a[activeArcs[i]] + p2a[activeArcs[i]];
								if (m2 < -q[t][j + 1] + p1a[activeArcs[i]] - p2a[activeArcs[i]])
									m2 = -q[t][j + 1] + p1a[activeArcs[i]] - p2a[activeArcs[i]];
								if (M2 > q[t][j + 1] + p1a[activeArcs[i]] - p2a[activeArcs[i]])
									M2 = q[t][j + 1] + p1a[activeArcs[i]] - p2a[activeArcs[i]];
								// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
							}
						}
					}

					Alp1[t] = (m1 + m2) / 2;
					Alp2[t] = (m1 - m2) / 2;

					////Record -----------------------------------------------------------------------

					cout << "Alpha1_" << t << ": " << Alp1[t] << ", Alpha2_" << t << ": " << Alp2[t] << endl;
				}

				C_added.clear();
				//----------------------------------------------------------------------------	

				//adding new C to MAIN problem constraint
				vector <double> addingC;
				int iter = 0;
				for (size_t i = 0; i < m; i++)
				{
					if (Yoptimal[i] > 0)
					{
						addingC.push_back(vals[iter]);
						iter++;
					}
					else
					{
						double CC = 0;
						for (size_t t = 0; t < NumberOfAttacks; t++)
							for (size_t j = 0; j < k[t]; j++)
								if (abs(p1a[i] - Alp1[t]) + abs(p2a[i] - Alp2[t]) > q[t][j] && abs(p1a[i] - Alp1[t]) + abs(p2a[i] - Alp2[t]) <= q[t][j + 1])
									if (d[t][j] > CC) {
										CC = d[t][j];
									}


						addingC.push_back(C0a[i] + CC);
					}
				}


				//-------------------------------------------------------
				C.push_back(addingC);             //adding the last C to LP constraint
				addingC.clear();
				activeArcs.clear();
				//DamagedArcsMatrix.clear();
				//------------------------------------------------------------------------------


				if (Wforcheck + 0.1 >= Vforcheck)  //Stop criteria
					flag = true;
			}

			//time limit:
			chrono::duration <double> TL = chrono::steady_clock::now() - begin;
			cout << "Elaspsed Time: " << TL.count() << "sec." << endl;
			if (TL.count() > 3600)
				flag = true;

			iteration++;
		}
		Objective.end();
		env1.end();

		cout << "--------------------------------------------------" << endl << "Number of iteration = " << iteration << endl;
		cout << "MIP iteration = " << MIP_iteration << endl;
		cout << endl;
		cout << "w = " << Wforcheck << endl;
		cout << "v = " << Vforcheck << endl;
		cout << "gap = " << Vforcheck - Wforcheck << endl;

		chrono::duration <double> elapsed_secs = chrono::steady_clock::now() - begin;
		cout << endl << "Total time: " << elapsed_secs.count() << " sec." << endl;
		cout << "Time to solve MIPs: " << Total_time << " sec." << endl;

		// EXCEL OUTPUT for "SINGLE" attack----------------------------------------------------------------
		MyExcelFile << dataset_name << ",";
		MyExcelFile << n << ",";
		MyExcelFile << m << ",";
		MyExcelFile << furthest_dist << ",";
		MyExcelFile << Wforcheck << ",";
		MyExcelFile << Vforcheck << ",";
		MyExcelFile << Vforcheck - Wforcheck << ",";
		MyExcelFile << iteration << ",";
		MyExcelFile << MIP_iteration << ",";
		MyExcelFile << elapsed_secs.count() << ",";
		MyExcelFile << Total_time << ",";
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			MyExcelFile << d[t].size() << ",";
			MyExcelFile << Alp1[t] << ",";
			MyExcelFile << Alp2[t] << ",";
		}
		MyExcelFile << endl;
	}

	MyExcelFile.close();
	// ---------------- // ----------------- END-----------------------------------
	return 0;
}