// Mehdi - Create 1/8/2020
/*
*Decomposition Formulations (DF) of subproblem
*Multiple disruptions
*L1-norm
*Active arcs used
*/

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <vector>
#include <math.h>       /* fmod */
#include <fstream>
#include <algorithm>
#include<chrono>

#include "GeneratingNetwork.h"
#include "PrintNetwork.h"
#include "InitialCosts.h"
#include "parameters.h"

#include "GeneratingCuts.h"	

#include <sstream>

namespace patch
{
	template < typename T > std::string to_string(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

using namespace std;


ILOSTLBEGIN

//Generating dataset names to open
vector <string> Inputnames(){
	vector <string> dataset_names(80);
	int dataset_counter = 0;
	for (size_t i = 0; i < 4; i++)
	{
		std::string num_Lay = patch::to_string((i + 1) * 5);
		for (size_t j = 0; j < 4; j++)
		{
			string num_H = patch::to_string((j + 1) * 5);
			if (j == 0)     // NOTICE: This is because of file names!
			{
				num_H = patch::to_string(0) + num_H;
			}
			for (size_t k = 0; k < 5; k++)
			{
				std::string sample = patch::to_string(k + 1);
				dataset_names[dataset_counter] = "Attack" + num_Lay + num_H + "Repeat" + sample + ".txt";
				dataset_counter++;
			}
		}
	}

	return dataset_names;
}

int main() {
	//OUTPUT
	ofstream MyExcelFile;
	string output_name = "Results_threads4_all.csv";
	MyExcelFile.open(output_name.c_str());

	vector <string> datasetss = Inputnames();
	for (size_t allDatasets = 0; allDatasets < 80; allDatasets++)
	{
		

		srand(time(0)); // !!!!!! It is necessary for generating different random variable!!!

		int Layer = 20;
		int H = 20;
		int ux = 70;
		int uy = 50; // Upper bounds for x and y
		int NumberOfAttacks = 3;
		int NumberofDivided = 1;

		//Uncomment the below out of the loop to generate new dataset
		/*cout << "Number of layer?  ";
		cin >> Layer;
		cout << "Number of nodes in each layer?  ";
		cin >> H;
		cout << endl;*/

		//cout << "Define and Upper bound for x's:  ";
		//cin >> ux;
		//cout << "Define and Upper bound for y's:  ";
		//cin >> uy;

		/*cout << "Number of attacks?  ";
		cin >> NumberOfAttacks;*/


		int n = Layer * H + 2; //number of nodes
		int m = (Layer - 1) * ((H - 2) * 3 + 4) + (H - 1) * Layer * 2 + (2 * H); //number of arcs

		cout << "n = " << n << "   m = " << m << endl;

		// Generating Network ------------------------------------------------------
		vector < vector <double> > MYnetworkXY(Layer + 3);
		MYnetworkXY = GeneratingNetwork(n, Layer, H, ux, uy);

		vector <double> x;
		for (size_t i = 0; i < MYnetworkXY[0].size(); i++)
			x.push_back(MYnetworkXY[0][i]);

		vector <vector <double> > y(Layer + 2);
		for (size_t i = 1; i < MYnetworkXY.size(); i++)
			for (size_t j = 0; j < MYnetworkXY[i].size(); j++)
				y[i - 1].push_back(MYnetworkXY[i][j]);
		// ----------------------------------------------------------------------------

		//Datas from file-----------------------------
		cout << "Do you want to input data from text file? (y/n) ";
		char textfile = 'y';		//change this if you would like to create new dataset
		//cin >> textfile;
		if (textfile == 'y')
		{
			cout << "Which file do you want to open? ";
			string Filename_input = datasetss[allDatasets];
			//cin >> Filename_input;
			cout << endl;
			string line;
			int counterXX = 0;

			ifstream file(Filename_input);
			if (file.is_open())
			{
				cout << Filename_input << endl;

				getline(file, line);
				getline(file, line);
				Layer = stoi(line);
				getline(file, line);
				getline(file, line);
				H = stoi(line);
				getline(file, line);
				getline(file, line);
				ux = stoi(line);
				getline(file, line);
				getline(file, line);
				uy = stoi(line);
				getline(file, line);
				getline(file, line);
				n = stoi(line);
				getline(file, line);
				getline(file, line);
				m = stoi(line);
				getline(file, line);
				for (size_t i = 0; i < n; i++)
				{
					getline(file, line);
					x[counterXX] = stod(line);
					counterXX++;
				}
				getline(file, line);
				getline(file, line);
				y[0][0] = stod(line);
				for (size_t i = 1; i < Layer + 1; i++)
				{
					for (size_t j = 0; j < H; j++)
					{
						getline(file, line);
						y[i][j] = stod(line);
					}
				}
				getline(file, line);
				y[Layer + 1][0] = stod(line);
				
				file.close();
			}
			else
			{
				cout << "file is not open" << endl;
			}
		}

		// Output results in text --------------------
		cout << "Do you want to save data to a text file? (y/n) ";
		char textfile1 = 'n';			//change this if you would like to save new dataset
		//cin >> textfile1;
		if (textfile1 == 'y') {
			cout << "Put a name for your text file: ";
			string Filename_save;
			cin >> Filename_save;
			cout << endl;
			ofstream myfile(Filename_save);
			if (myfile.is_open())
			{
				myfile << "Layer\n";
				myfile << Layer << endl;
				myfile << "H\n";
				myfile << H << endl;
				myfile << "Ux\n";
				myfile << ux << endl;
				myfile << "Uy\n";
				myfile << uy << endl;
				myfile << "n\n";
				myfile << n << endl;
				myfile << "m\n";
				myfile << m << endl;
				myfile << "x\n";
				for (size_t i = 0; i < n; i++)
				{
					myfile << x[i] << endl;
				}
				myfile << "y\n";
				myfile << y[0][0] << endl;
				for (size_t i = 1; i < Layer + 1; i++)
				{
					for (size_t j = 0; j < H; j++)
					{
						myfile << y[i][j] << endl;
					}
				}
				myfile << y[Layer + 1][0] << endl;
				myfile << "end";
				myfile.close();
			}
			else cout << "Unable to open file" << endl;
		}

		//BEGINNING TIME--------------------------------------------
		auto begin = chrono::steady_clock::now();

		// Printing -------------------------------------------------------------------
		cout << endl << "Ordered Pair:" << endl;
		PrintNetwork(x, y, Layer, H);
		// ----------------------------------------------------------------------------

		// Distances --------------------------------------------------	
		vector < vector <double> > CallC0(n + 1);
		CallC0 = InitialCosts(n, m, Layer, H, x, y);

		vector<double> C0a;
		for (size_t i = 0; i < m; i++)
			C0a.push_back(CallC0[n][i]);

		// First Cut --------------------------------------------------	
		vector <vector <double> > C;
		C.push_back(C0a);

		// ----------------------------------------------------------------------------
		// ----------------------------------------------------------------------------
		// OPTIMIZATION - MATHEMATICAL MODELING ! -------------------------------------
		//clock_t MIP_begin; //Defining variables for MIP TIME
		//clock_t MIP_end;
		//clock_t TimeLim;
		double Total_time = 0;

		//PARAMETERS -----------------------------------------------------------------
		vector <int> k(NumberOfAttacks);	 //# of Regions 
		k[0] = 5;
		k[1] = 7;
		k[2] = 9;

		vector <vector <double> > MYparameters = parameters(n, m, Layer, H, x, y, k[0]);

		//vector <vector <double>> B(n); //Adjacency Matrix
		//for (size_t i = 0; i < n; i++)
		//	B[i] = MYparameters[i];

		//vector <double> d = MYparameters[n]; //damage level
		vector <vector <double> > d(NumberOfAttacks);
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			for (size_t j = 0; j < k[t]; j++)
			{
				d[t].push_back((40.00 / pow(1.5, (j + 1))));
			}
		}

		double Ymax = MYparameters[n + 2][0];
		double Ymin = MYparameters[n + 2][1];
		double Qarea = MYparameters[n + 2][2];
		vector <double> BIGvalue = MYparameters[n + 5];

		//vector <double> q = MYparameters[n + 1]; //damage area
		vector <vector <double> > q(NumberOfAttacks);
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			q[t].push_back(-0.001);
			for (size_t j = 1; j < k[t] + 1; j++)
			{
				q[t].push_back((j)*Qarea / k[t]);
			}
		}

		vector <double> p1a = MYparameters[n + 3];
		vector <double> p2a = MYparameters[n + 4];

		vector <double> Alp1(NumberOfAttacks);
		vector <double> Alp2(NumberOfAttacks);

		// Initial Solutions -----------------------------------------------------------
		double Vforcheck = 0;
		for (size_t i = 0; i < m; i++) //better cutting palne - look at the obj of MIP
			Vforcheck += C0a[i] + NumberOfAttacks * d[0][0];
		double Wforcheck = 0;

		//----******************************************************************____________-------------------------------
		int MIP_iteration = 0;
		int iteration = 0;
		int flag = 0;

		//DEfining var and const. for MASTER Problem
		IloEnv env1;
		IloModel model1(env1);

		// VARIABLES ---------------------------------------------------------------
		IloNumVarArray Yvar(env1, m);
		for (size_t i = 0; i < m; i++)
			Yvar[i] = IloNumVar(env1, 0, 1);

		IloNumVar w(env1);
		model1.add(Yvar);		//Never remove it! why?! I do not know!!
		// CONSTRAINT --------------------------------------------------------------
		//for (size_t i = 1; i < n - 1; i++) // MAin Constraint!!!
		//{
		//	IloExpr MainConstraint(env1);
		//	for (size_t j = 0; j < m; j++)
		//	{
		//		MainConstraint += B[i][j] * Yvar[j];
		//	}
		//	model1.add(MainConstraint == 0);
		//	MainConstraint.end();
		//}
		int arcs_counter = 0;
		vector <vector <int> > adj_list(n);
		vector <vector <int> > arc_adj_list(n);

		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n; j++)
			{
				if (CallC0[i][j] > 0)
				{
					adj_list[i].push_back(j);
					arc_adj_list[i].push_back(arcs_counter);
					arcs_counter++;
				}
			}
		}

		if (arcs_counter != m)
		{
			cout << "Error in numbering!" << m << "  " << arcs_counter << endl;
			int er_numb;
			cin >> er_numb;
		}

		IloExprArray MainConstraint(env1, n);
		for (size_t i = 0; i < n; i++)
		{
			MainConstraint[i] = IloExpr(env1);
		}
		for (size_t i = 0; i < adj_list.size(); i++)
		{
			for (size_t j = 0; j < adj_list[i].size(); j++)
			{
				MainConstraint[i] += Yvar[arc_adj_list[i][j]];
				MainConstraint[adj_list[i][j]] += -Yvar[arc_adj_list[i][j]];
			}
		}
		model1.add(MainConstraint[0] == 1);
		for (size_t i = 1; i < n - 1; i++)
		{
			model1.add(MainConstraint[i] == 0);
		}
		model1.add(MainConstraint[n - 1] == -1);

		MainConstraint.end();

		//IloExpr MainConstraint0(env1);
		//IloExpr MainConstraint00(env1);
		//for (size_t i = 0; i < m; i++)
		//{
		//	MainConstraint0 += B[0][i] * Yvar[i];
		//	MainConstraint00 += B[n - 1][i] * Yvar[i];
		//}
		//model1.add(MainConstraint0 == 1);
		//model1.add(MainConstraint00 == -1);

		//MainConstraint0.end();
		//MainConstraint00.end();

		IloExpr Objective(env1);

		model1.add(IloMinimize(env1, w));
		IloCplex cplex1(model1);

		while (flag != 1)
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
			vector <double> LP;

			for (size_t i = 0; i < m; i++) {
				LP.push_back(cplex1.getValue(Yvar[i]));
			}

			LP.push_back(cplex1.getObjValue());

			Wforcheck = LP[m];
			cout << "The objective value (MAIN) = " << Wforcheck << endl;
			cout << endl;

			//----- Getting Y for next MIP
			//----- Recognizing Indicies of active arcs
			vector<double> Yoptimal;
			vector<int> activeArcs;

			for (size_t i = 0; i < m; i++)
			{
				Yoptimal.push_back(LP[i]);
				if (Yoptimal[i] > 0)
				{
					activeArcs.push_back(i);
				}
				//cout << "y" << i << "  " << Yoptimal[i] << "     ";
			}
			//cout << endl;

			// Heuristics -----------------------------------------------------------------------------
			int MIP_needed = 1;


			if (MIP_needed == 1)
			{
				vector <double> C_added(activeArcs.size());
				for (size_t i = 0; i < activeArcs.size(); i++)
				{
					C_added[i] = C0a[activeArcs[i]];
				}
				NumberofDivided++;

				vector <double> vals(activeArcs.size());
				// SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP --- SUBPROBLEM MIP ---
				for (size_t t = 0; t < NumberOfAttacks; t++)
				{
					IloEnv env;
					IloModel model(env);

					//Varibales -------------------------------------------------------------------
					typedef IloArray<IloNumVarArray> MatrixLambda;
					//typedef IloArray<FloatMatrix> Float3Matrix;
					MatrixLambda Lambda(env, activeArcs.size());

					IloNumVarArray Ca(env, activeArcs.size());

					IloNumVar Alpha1(env, 0);
					IloNumVar Alpha2(env, 0);

					IloNumVarArray BETAminus1(env, activeArcs.size());
					IloNumVarArray BETAplus1(env, activeArcs.size());
					IloNumVarArray BETAminus2(env, activeArcs.size());
					IloNumVarArray BETAplus2(env, activeArcs.size());

					for (size_t i = 0; i < activeArcs.size(); i++) {
						Lambda[i] = IloNumVarArray(env, k[t], 0, 1, ILOINT);

						BETAminus1[i] = IloNumVar(env, 0, IloInfinity);
						BETAplus1[i] = IloNumVar(env, 0, IloInfinity);
						BETAminus2[i] = IloNumVar(env, 0, IloInfinity);
						BETAplus2[i] = IloNumVar(env, 0, IloInfinity);

						Ca[i] = IloNumVar(env);
					}

					//constraints ------------------------------------------------------------------	
					for (size_t i = 0; i < activeArcs.size(); i++)
					{
						IloExpr constraint4a(env); // Constraint 4a ---------------------------
						constraint4a += C_added[i];
						IloExpr constraint4b(env); // Constraint 4b ---------------------------
						for (size_t j = 0; j < k[t]; j++) {
							constraint4a += d[t][j] * Lambda[i][j];
							constraint4b += Lambda[i][j];
						}

						model.add(constraint4b == 1);
						constraint4b.end();

						model.add(constraint4a == Ca[i]);
						constraint4a.end();

						model.add(p1a[activeArcs[i]] + BETAminus1[i] - BETAplus1[i] == Alpha1); // Constraint 4c ------------------------
						model.add(p2a[activeArcs[i]] + BETAminus2[i] - BETAplus2[i] == Alpha2); // Constraint 4d ------------------------

						for (size_t j = 0; j < k[t]; j++) // Constraint 4e ------------------------
							model.add((BETAminus1[i] + BETAminus2[i] + BETAplus1[i] + BETAplus2[i]) - BIGvalue[activeArcs[i]] * (1 - Lambda[i][j]) <= q[t][j + 1]);
					}


					//Objective -------------------------------------------------------------------
					IloExpr obj(env);
					for (size_t i = 0; i < activeArcs.size(); i++)
						obj += Ca[i] * Yoptimal[activeArcs[i]];

					model.add(IloMaximize(env, obj));

					auto MIP_begin = chrono::steady_clock::now();

					//SOLVE ----------------------------------------------------------------------
					IloCplex cplex0(model);
					cplex0.setParam(IloCplex::Param::Threads, 4);
					cplex0.setParam(IloCplex::TiLim, 600); //TIMElimit of solving
					cplex0.solve();
					MIP_iteration++;

					//TIME OF SOLVING-------------------------------------------------------------
					chrono::duration <double> elapsed_secs_MIP = chrono::steady_clock::now() - MIP_begin;
					Total_time += elapsed_secs_MIP.count();

					//COut -----------------------------------------------------------------------
					cout << "The objective value (MIP) = " << cplex0.getObjValue() << "  " << cplex0.getStatus() << endl;

					Vforcheck = cplex0.getObjValue();

					for (size_t i = 0; i < activeArcs.size(); i++)
					{
						vals[i] = cplex0.getValue(Ca[i]);
						//cout << "Ca" << i << " : " << vals[i] << "   " << CAsMIPStart[i] << "   " << C_H[activeArcs[i]] <<endl;
						C_added[i] = vals[i];
					}

					////Record -----------------------------------------------------------------------
					Alp1[t] = cplex0.getValue(Alpha1);
					Alp2[t] = cplex0.getValue(Alpha2);
					cout << "Alpha1_" << t << ": " << Alp1[t] << ", Alpha2_" << t << ": " << Alp2[t] << endl;

					env.end();
				}

				C_added.clear();
				//----------------------------------------------------------------------------	

				//adding new C to MAIN problem constraint
				vector<double> addingC = MYaddingC(m, k, NumberOfAttacks, vals, Yoptimal, p1a, p2a, Alp1, Alp2, d, q, C0a);


				//-------------------------------------------------------
				C.push_back(addingC);             //adding the last C to LP constraint
				addingC.clear();
				activeArcs.clear();
				//DamagedArcsMatrix.clear();
				//------------------------------------------------------------------------------


				if (Wforcheck + 0.1 >= Vforcheck)  //Stop criteria
					flag = 1;
			}

			//time limit:
			chrono::duration <double> TL = chrono::steady_clock::now() - begin;
			/*TimeLim = clock();
			double TL = double(TimeLim - begin) / CLOCKS_PER_SEC;*/
			cout << "Elaspsed Time: " << TL.count() << "sec." << endl;
			if (TL.count() > 3600)
				flag = 1;

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

		/*clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;*/
		chrono::duration <double> elapsed_secs = chrono::steady_clock::now() - begin;
		cout << endl << "Total time: " << elapsed_secs.count() << " sec." << endl;
		cout << "Time to solve MIPs: " << Total_time << " sec." << endl;

		// EXCEL OUTPUT for "SINGLE" attack----------------------------------------------------------------
		
		MyExcelFile << Wforcheck << ",";
		MyExcelFile << Vforcheck << ",";
		MyExcelFile << Vforcheck - Wforcheck << ",";
		MyExcelFile << iteration << ",";
		MyExcelFile << MIP_iteration << ",";
		MyExcelFile << elapsed_secs.count() << ",";
		MyExcelFile << Total_time << ",";
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			MyExcelFile << Alp1[t] << ",";
			MyExcelFile << Alp2[t] << ",";
		}
		MyExcelFile << endl;
		

		// Checking -------------------------------------------------------------------
		/*ofstream Checking;
		Checking.open("AA-Check.csv");
		for (size_t i = 0; i < m; i++)
		{
			Checking << p1a[i] << "," << p2a[i] << "," << C0a[i] << endl;
		}
		Checking << endl;
		for (size_t i = 0; i < NumberOfAttacks; i++)
		{
			Checking << Alp1[i] << "," << Alp2[i] << endl;
		}
		Checking << endl << "q,d" << endl;
		for (size_t t = 0; t < NumberOfAttacks; t++)
		{
			for (size_t i = 0; i < q.size(); i++)
			{
				Checking << q[t][i] << "," << d[t][i] << endl;
			}
		}

		Checking << endl;
		Checking.close();*/
	}
	MyExcelFile.close();
	// ---------------- // ----------------- END-----------------------------------
	return 0;
}