// This is the code for Single Allocation Hub and Spoke Network Design for the model described in:
// Pirkul, H.,  Schilling, D.A. 1998. An efficient procedure for designing single allocation hub 
// and spoke systems. Management Science, 44(1), 235-242.

#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
//#include 

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

ILOSTLBEGIN

typedef IloArray<IloNumArray> TwoDMatrix;
typedef IloArray<TwoDMatrix> ThreeDMatrix;
typedef IloArray<ThreeDMatrix> FourDMatrix;

int main(int argc, char **argv)   
{
	IloEnv env;
	int N, n_cap;//N=no. of nodes, n_cap=n. of max capacity levels at a hub
	IloNum c_collect=1, c_distribute=1, c_tranship=0.90;//per unit costs of shipping

	TwoDMatrix lambda(env);//flow data between pairs of nodes
	TwoDMatrix dist(env);//distance data between pairs of nodes
	IloNum p = 3;//Number of hubs to open

	IloNum eps;//later assigned as eps = cplex.getParam(IloCplex::EpInt); EpInt is the tolerance gap for integer variables

	try 
	{
		///////// DATA FILE READING //////////
	
		//const char* data_filename  = "Data_HS_CAB_5node.dat";
		const char* data_filename  = "Data_HS_CAB_10node.dat";
		if (argc > 1)
		{
			data_filename = argv[1];
		}
		fstream datafile;
		datafile.open(data_filename,ios::in);

		if (!datafile) 
		{
			cerr << "ERROR: could not open file " << data_filename << " for reading" << endl;
			cerr << "usage:   " << argv[0] << " <datafile>" << endl;
			throw(-1);
		}

		datafile >> lambda >> dist;
		
		N = lambda.getSize();
		cout<<"Number of nodes = "<<N<<endl;
		
		IloBool consistentData = (lambda.getSize() == dist.getSize());
		if (!consistentData) 
		{
			cerr << "ERROR: data file '" << data_filename << "' contains inconsistent data" << endl;
			throw(-1);
		}
		datafile.close();

		// READING DONE..................................
		//=================Decleare Variables===============================
		IloModel model(env);
		typedef IloArray<IloBoolVarArray> array2d;//Creating a 2d array of x variables
		typedef IloArray<array2d> array3d;//Creating a 3d array of x variables
		typedef IloArray<array3d> array4d;//Creating a 4d array of x variables
		array4d x(env, N);//xijkm
		array2d z(env, N); //zik = 1 if node i is allocated to hub k; zkk = 1 means node k is a hub
		//=======================================================================
		FourDMatrix x_val(env, N);//to hold values of x
		TwoDMatrix z_val(env, N);//to hold values of z
		//=======================================================================
		for (int i=0; i<N; i++)
		{
			x[i]=array3d(env, N);
			x_val[i]=ThreeDMatrix(env, N);
			for (int j=0;j<N;j++)
			{ 
				x[i][j]=array2d(env, N);
				x_val[i][j]=TwoDMatrix(env, N);
				for (int k=0;k<N;k++)
				{
					x[i][j][k]=IloBoolVarArray(env, N);
					x_val[i][j][k]=IloNumArray(env, N);
				}
			}
		}

		for (int i=0;i<N;i++)
		{
			for (int k=0;k<N;k++) 
			{
				z[i] = IloBoolVarArray(env, N);
				z_val[i] = IloNumArray(env, N);
			}
		}

		// Objective Function: Minimize: sum {i in 1..N}{j in i..N}{k in 1..N}{m in 1..N}
		//lambda[i][j][k][m]*(c_collect*dist[i][k]+c_tranship*dist[k][m]+c_distribute*dist[m][j])*x[i][j][k][m]
		IloExpr Obj(env); // Creates an expression with the name Obj (Objective)
		for (int i=0;i<N;i++)
		{
			for (int j=0;j<N;j++)
			{
				for (int k=0;k<N;k++)
				{
					for (int m=0;m<N;m++)
					{
						Obj+=lambda[i][j]*(c_collect*dist[i][k]+c_tranship*dist[k][m]+c_distribute*dist[m][j])*x[i][j][k][m];
					}
				}
			}
		}


		// model.add is a function to add constraints and objectives in the CPLEX environment
		model.add(IloMinimize(env,Obj)); // IloMinimize is used for minimization problems
		Obj.end(); // Clear out the memory allocated for the expression 


		//Constraint 1: for {i in 1..N}: sum {k in 1..N}(z[i][k])=1;
		for(int i=0;i<N;i++)
		{
			IloExpr SumZ(env);
			for(int k=0;k<N;k++)
			{
				SumZ+=z[i][k];
			}
			model.add(SumZ==1);
			SumZ.end();
		}


		//Constraint 2: for {i in 1..N}{k in 1..N}: z[i][k] <= z[k][k] 
		for(int i=0;i<N;i++)
		{
			for(int k=0;k<N;k++)
			{
				model.add(z[i][k] <= z[k][k]);
			}
		}

		//Constraint 3: sum {k in 1..N}(z[k][k])=p;
		IloExpr SumZkk(env);
		for(int k=0;k<N;k++)
		{
			SumZkk+=z[k][k];
		}
		model.add(SumZkk==p);
		SumZkk.end();

		//Constraint 4: for {i in 1..N}{j in 1..N}{k in 1..N}:sum {m in 1..N}(x[i][j][k][m])=z[i][k];
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				for (int k=0;k<N;k++)
				{
					IloExpr SumX(env);
					for(int m=0;m<N;m++)
					{
						SumX+=x[i][j][k][m];
					}
					model.add(SumX==z[i][k]);
					SumX.end();
				}
			}
		}

		//Constraint 5: for {i in 1..N}{j in 1..N}{m in 1..N}:sum {k in 1..N}(x[i][j][k][m])=z[j][m];
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				for (int m=0;m<N;m++)
				{
					IloExpr SumX(env);
					for(int k=0;k<N;k++)
					{
						SumX+=x[i][j][k][m];
					}
					model.add(SumX==z[j][m]);
					SumX.end();
				}
			}
		}
		//==============================================================================

		//Optimize
		IloCplex cplex(model);
		cplex.exportModel("ModelHsSingleAlloc10Node.lp");
		//cplex.setOut(env.getNullStream()); //This is to supress the output of Branch & Bound Tree on screen
		//cplex.setWarning(env.getNullStream()); //This is to supress warning messages on screen
		eps = cplex.getParam(IloCplex::EpInt);

		cplex.solve();//solving the MODEL

		if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
		{
			cout << "Problem is Infeasible" << endl; 
		}

		// Print results
		cout << "Minimum Cost = " << cplex.getObjValue() << endl;

		cout<<"Nodes \t Hubs: "<<endl;

		for (int i=0;i<N;i++)
		{
			for (int k=0;k<N;k++)
			{
				if (cplex.getValue(z[i][k]) > 1-eps)
				{
					cout<<i+1<<"\t"<<k+1<<endl;
					z_val[i][k] = 1;
				}
				else
				{
					z_val[i][k] = 0;
				}
			}
		}

		cout<<"i j k m x[i][j][k][m]" << endl;
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				for (int k=0;k<N;k++)
				{
					for (int m=0;m<N;m++)
					{
						if (cplex.getValue(x[i][j][k][m]) > 1-eps)
						{
							x_val[i][j][k][m] = 1;
							cout<<i+1<<" "<<j+1<<" "<<k+1<<" "<<m+1<<" "<<cplex.getValue(x[i][j][k][m])<<endl;
						}
						else
						{
							x_val[i][j][k][m] = 0;
						}
					}
				}
			}
		}

		env.end();

		}//end of of try block
		catch (IloException& ex) 
		{
			cerr << "Error: " << ex << endl;
		}
		catch (...) 
		{
			cerr << "Error" << endl;
		}
		return 0;
}