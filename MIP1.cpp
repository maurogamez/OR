/************************************************************************************
 * This example is an example on IP from Application oriented Guide to LR by Fisher, Interfaces, 1985.
 * Minimize 5x - 3y
 * s.t.: 
 *			x + 2y >= 10
 *			 2x - y >= 0
 *			  x - 3y >= -13
 * 0 <= y <= 10 integers
 ***********************************************************************************/
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <deque>
#include <sstream>
#include <time.h>
//#include <iomanip.h>
#include <stdlib.h>
#include <vector>//for vectors
#include <math.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

using namespace std;

ILOSTLBEGIN

//typedef IloArray<IloNumArray> TwoDMatrix;

int main(int argc, char **argv) 
{
	IloEnv env;
	try 
	{
		IloNumVar X(env, 0, IloInfinity, ILOFLOAT);
		IloNum X_val;
		IloNumVar Y(env, 0, 10, ILOINT);
		IloNum Y_val;
		IloModel model(env);
		IloExpr Obj(env);
		Obj = 5*X - 3*Y;
		model.add(IloMinimize(env,Obj)); // IloMinimize is used for minimization problems
		//Alternatively, model.add(IloMinimize(env,Obj)); can be replaced by the following three lines.
		//This is useful while iterating in Decomposition Techniques where the objective function is redefined at each iteration
		//IloObjective Objective = IloMinimize(env); 
		//model.add(Objective); 
		//Objective.setExpr(Obj);

		Obj.end();
		model.add(X + 2*Y >= 10);
		model.add(2*X - Y >= 0);
		model.add(X - 3*Y >= -13);
		// Optimize
		IloCplex cplex(model);
		//cplex.setOut(env.getNullStream()); // This is to supress the output of Branch & Bound Tree on screen
		//cplex.setWarning(env.getNullStream()); //This is to supress warning messages on screen
		cplex.solve();//solving the MODEL
		if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
		{
			env.out() << "Problem Infeasible" << endl; 
		}
		X_val = cplex.getValue(X);
		Y_val = cplex.getValue(Y);
		// Print results
		cout<< "Objective Value = " << cplex.getObjValue() << endl;
		cout<<"X = "<<X_val<<endl;
		cout<<"Y = "<<Y_val<<endl;
	}
	catch(IloException &e) 
	{
		env.out() << "ERROR: " << e << endl;
	}
	catch(...)
	{
		env.out() << "Unknown exception" << endl;
	}
	env.end();
	getch();
	return 0;
}