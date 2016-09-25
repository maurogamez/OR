#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <queue>
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <typeinfo>
#include <limits.h>

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>

ILOSTLBEGIN

typedef IloArray<IloNumArray> Num2DMatrix;
typedef IloArray<Num2DMatrix> Num3DMatrix;
typedef IloArray<IloNumVarArray> NumVar2DMatrix;
typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;

#define M INT_MAX
#define RC_EPS 1.0e-6

int main(IloInt argc, char **argv)
{
	IloEnv env;
	
	//Parameters Read from file 
	IloInt max_width;
	IloIntArray item(env), demand(env);
	
	try	
	{
		/*Reading data file*/
		const char* data_filename  = "cutstock.dat";
		
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
		else
		{
			cout<<"Reading"<<endl;
			datafile>>max_width;
			datafile>>item>>demand;
		}
		cout<<"Reading Done"<<endl;

		//Probems
		IloModel CutStock(env);
		IloModel PatGen(env);
	
		IloInt size, pat;
		size = item.getSize();

		//Variables
		IloNumVarArray X(env);			//Number of sheets of a particular pattern for CutStock ,0,IloInfinity
	
		//IloCplex cutSolver(env);
		//IloCplex patSolver(env);

		IloNumArray price(env, size);
		IloNumArray newPatt(env, size);
		
		/*Problem Cutting stock*/
		
		IloObjective totalStocks = IloMinimize(env);
		CutStock.add(totalStocks);

		IloRangeArray constr = IloAdd(CutStock, IloRangeArray(env, demand, IloInfinity));
		for (IloInt j = 0; j < size; j++) 
		{
			X.add(IloNumVar(totalStocks(1) + constr[j](int(max_width / item[j]))));
		}
		IloCplex cutSolver(CutStock);
		//cutSolver.extract(CutStock);
		cout<<"Cutting Defined"<<endl;
		/*Problem Pattern Generation*/
		
		IloObjective reducedCost = IloAdd(PatGen, IloMinimize(env, 1));
		IloNumVarArray Y(env, size, 0, IloInfinity, ILOINT);			//Number of sheets of a particular pattern for PatGen
		PatGen.add(IloScalProd(item, Y) <= max_width);
		IloCplex patSolver(PatGen);
		//patSolver.extract(PatGen);
		cout<<"Pattern defined"<<endl;

		/*Column Generation for addition of columns(iterations)*/
		
		while(1)
		{
			//cout<<"Start solve"<<endl;
			cutSolver.solve();
			//cout<<"solved once"<<endl;
			for (IloInt i = 0; i < size; i++)
			{
			   price[i] = -cutSolver.getDual(constr[i]);
			}
			reducedCost.setLinearCoefs(Y, price);
			patSolver.solve();
			if (patSolver.getValue(reducedCost) > -RC_EPS) break;
			patSolver.getValues(newPatt, Y);
			X.add( IloNumVar(totalStocks(1) + constr(newPatt)) );
		}
		//CutStock.add(IloConversion(env, X, ILOINT));
		//cutSolver.solve();
		cout << "Best solution uses " 
        << cutSolver.getObjValue() << " rolls" << endl;
	   cout << endl;
	   for (IloInt j = 0; j < X.getSize(); j++)
	   {
		   cout << "  Cut" << j << " = " << cutSolver.getValue(X[j]) << endl;
	   }
	   
	}//End try block
	catch (IloException& ex) 
	{
		cerr << "Error: " << ex << endl;
	}
	catch (...) 
	{
		cerr << "Error" << endl;
	}
	
	getch();
	return 0;
}