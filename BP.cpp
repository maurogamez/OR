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

IloEnv env;

IloObjective   obj;
IloObjective totalStocks;
IloNumVarArray X(env); //Number of sheets of a particular pattern for CutStock
//IloNumVarArray var(env);
IloNumArray vals(env);
IloRangeArray constr(env); 
IloCplex cutSolver(env);
IloCplex cplex(env);

IloNum eps;


IloNum incumbent;
IloNumArray intSol(env);
int min_obj=1;
IloInt size;

//Parameters Read from file 
IloInt max_width;
IloIntArray item(env), demand(env);
long int ctr=0;
void branch(queue<IloModel> mod, IloNum &incumbent,IloCplex cplex, IloNumVarArray vars);
void solve(IloModel);
int res;

int main(IloInt argc, char **argv)
{	
	//clock_t t;
	//t = clock();
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
			datafile>>max_width;
			datafile>>item>>demand;
		}
		env.out()<<"Reading Done"<<endl;
		size = item.getSize();
		//Probems
		IloModel CutStock(env);
		
		
		/*Problem Cutting stock*/
		
		totalStocks = IloMinimize(env);
		CutStock.add(totalStocks);

		constr = IloAdd(CutStock, IloRangeArray(env, demand, IloInfinity));
		for (IloInt j = 0; j < size; j++) 
		{
			X.add(IloNumVar(totalStocks(1) + constr[j](int(max_width / item[j]))));
		}
		
		
		env.out()<<"Cutting Defined"<<endl;
		solve(CutStock);
		env.out() << "Best solution uses " << cutSolver.getObjValue() << " rolls" << endl;
		env.out() << endl;
		/*Branch and bound*/
		if(vals.areElementsInteger())
		{
			
			env.out()<<"\n\nFinal incumbent solution: "<<totalStocks<<endl;
			//env.out()<<"Final incumbent solution vector: "<<vals<<endl;
			//t = clock() - t;
			//env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
			getch();
			//env.out() << "Solution vector = " << vals << endl;
			//getch();
			return 0;
		}
		
		eps = cplex.getParam(IloCplex::EpInt);
		
		/*Queue of problems*/
		queue<IloModel> mod;
		mod.push(CutStock);
		if(min_obj==1)
			incumbent = M;
		else incumbent = -M;
		//var.add(X);
		//cplex.getValues(vals,var);
		
		//env.out()<<vals<<endl;
		cutSolver.exportModel("relaxed.lp");
		getch();
		branch(mod, incumbent, cplex, X);
		//env.out()<<incumbent<<endl;
		
		env.out()<<"---------------------------------------------"<<endl;
		env.out()<<"\n\nFinal incumbent solution: "<<incumbent<<endl;
		env.out()<<"Final incumbent solution vector: "<<intSol<<endl;
		//t = clock() - t;
		//env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
		getch();


	}//End try block
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


void solve(IloModel CutStock)
{
	cutSolver.extract(CutStock);
	//IloNumArray temp(env);
	//cutSolver.getValues(temp,X);
	
	/*Problem Pattern Generation*/
	IloModel PatGen(env);
	IloObjective reducedCost = IloAdd(PatGen, IloMinimize(env, 1));
	IloNumVarArray Y(env, size, 0, IloInfinity, ILOINT);			//Number of sheets of a particular pattern for PatGen
	PatGen.add(IloScalProd(item, Y) <= max_width);
		
	IloNumArray price(env, size);
	IloNumArray newPatt(env, size);
		
	IloCplex patSolver(PatGen);
	//patSolver.extract(PatGen);
	cout<<"Pattern defined"<<endl;

	/*Column Generation for addition of columns(iterations)*/
		
	while(1)
	{
		//cout<<"Start solve"<<endl;
		if(cutSolver.solve())
			res=1;
		else res=0;
		//cout<<"solved once"<<endl;
		for (IloInt i = 0; i < size; i++)
		{
			price[i] = -cutSolver.getDual(constr[i]);
		}
		reducedCost.setLinearCoefs(Y, price);
		if(patSolver.solve())
			res=1;
		else res=0;
		if (patSolver.getValue(reducedCost) > -RC_EPS) break;
		patSolver.getValues(newPatt, Y);
		X.add( IloNumVar(totalStocks(1) + constr(newPatt)) );
	}
	//CutStock.add(IloConversion(env, X, ILOINT));
	//cutSolver.solve();
	//cout << "Best solution uses " 
    //<< cutSolver.getObjValue() << " rolls" << endl;
	cout << endl;
	cutSolver.getValues(vals,X);
	env.out()<<"Vals is here: "<<vals<<endl;
	
	for (IloInt j = 0; j < X.getSize(); j++)
	{
		cout << "  Cut" << j << " = " << cutSolver.getValue(X[j]) << endl;
	}
	
}
void branch(queue<IloModel> mod, IloNum &incumbent, IloCplex cplex, IloNumVarArray var)
{		
	/*Load current model(branch) with non-Integer solutions*/
	string left = "left", right="right";
	ctr++;
	left=left+to_string(static_cast<long long>(ctr))+".lp";
	right=right+to_string(static_cast<long long>(ctr))+".lp";
	//cout<<left<<endl;
	if(mod.size()==0)
		return;
	IloModel model(mod.front());
	mod.pop();
	//cplex.extract(model);
	//cplex.solve();
	solve(model);
	//IloNumArray vals(env);
	//cplex.getValues(vals, var);
	//cutSolver.exportModel("tr.lp");
	env.out()<<vals<<endl;
	//getch();
	/*Locate position of non-integrality*/
	IloInt i;
	IloNum value;
	for( i=0;i<vals.getSize();i++)
	{
		if(abs(vals[i]-floor(vals[i]))>=eps)
			break;
	}
	value=vals[i];

	IloNum objVal;

	/*Branch Left->lower bound*/
	IloModel lb(model);	
	IloConstraint con1 = var[i]<=floor(value);
	lb.add(con1);			
	//cplex.extract(lb);
	solve(lb);
	cutSolver.exportModel(left.c_str());
	//patSolver.exportModel(left.c_str());
	env.out()<<"\nleft BRANCHED PROBLEM"<<endl;
	env.out()<<"-----------------------------------------------------------"<<endl;
	cout<<"HEY"<<endl;
	//cplex.exportModel(left);
	//env.out()<<lb<<endl;
	if(res)
	{
		mod.push(lb);
		//cplex.getValues(vals,var);
		objVal = cutSolver.getObjValue();
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << objVal << endl;
		//env.out() << "Solution vector = " << vals << endl;
		getch();
	
		if(vals.areElementsInteger())
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
			{
				incumbent = objVal;
				//cplex.getValues(intSol,var);
				//env.out()<<vals;
			}
		}
		else if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
		{
			branch(mod,incumbent,cplex,var);
		}
	}
	else env.out()<<"Lower Bound Infeasible!"<<endl;
	
	/*Branch right->upper bound*/
	IloModel ub(model);
	IloConstraint con2 = var[i]>=ceil(value);
	ub.add(con2);
	ub.remove(con1);
	//cplex.extract(ub);
	solve(ub);
	env.out()<<"\nright BRANCHED PROBLEM"<<endl;
	env.out()<<"-----------------------------------------------------------"<<endl;
	cout<<"YO"<<endl;
	cutSolver.exportModel(right.c_str());
	//env.out()<<ub<<endl;
	//cplex.exportModel("out2.lp");
	if(res)
	{
		mod.push(ub);
		//cplex.getValues(vals,var);
		objVal = cutSolver.getObjValue();;
		//env.out() << "Solution status = " << objVal << endl;
		env.out() << "Solution value  = " << cutSolver.getObjValue() << endl;
		//env.out() << "Solution vector = " << vals << endl;
		getch();
		if(vals.areElementsInteger())
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
			{
				incumbent = objVal;
				//cplex.getValues(intSol,var);
				//env.out()<<vals;
			}
		}
		else if((objVal<=incumbent && min_obj ==1) || (objVal>=incumbent && min_obj ==0))
		{
			branch(mod,incumbent,cplex,var);
		}
	}
	else env.out()<<"Upper Bound Infeasible!"<<endl;
	ub.remove(con2);
	env.out()<<"Incum: "<<incumbent<<endl;
	return;
}