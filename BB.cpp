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

IloEnv env;
IloObjective   obj;
IloNum eps;

IloNum incumbent;
IloNumArray intSol(env);
int min_obj;
//long int ctr=0;
void branch(queue<IloModel> mod, IloNum &incumbent,IloCplex cplex, IloNumVarArray vars);

int main(IloInt argc, char **argv)
{
	clock_t t;
	t = clock();
	char min[] = "Minimize", max[] = "Maximize";

	try
	{
		/*Read Master problem*/
		IloModel master(env);
        IloCplex cplex(env);
		
        IloNumVarArray var(env);
        IloRangeArray  rng(env);
		const char* data_filename  = "ModelHsSingleAlloc10Node.lp";
		if (argc > 1)
		{
			data_filename = argv[1];
		}
        cplex.importModel(master, data_filename, obj, var, rng);
		
		ifstream fin(data_filename);
		string line;
		min_obj=0;
		while (getline(fin, line))
		{
			if (line.find(min) != string::npos)
			{
				min_obj=1;
				break;
			}
			if (line.find(max) != string::npos)
			{
				break;
			}
		}
		env.out()<<"MASTER PROBLEM"<<endl;
		env.out()<<"-----------------------------------------------------------"<<endl;
		//env.out()<<master<<endl;
		/*Create Relaxed Master Problem*/
		IloModel relaxed(env);
		relaxed.add(master);
		relaxed.add(IloConversion(env, var, ILOFLOAT));
			
		/*Solve relaxed LP*/
		env.out()<<"\nRELAXED PROBLEM"<<endl;
		env.out()<<"-----------------------------------------------------------"<<endl;
		cplex.extract(relaxed);
		int result = cplex.solve();
		cplex.exportModel("Relaxed.lp");
		//env.out()<<relaxed<<endl;
		/*Check Infeasibility*/
		if ( !result ) {
           env.error() << "Infeasible solution" << endl;
           return 0;
        }
		IloNumArray vals(env);
		cplex.getValues(vals, var);
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		//env.out() << "Solution vector = " << vals << endl;
		//getch();
		/*Check Integrality*/
		if(vals.areElementsInteger())
		{
			
			env.out()<<"\n\nFinal incumbent solution: "<<cplex.getObjValue()<<endl;
			//env.out()<<"Final incumbent solution vector: "<<vals<<endl;
			t = clock() - t;
			env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
			getch();
			//env.out() << "Solution vector = " << vals << endl;
			//getch();
			return 0;
		}
		
		eps = cplex.getParam(IloCplex::EpInt);
		
		/*Queue of problems*/
		queue<IloModel> mod;
		mod.push(relaxed);
		if(min_obj==1)
			incumbent = M;
		else incumbent = -M;
		branch(mod, incumbent, cplex, var);
		//env.out()<<incumbent<<endl;
		
	}//end of of try block
	catch (IloException& ex) 
		{
			cerr << "Error: " << ex << endl;
		}
		catch (...) 
		{
			cerr << "Error" << endl;
		}
	env.out()<<"---------------------------------------------"<<endl;
	env.out()<<"\n\nFinal incumbent solution: "<<incumbent<<endl;
	env.out()<<"Final incumbent solution vector: "<<intSol<<endl;
	t = clock() - t;
	env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
	getch();
	return 0;
}
void branch(queue<IloModel> mod, IloNum &incumbent, IloCplex cplex, IloNumVarArray var)
{		
	/*Load current model(branch) with non-Integer solutions*/
	/*string left = "left", right="right";
	ctr++;
	left=left+to_string(static_cast<long long>(ctr))+".lp";
	right=right+to_string(static_cast<long long>(ctr))+".lp";*/
	//cout<<left<<endl;
	if(mod.size()==0)
		return;
	IloModel model(mod.front());
	mod.pop();
	cplex.extract(model);
	cplex.solve();
	IloNumArray vals(env);
	cplex.getValues(vals, var);
	//cplex.exportModel("tr.lp");
	
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
	cplex.extract(lb);
	//cplex.exportModel(left.c_str());
	env.out()<<"\nleft BRANCHED PROBLEM"<<endl;
	env.out()<<"-----------------------------------------------------------"<<endl;
	cout<<"HEY"<<endl;
	//cplex.exportModel(left);
	//env.out()<<lb<<endl;
	if(cplex.solve())
	{
		mod.push(lb);
		cplex.getValues(vals,var);
		objVal = cplex.getObjValue();
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << objVal << endl;
		//env.out() << "Solution vector = " << vals << endl;
		//getch();
	
		if(vals.areElementsInteger())
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
			{
				incumbent = objVal;
				cplex.getValues(intSol,var);
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
	cplex.extract(ub);
	env.out()<<"\nright BRANCHED PROBLEM"<<endl;
	env.out()<<"-----------------------------------------------------------"<<endl;
	cout<<"YO"<<endl;
	//cplex.exportModel(right.c_str());
	//env.out()<<ub<<endl;
	//cplex.exportModel("out2.lp");
	if(cplex.solve())
	{
		mod.push(ub);
		cplex.getValues(vals,var);
		objVal = cplex.getObjValue();
		//env.out() << "Solution status = " << objVal << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;
		//env.out() << "Solution vector = " << vals << endl;
		//getch();
		if(vals.areElementsInteger())
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
			{
				incumbent = objVal;
				cplex.getValues(intSol,var);
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
	
	return;
}		