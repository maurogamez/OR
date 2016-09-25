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

ILOSTLBEGIN

typedef IloArray<IloNumArray> Num2DMatrix;

typedef IloArray<Num2DMatrix> Num3DMatrix;//3D array of Num

typedef IloArray<IloNumVarArray> NumVar2DMatrix; // 2D array of NumVar

typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;//3D array of Var

int main(IloInt argc, char **argv) {

IloEnv env;

cout << "Hello World!"<<endl;

}