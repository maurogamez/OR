#include <ilcplex/ilocplex.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	IloEnv env;
	IloModel model(env);
	try
	{
		ifstream file;
		file.open("Data_HS_CAB_5node.dat",ios::in);

	}
	catch(IloException &ex)
	{
		cerr << "Error: " << ex << endl;
	}
	catch(...)
	{
		cerr << "Error" << endl;
	}
	getch();
	return 0;
}
