/*Hardcoded: S,R,M,N, node_number*/
/*Branch and Price for cross Dock problem*/

//STANDARD C++ LIBRARIES
#include<stdio.h>
#include<conio.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include <queue>		//For implementing branching
#include <deque>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <typeinfo>
#include <limits.h>		//For using INT_MAX

//CONCERT TECHNOLOGY LIBRARIES
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>
#include <ilconcert/ilocsvreader.h>

ILOSTLBEGIN


typedef IloArray<IloBoolArray> Bool2DMatrix;		//2D array of Bool
typedef IloArray<Bool2DMatrix> Bool3DMatrix;		//3D array of Bool
typedef IloArray<IloIntArray> Int2DMatrix;			//2D array of Int
typedef IloArray<Int2DMatrix> Int3DMatrix;			//3D array of Int
typedef IloArray<IloNumArray> Num2DMatrix;			//2D array of Num
typedef IloArray<Num2DMatrix> Num3DMatrix;			//3D array of Num
typedef IloArray<Num3DMatrix> Num4DMatrix;			//4D array of Num
typedef IloArray<IloBoolVarArray> BoolVar2DMatrix;	// 2D array of Bool Var
typedef IloArray<BoolVar2DMatrix> BoolVar3DMatrix;	//3D array of Var
typedef IloArray<IloIntVarArray> IntVar2DMatrix;	// 2D array of Int Var
typedef IloArray<IntVar2DMatrix> IntVar3DMatrix;	//3D array of Int Var
typedef IloArray<IloNumVarArray> NumVar2DMatrix;	// 2D array of Num Var
typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;	//3D array of Num Var
typedef IloArray<NumVar3DMatrix> NumVar4DMatrix;	//4D array of Num Var
typedef IloArray<IloRangeArray> Range2DMatrix;		//2D Arrays of Ranges


#define MAX INT_MAX
#define RC_EPS 1.0e-6		//if shadow_price or Reduced Cost < RC_EPS, treat it as 0
#define INF (0x7FFFFFFF)
#define verbose (1)
#define NUM_SOLUTIONS 100
#define TOL 1.0e-6			/*define tolerance level to check improving obj value*/
#define LN 99999999.0		/* Large number */
#define SETSIZE 50
#define feasibility (RC_EPS*LN)

IloEnv env;					//ENVIRONMENT DECLARED

/*STRUCTURE FOR STORING (r,j,k) FOR A PARTICULAR Y-COLUMN*/
struct YHEADER
{
	int planid;
	int r;
	int j;
} yheader[50];

/*STRUCTURE FOR RETRIEVAL OF Y-COLUMN VALUES i.e PLAN, FOR EACH (r,j)*/
struct YASSIGNMENT
{
	int nofplans; 
	int assignment[20][10]; //assignments for each plan
} yassignment[50][50]; 

/*STRUCTURE FOR STORING NUMBER OF SUPPLIERS AND THEIR IDS TO WHICH EACH r IS ASSOCIATED*/
struct FLOWDATA
{
	int numsuppliers, supplierids[20];
	double flow[20];
} flowretailer[50];

/*CLASS STORING DATA AT EACH NODE OF THE BRANCH AND PRICE TREE*/
class NODE_DATA  //template for data to be stored in each B & P node
{
public:
	int status; //status=1 for active, 0 for pruned
    double lowerbound;
    double upperbound;
    double dual[100];
	int xcolumnindices [100],ycolumnindices[100];				//indices of variables in basis
	int numxinclude,numxexclude, xexclude[100];					// xexclude: indices of the columns excluded
	IloNum xinclude[100];										//xinclude: actual column included
	int numyinclude,yinclude[100],numyexclude,yexclude[100];	//indices of variables to be included or excluded
	int is_xincluded;
	
	NODE_DATA()
	{
		numxinclude=0;
		numyinclude=0;
		numxexclude=0;
		numyexclude=0;
		is_xincluded=0;
		//xexclude = NULL;
	}

} node[1000];

/*USED FOR ASSIGNMENT*/
struct SET
{
	int status;
	double a[100][50],z;
	int assignment[50];
} set[SETSIZE];

/*STRUCTURE FOR STORING K-BEST SOLUTION ACCORDING TO MURTHY'S ALGORITHM*/
struct YSOLUTIONS
{
	double z;
	int assignment[50];
} ksolutions[100];

/*FUNCTION PROTOTYPE DECLARATIONS*/
int getbestsol4semiassnprob(int setnum, int nr, int nc);	// solve semiassignmentproblem
int setpartition(int setnum,int nr, int nc);				//partition from node, return index of best set
void cleanspace();

static void readData (const char* filename, Num2DMatrix& flow, Num2DMatrix& distance);
static void sort(IloEnv env, Num2DMatrix Unsorted_matrix, Num2DMatrix& Sorted_array, IloBool order);//Bubble Sort; 
// order = 1 for ascending; order = 0 for descending
static void report1 (IloCplex& cplex_master, IloNumVarArray Col_select, IloRangeArray Constraint_master, vector <int> Column_index_r, vector <int> Column_index_j, int Num_source, int Num_destination, int Num_idoors, int Num_odoors);
static void report2 (IloAlgorithm& cplex_sub, IloNumArray newCol_val, IloObjective Objective_sub);
static void report3 (IloCplex& cplex_master, IloNumVarArray Col_select);
static void gen_col_matrix(Num2DMatrix comb, int Num_source, int Num_destination, Int2DMatrix SI, Int2DMatrix RJ, Int2DMatrix IJ, Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix);
static void gen_cols(Num2DMatrix comb, int Num_source, int Num_destination, Int2DMatrix SI, Int2DMatrix RJ, Int2DMatrix IJ, Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix);
void Millers_Code(int S, int R, Num2DMatrix kbestsols, IloNumArray koptvals, IloNumArray shadow_price);

/*GLOBAL VARIABLE DECLARATIONS*/
int setsgenerated,setsremaining;
double ycolumns_semiassignment[30][30];
int basissize, numycolumns, numnodegen, numnodeactive;
double upperbound, lowerbound;

/*
const int S=4;	//STORES TOTAL NUMBER OF SUPPLIERS
const int R=3;	//STORES TOTAL NUMBER OF RETAILERS
const int N=3;	//STORES TOTAL NUMBER OF OUTBOUND DOORS
const int M=4;	//STORES TOTAL NUMBER OF INBOUND DOORS
*/
int S=3, R=2, M=3, N=2;
int i,s;

IloModel model_master(env);
IloNumVarArray Col_select(env);
IloObjective Col_cost = IloAdd(model_master, IloMinimize(env));
IloRangeArray Constraint_master(env);
IloCplex cplex_master(model_master);

IloModel model_sub_y(env);
//IloBoolVarArray ycolumns(env);
IloObjective Objective_sub_y = IloMinimize(env);
IloRangeArray sub_y_constraint(env);//Constraints added later according to r
IloCplex cplex_sub_y(model_sub_y);

IloModel model_sub_x(env);
//IloBoolVarArray xcolumns(env);
IloObjective Objective_sub_x = IloMinimize(env);
IloRangeArray sub_x_constraint(env);
IloCplex cplex_sub_x(model_sub_x);

IloNumArray reduced_cost_y_column(env, 100);
vector <int> Column_index_r;//This is to map the retailer (r) associated with the column added; r=Num_destination for x column
vector <int> Column_index_j;//This is to map the outbound door (j) associated with the the column added; j=Num_odoors for x column
//vector <int> Column_index_plan;//This is to map the plan (k) associated with the the column added;
IloNumArray shadow_price(env, R+(S*M)+1);//To store shadow prices
IloNumArray ycolumns_val(env);//To store the new column (y or x) obtained
IloNumArray xcolumns_val(env);//To store the new column (y or x) obtained
int isY[10000];					//To check if column added is a y column(=1) or an x column(=0)
IloNum xcol[1000][50];		  //Stores all the xcolumns added into the basis
int xcol_ctr=0;			  //Number of xcolumns added
int xcol_check=-1;			  //To check which which xcol is fractional
IloNumArray obj_init(env,100);		//To store the objective function values to be added in the initial columns

IloNum eps;					//value below which cplex considers number to be 0
IloNumArray vals(env);		//for storing value of the solution vector
IloNum RC_FLAG = 1;			//Flag indicating negative reduced cost for at least one y or x column
int Col_Counter=3;
int iter=0;					//Forces x_column to be added only once if branched right	  
int frac=0;				//Index of fractional solution

IloNum incumbent;			//Stores the incumbent integer value uptill that point in the code
IloNumArray intSol(env);	//Final incumbent solution is stored here
int min_obj=1;				//This is to imply the branch and bound works for minimization problem
int size;
int res=0;					//If the sove function is infeasible=0 else =1
int nodenum = 0;			//Keeps count of the nodes in the branches

int ctr=0;	//For debugging

//Parameters Read from file
Num2DMatrix flow(env);
Num2DMatrix dist(env);


void generateycolumn(int nodenum, int S, int M, int N, int R, int r, int j, Num2DMatrix dist, IloNumArray reduced_cost_y_column); //generate y column by solving assignment subproblem 1. This is only a definition.
void branch(queue<IloModel> mod, IloNum &incumbent,IloCplex &cplex, IloNumVarArray vars);	//This function creates branching at each node
void solve(IloCplex &cplex_master); //res=1 implies solution obtained whereas res=0 implies infeasible
//This is the custom solve function that implements the column generation

using namespace::std;

int main(int argc, char **argv)
{	
	//clock_t t;
	//t = clock();
	try	
	{
		/*-------------------Reading data file-------------------------*/
		//const char* data_filename  = "cutstock.dat";		
		if ( argc > 1 )
			readData(argv[1], flow, dist);
		else
			readData("Ex_Illusttrative.txt", flow, dist);

		cout<<flow<<endl;
		cout<<dist<<endl;
		//S=4;
		/*S = 4;//flow.getSize();
		R = 3;//flow[0].getSize();
		M = 4;//dist.getSize();
		N = 3;//dist[0].getSize();*/

		cout<<"No. of Sources = "<<S<<endl;
		cout<<"No. of Destinations = "<<R<<endl;
		cout<<"No. of Inbound Doors = "<<M<<endl;
		cout<<"No. of Outbound Doors = "<<N<<endl;
				
		env.out()<<"Reading Done"<<endl;
		//system("pause");
		/*---------------------------Sorting of Flow and distance matrices-------------------------------------*/
		cout<<"\nSorting of Flow and distance matrices "<<endl;
		Num2DMatrix flow_des_sort(env, S*R);		//flow matrix as a sorted array in descending order of flows
		for (int s=0; s<S*R; s++)
		{
			flow_des_sort[s] = IloNumArray(env, 3);		//To store flow, source, destination
		}
		sort(env, flow, flow_des_sort, 0);
		cout<<"Flow Array SORTED in Decending Order"<<flow_des_sort<<endl;

		Num2DMatrix distance_des_sort(env, M*N);		//distance matrix as a sorted array in descending order of flows
		for (int i=0; i<M*N; i++)
		{
			distance_des_sort[i] = IloNumArray(env, 3);		//To store flow, source, destination
		}
		sort(env, dist, distance_des_sort, 1);
		cout<<"Distance Array SORTED in Ascending Order"<<distance_des_sort<<endl;
		//system("pause");
		/*---------------------------------INITIALISATION VIA HEURISTIC--------------------------------*/

		/* ------------------------------- START HEURISTIC TO GENERATE COLUMNS ------------------------------- */	
		cout<<"Start Heuristic"<<endl;
	  // (i) Declare required matrices:
			// a. combined matrix [S R M N Flow Distance]
			Num2DMatrix comb(env, S*R);
			for(int sr=0; sr<S*R; sr++){
				comb[sr] = IloNumArray(env, 6);
			}

			// b. SI matrix
			Int2DMatrix SI(env, S);
			for(int s=0; s<S; s++){
				SI[s] = IloIntArray(env,2);
			}
			// Initialize each element of SI matrix with -100
			for(int s=0; s<S; s++){
				for(int col=0; col<2; col++){
					SI[s][col]=-100;
				}
			}
			
			// c. RJ matrix
			Int2DMatrix RJ(env, R);
			for(int r=0; r<R; r++){
				RJ[r] = IloIntArray(env,2);
			}
			for(int r=0; r<R; r++){
				for(int col=0; col<2; col++){
					RJ[r][col]=-100;
				}
			}
			cout<<"\nI'm at ij"<<endl;
			// d. IJ matrix						/*---------------------New Addition here------------------------------------*/
			Int2DMatrix IJ(env, M);
			for(int i=0; i<M;i++)
			{
				IJ[i] = IloIntArray(env,2);
			}
			for(int i=0; i<M; i++){
				for(int col=0; col<2; col++){
					IJ[i][col]=-100;
				}
			}
			cout<<"\nij worked"<<endl;
			// e. Initial Column Matrix: Dimension = Col_Matrix[R + S*I +1]by[J+1]
				Num2DMatrix Col_Matrix(env, R + S*M + 1);
				for(int row=0; row<R + S*M + 1; row++){
					Col_Matrix[row] = IloNumArray(env, N + 1);
				}
      
		// (ii) Call required function
				cout<<"\nCalling new function"<<endl;
				gen_cols(comb, S, R, SI, RJ, IJ, flow_des_sort, distance_des_sort, Col_Matrix);
				
		// (iii)Set initial objective function values to the Col_Matrix columns
				IloInt s1,m1,n1,r1;
				for(int j=0;j<Col_Matrix[0].getSize()-1;j++)
				{
					obj_init[j]=0;
					r1 = RJ[j][0];
					n1 = RJ[j][1];						//To find the j of the corresponding unique r value
					/*for(int ij =0;ij<RJ.getSize();ij++)
					{
						if(RJ[ij][0] == r1)				//To find the j of the corresponding unique r value
							n1 = RJ[ij][1];
					}*/
					for(int k=0;k<Col_Matrix.getSize()-1;k++)
					{
						if(k<R)
							continue;
						s1 = (k-R)/S;
						m1 = (k-R)%S;
						//cout<<"S:M:N:R= "<<s1<<":"<<m1<<":"<<n1<<":"<<r1<<endl;
						if(Col_Matrix[k][j]!=0)
						{
							obj_init[j] += flow[s1][r1]*dist[m1][n1];
							//cout<<"Flow: "<<flow[s1][r1]<<" Dist: "<<dist[m1][n1]<<endl;
						}
					}
					Column_index_r.push_back(r1);					//For initial y-columns added
					Column_index_j.push_back(n1);					//For initial y-columns added
					cout<<"obj_init[j]: "<<obj_init[j]<<endl;
				}
				Column_index_r.push_back(R);						//For initial x-column added
				Column_index_j.push_back(N);						//For initial x-column added
				obj_init[Col_Matrix[0].getSize()-1] = 0;			//For initial x-column added
				system("pause");
		// OUTPUTS of gen_col_matrix function:
				cout<<'\n';
				cout<<'\n';
				cout<<"Combined matrix"<<comb<<endl;
				cout<<'\n';
	      		cout<<"SI Matrix"<<SI<<endl;	
				cout<<'\n';
				cout<<"RJ Matrix"<<RJ<<endl;
				cout<<'\n'<<"|********* END: comb matrix update *********|";
				cout<<'\n';
				////system("pause");
				cout<<'\n';
				cout<<"Column Matrix"<<Col_Matrix<<endl;
				cout<<"Sizes: "<<Col_Matrix.getSize()<<":"<<Col_Matrix[0].getSize()<<endl;
				cout<<'\n'<<"|********* END: Column matrix update *********|";
				cout<<'\n'<<"|--------------- END: HEURISTIC TO GENERATE INITIAL COLUMNS ---------------|";
				cout<<'\n';
				cout<<"Heuristic Generated columns"<<endl;
				system("pause");
				
/* ------------------------------- END: HEURISTIC TO GENERATE COLUMNS ------------------------------- */

 /*--------------------------------MASTER PROBLEM: SELECTING AMONG GIVEN COLUMNS---------------------------------*/
	  cout<<"Master Problem"<<endl;
	  
	  IloBoolArray UB_master_constraint(env, R+(S*M)+1);
	  IloBoolArray LB_master_constraint(env, R+(S*M)+1);
	  
	  for (int r=0; r<R; r++)
	  {
		  UB_master_constraint[r] = 1;
		  LB_master_constraint[r] = 1;
	  }

	  for (int s=0; s<S; s++)
	  {
		  for (int i=0; i<M; i++)
		  {
			  UB_master_constraint[R+s*S+i] = 0;
			  LB_master_constraint[R+s*S+i] = -MAX;
		  }
	  }

	  UB_master_constraint[R+(S*M)] = 1;
	  LB_master_constraint[R+(S*M)] = 1;
	  cout<<"UB_master_constraint: "<<UB_master_constraint<<endl;
	  cout<<"LB_master_constraint: "<<LB_master_constraint<<endl;
	  
	  Constraint_master = IloAdd(model_master, IloRangeArray(env, LB_master_constraint, UB_master_constraint));
	  /*Col_select.add(IloNumVar(Col_cost(27.5) + Constraint_master[0](1) + Constraint_master[2](flow[0][0]) + Constraint_master[10](flow[S-1][0])));
	  Col_select.add(IloNumVar(Col_cost(22.5) + Constraint_master[1](1) + Constraint_master[6](flow[1][1]) + Constraint_master[10](flow[S-1][1])));
	  Col_select.add(IloNumVar(Col_cost(0) + Constraint_master[2](-1) + Constraint_master[6](-1) + Constraint_master[10](-1) + Constraint_master[11](1)));
	  */
	 
	  for(int j=0; j<Col_Matrix[0].getSize();j++)
	  {
		IloNumColumn col = Col_cost(obj_init[j]);
		for (int i = 0; i < Constraint_master.getSize(); ++i)
		{
			col += Constraint_master[i](Col_Matrix[i][j]);
		}
		IloNumVar var(col);
		Col_select.add(var);
		//Col_select.add(IloNumVar(Col_cost(27.5) + Constraint_master[0](1) + Constraint_master[2](flow[0][0]) + Constraint_master[10](flow[S-1][0])));
		//Col_select.add(IloNumVar(Col_cost(22.5) + Constraint_master[1](1) + Constraint_master[6](flow[1][1]) + Constraint_master[10](flow[S-1][1])));
		//Col_select.add(IloNumVar(Col_cost(0) + Constraint_master[2](-1) + Constraint_master[6](-1) + Constraint_master[10](-1) + Constraint_master[11](1)));
	  }
	  cout<<"Master problem declaration end"<<endl;
	  system("pause");
/*----------------------------------------END OF MASTER PROBLEM-------------------------------------------------------------------------*/

	  model_sub_y.add(Objective_sub_y);
	  model_sub_x.add(Objective_sub_x);
	  
	  IloBoolVarArray xcolumns(env,S*M);
	  IloBoolVarArray ycolumns(env,S*M);

	  for(int s=0; s<S; s++)
	  {
		  IloExpr LHS_x(env);
		  for(int i=0; i<M; i++)
		  {
			  LHS_x+= xcolumns[s*M + i];
		  }
		  sub_x_constraint.add(IloRange(env, 1, LHS_x, 1));
		  model_sub_x.add(sub_x_constraint[sub_x_constraint.getSize()-1]);
	  }
	  
	  //system("pause");
	  for(int i=0; i<M; i++)
	  {
		  IloExpr LHS_x(env);
		  for(int s=0; s<S; s++)
		  {
			  LHS_x+= xcolumns[s*M + i];
		  }
		  sub_x_constraint.add(IloRange(env, 1, LHS_x, 1));
		  model_sub_x.add(sub_x_constraint[sub_x_constraint.getSize()-1]);
		  LHS_x.end();
	  }
	  
	  //system("pause");
	  model_sub_x.add(sub_x_constraint);
	  	  
	  shadow_price.setSize(R+(S*M)+1);
	  ycolumns_val.setSize(R+(S*M)+1);
	  xcolumns_val.setSize(R+(S*M)+1);
	  
	  for(int r=0; r<R; r++)
	  {
		ycolumns_val[r] = 0;//Later for column corresponding to r, newCol_val[r] set to 1 for y column
		xcolumns_val[r] = 0;
	  }
	  ycolumns_val[ycolumns_val.getSize()-1] = 0;//Last element in y column is always 0
	  xcolumns_val[xcolumns_val.getSize()-1] = 1;//Last element in x column is always 1
	  cout<<"newCol_y_val "<<ycolumns_val<<endl;
	  cout<<"newCol_x_val "<<xcolumns_val<<endl;

/*---------------------------------------------SOME INITIAL DECLARATIONS FOR YCOLUMN - BEGIN-------------------------------------------------*/
	  /*Column_index_r.push_back(0);//Initial Columns Y
	  Column_index_j.push_back(0);//Initial Columns Y
	  //Column_index_plan.push_back(1);//Initial Columns Y
	  Column_index_r.push_back(1);//Initial Columns Y
	  Column_index_j.push_back(1);//Initial Columns Y
	  //Column_index_plan.push_back(1);//Initial Columns Y
	  Column_index_r.push_back(R);//Initial Columns (X)
	  Column_index_j.push_back(N);//Initial Columns (X)
	  //Column_index_plan.push_back(R+N);//Initial Columns (X)
	  */
	  /// COLUMN-GENERATION PROCEDURE ///
	  
	  
	  // Below values are also hard coded; however these should be extracted from the flow matrix. 
	  // Change this once the code is finalized
	  // 1. r=1;
		/*flowretailer[1].numsuppliers = 2;
		// 1.1 s=1;
		flowretailer[1].supplierids[1] = 1;
		flowretailer[1].flow[1] = 1;//0.5;
		
		// 1.2 s=2;
		flowretailer[1].supplierids[2] = 3;
		flowretailer[1].flow[2] = 0.5;//0.25;

		// 2. r=2;
		flowretailer[2].numsuppliers = 2;
		// 2.1 s=1;
		flowretailer[2].supplierids[1] = 2;
		flowretailer[2].flow[1] = 1;//0.5;
		
		// 2.2 s=2; 
		flowretailer[2].supplierids[2] = 3;
		flowretailer[2].flow[2] = 0.5;//0.25;*/

		/**Ex_Illustrative1**/
		// 1. r=1;
		/*flowretailer[2].numsuppliers = 2;
		// 1.1 s=1;
		flowretailer[2].supplierids[1] = 1;
		flowretailer[2].flow[1] = 1;
		
		// 1.2 s=2;
		flowretailer[2].supplierids[2] = 3;
		flowretailer[2].flow[2] = 0.5;

		// 2. r=2;
		flowretailer[1].numsuppliers = 2;
		// 2.1 s=1;
		flowretailer[1].supplierids[1] = 2;
		flowretailer[1].flow[1] = 1;
		
		// 2.2 s=2; 
		flowretailer[1].supplierids[2] = 3;
		flowretailer[1].flow[2] = 0.5;

		*/
	  //env.out()<<flow<<endl;
		for(int j=1;j<=R;j++)
		{
			flowretailer[j].numsuppliers=0;
			for(int k=1;k<=S;k++)
			{
				//cout<<flowretailer[j].numsuppliers<<":"<<flow[k-1][j-1]<<endl;
				if(flow[k-1][j-1]==0)
					continue;
				flowretailer[j].numsuppliers++;
				flowretailer[j].supplierids[flowretailer[j].numsuppliers]=k;
				flowretailer[j].flow[flowretailer[j].numsuppliers]=flow[k-1][j-1];
			}
		}
		//return 0;
	  basissize = R + S*M + 1;
	  //Initial solution will have all but 1 y solutions and 1 x solution in the same order 
	  for(int j=0; j<Col_Matrix[0].getSize()-1;j++)
	  {
		  isY[j] = 1;			//Since these are y-columns
	  }
	  isY[Col_Matrix[0].getSize()-1] = 0;
	  /*isY[0]=1;
	  isY[1]=1;
	  isY[2]=0;*/

	  cout<<"\n\nMade initial declaration for y columns\n"<<endl;
	  //system("pause");
	  
/*---------------------------------------------SOME INITIAL DECLARATIONS FOR YCOLUMN - END -------------------------------------*/
		//cplex_master.setParam(IloCplex::EpOpt, 1e-2);
		//cplex_master.setParam(IloCplex::EpAGap, 1e-2);
		//cplex_master.setParam(IloCplex::EpGap, 1e-2);
		//cplex_master.setParam(IloCplex::EpInt, 1e-2);
		//cplex_master.setParam(IloCplex::EpRelax, 1e-2);
		//cout<<cplex_master.getParam(IloCplex::EpOpt)<<endl;
		system("pause");
		cplex_master.extract(model_master);		//Read master problem
		solve(cplex_master);
		cout<<res<<endl;
		if(res==0)				//If relaxed master problem is infeasible then original problem is infeasible
		{
			cout<<"Problem infeasible!!"<<endl;
			return 0;
		}
		system("pause");
		//cplex_master.exportModel("master.lp");
		cplex_master.getValues(vals,Col_select); //Copy values of solution vector
		cout<<"Solution value: "<<cplex_master.getObjValue()<<endl;
		cout<<"Solution vector: "<<vals<<endl;
		//report1 (cplex_master, Col_select, Constraint_master, Column_index_r, Column_index_j, S, R, M, N);
		cout<<"\n\nSolved relaxed lp once completely!\n"<<endl;
		cout<<"Start Branching"<<endl;
		cout<<"Iter:S,M,N,R"<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
		/*------------------------------------------------------Branch and bound-----------------------------------------------------*/
		system("pause");
			
		if(vals.areElementsInteger())				//Are the values in the solution vector integral
		{
			env.out()<<"Relaxed LP gives integral solution"<<endl;
			env.out()<<"\n\nFinal incumbent solution: "<<cplex_master.getObjValue()<<endl;
			env.out()<<"Final incumbent solution vector: "<<vals<<endl;
			//t = clock() - t;
			//env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
			
			//env.out() << "Solution vector = " << vals << endl;
			//getch();
			return 0;
		}
		
		eps = cplex_master.getParam(IloCplex::EpInt);	//Take epsilon value from cplex
		
		//return 0;
		/*Queue of problems*/
		queue<IloModel> mod;		//Create a queue of models
		mod.push(model_master);		//Push the master problem
		if(min_obj==1)				//Is it minimization or maximization problem
			incumbent = MAX;
		else incumbent = -MAX;
		
		cout<<"Queue created!"<<endl;
		//system("pause");
		
		branch(mod, incumbent, cplex_master, Col_select);	//Start branching
		
		env.out()<<"---------------------------------------------"<<endl;
		env.out()<<"\n\nFinal incumbent solution: "<<incumbent<<endl;
		env.out()<<"Final incumbent solution vector: "<<intSol<<endl;
		//t = clock() - t;
		//env.out()<<"Time of execution in seconds: "<<((float)t)/CLOCKS_PER_SEC;
		system("pause");


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


void solve(IloCplex &cplex_master)
{
	res=0;				//If a feasible solution exists then 1 else 0
	RC_FLAG=1;			//Until atleast 1 x-column or y-column is added RC_FLAG=1 else 0
	cout<<"Node number: "<<nodenum<<endl;	//Solve called at which node
	cplex_master.exportModel("model.lp");
	//cout<<cplex_master.solve()<<endl;
	//system("pause");
	while (RC_FLAG)
	{
		  RC_FLAG = 0;
		  /*--------------------------------------OPTIMIZE OVER CURRENT COLUMNS -------------------------------------------------*/
		  
		  if(!cplex_master.solve())						//Check infeasibility of solution by cplex
		  {
			  res=0;
			  cout<<"Infeasible"<<endl;
			  return;
		  }
		  else res=1;
		  //cout<<"check2: "<<res<<endl;
		  //report1 (cplex_master, Col_select, Constraint_master, Column_index_r, Column_index_j, S, R, M, N);
		  
		  
		  //system("pause");
		  /*---------------------------------------------FIND AND ADD y/x COLUMNS------------------------------------------------*/
		  for(int n=0; n<shadow_price.getSize(); n++)
		  {
			  shadow_price[n] = cplex_master.getDual(Constraint_master[n]);
			  if (IloAbs(shadow_price[n]) < RC_EPS)
			  {
				  shadow_price[n] = 0.0;
			  }
		  }
		  
		  /// FIND AND ADD A NEW Y COLUMN WITH NEGATIVE REDUCED COST ///
		  // 1. Read shadow prices
		 cout<<"Reading shadow prices!"<<endl;
		 for(int n=0; n<shadow_price.getSize(); n++)
		 {
			  node[nodenum].dual[n+1]= shadow_price[n];
			  printf("\t[%d, %f]\n",n,node[nodenum].dual[n+1]);
		 }
		 cout<<'\n';
		 //cout<<"Solve()::R,N: "<<R<<","<<N<<endl;
		 //system("pause");

		 for (int r=0; r<R; r++)
		 {
			  for(int j=0; j<N; j++)
			  {//start y

				  for(int i=0; i<100; i++)
				  {	  // Before begining every r-j pair, makes every entry of reduced cost = 0;
					  reduced_cost_y_column[i] = 0;
				  }

				  numycolumns = 0; // For now hard code this value; otherwise this should come from the code depending upon how many columns 
					               // have previously been added, and increment this value as and when a ycolumn is added

/*------------------------------------------------------EXPERIMENTS HERE-----------------------------------------------------------------------------------------------*/
				  //if(nodenum==0)
				  //{		//For testing x-gen.....remove later
					cout<<"Solving y column generation problem for Destination "<<r+1<<", Outbound Door "<<j+1<<endl;
					//cout<<"Iter::r,n,N: "<<r<<","<<j<<","<<N<<endl;
					cout<<S<<":"<<M<<":"<<N<<":"<<R<<":"<<r<<":"<<j<<":";//<<dist<<":"<<reduced_cost_y_column<<endl;
					//system("pause");
					//S=4;
					//M=4;
					//N=3;
					//R=3;
					//system("pause");
					generateycolumn(nodenum, S, M, N, R, r, j, dist, reduced_cost_y_column);
					//cout<<"Iter: "<<S<<":"<<M<<":"<<N<<":"<<R<<":"<<r<<":"<<j<<":";
					cout<<"Iter:S,M,N,R"<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
					//system("pause");
					//cout<<"check3: "<<res<<endl;
					// Reading values of ycolumn[][] into ycolumns_val[] to add to restricted master problem
					int y_node_output = 0; // r-j semi-assignment solution number(p) for (p+1)th best
				  
					for(int r1=0; r1<shadow_price.getSize(); r1++)
					{
						ycolumns_val[r1] = ycolumns_semiassignment[r1+1][y_node_output+1];
					}
				  
					IloNum current_sub_red_cost = -reduced_cost_y_column[y_node_output+1];;
					cout<<'\n'<<'\n'<<"y_sub_problem_objective value: "<<current_sub_red_cost<<'\n';
					//system("pause");
					if (current_sub_red_cost< -RC_EPS)
					{
						IloNum Col_cost_val = ycolumns_semiassignment[0][y_node_output+1];
		  				cout<<"new_y_Col_val "<<ycolumns_val<<'\n';
						cout<<"Cost of the new Column: "<<Col_cost_val<<endl;
						RC_FLAG = 1;
	 					cout<<'\n'<<"RC_FLAG is 1 because there is a new y-column for r="<<r+1<<" and j="<<j+1<<'\n'<<'\n';
						//system("pause");
						Col_select.add(IloNumVar(Col_cost(Col_cost_val) + Constraint_master(ycolumns_val)));
						//cout<<"check1"<<endl;
						isY[Col_Counter]=1;			//Keeps track of the column added is y column or xcolumn
						Col_Counter++;				//Number of columns added
						
						Column_index_r.push_back(yheader[numycolumns].r-1);
						Column_index_j.push_back(yheader[numycolumns].j-1);
					  
					}//if (current_sub_red_cost< -RC_EPS)
					else
					{
						cout<<'\n'<<"There is no new y-column for r="<<r+1<<" and j="<<j+1<<'\n'<<'\n';
						//system("pause");
					}
					  
				  //}
				  ycolumns_val[r] = 0;//Change it back to 0 by default. Will be changed to 1 for r inside the loop
			  }//for(int j=0; j<Num_odoors; j++) - endy
		 }//for (int r=0; r<Num_destination; r++)
		  cout<<"Solving X column generation problem "<<endl;
		  //cout<<"S,M,N,R: "<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
		  //system("pause");
		  /// FIND AND ADD A NEW X COLUMN WITH NEGATIVE REDUCED COST ///
		  // 1. Declare 2D Array, 2D Result Matrix and reduced cost place holder that will go as an input to hungarian
		  //double node_number = 0.0;				//kth best
		  IloInt node_number = 0;					//kth best
		  if(nodenum!=0)
		  {
			  cout<<"k-best: "<<node_number<<endl;
			  //system("pause");
		  }
		  if(node[nodenum].is_xincluded==0)			//Is x column already included? If not then
		  {
			  Num2DMatrix kbestsols(env, NUM_SOLUTIONS*S*S);
			  for (int s=0; s<NUM_SOLUTIONS*S*S; s++)
			  {
				  kbestsols[s] = IloNumArray(env, 2);
			  }
			  //system("pause");	  
			  IloNumArray koptvals(env, NUM_SOLUTIONS);
		  
			  // 2. Call Miller function

			  Millers_Code(S, R, kbestsols, koptvals, shadow_price);
			  cout<<'\n'<<"Objective function from assignment problem: "<<koptvals[node_number]<<'\n';
			  koptvals[node_number] += -shadow_price[shadow_price.getSize()-1];//shadow price of last constraint (w)
			  cout<<'\n'<<"Final reduced cost of x sub-problem: "<<koptvals[node_number]<<'\n';
			  k_best:								//label k_best used by goto below		  
			  //system("pause");
			  /* Consider root node only*/
			  for(int s=0; s<S; s++)
			  {
				  for(int i=0; i<M; i++)
				  {
					  if(kbestsols[node_number*S*S+s*S+i][0] == node_number && kbestsols[node_number*S*S+s*S+i][1] == 1.0)
					  {
						  xcolumns_val[R + s*M+i] = -1.0;
					  }	
					  else
					  {
						  xcolumns_val[R + s*M+i] = 0.0;
					  }
				  }
			  }
		  
			  xcolumns_val[xcolumns_val.getSize()-1] = 1;//Last element in y column is always 0 for y column (1 for x column)
			  cout<<"newCol_x_val "<<xcolumns_val<<endl;
			  //cout<<"newCol_x_val "<<newCol_x_val<<endl;
			  //report2 (cplex_sub_x, xcolumns_val,  Objective_sub_x);
			  ////system("pause");
			  cout<<"xcol_check: "<<xcol_check<<endl;
			  /*for(int i=0;i<xcolumns_val.getSize();i++)
				cout<<xcolumns[xcol_check][i]<<" ";
			  */cout<<endl;
			  //system("pause");
			  //Col_select.remove(IloNumVar(Col_cost(0) + Constraint_master(xcolumns[xcol_check])));
			  /*int is_excl=1;
			  for(int i=0;i<node[nodenum].numxexclude;i++)
			  {
				  if(xcolumns_val[i]!=xcolumns[xcol_check][i])
					  is_excl=0;
			  }*/ //Check excluded here
			  //if(is_excl==1 && node[nodenum].numxexclude!=0)
			  if((koptvals[node_number]< -RC_EPS))//&&node_number==0.0 || node_number!=0.0)
			  {
				  int dummy=1;			//ADD X-COLUMN IF DUMMY IS 1
				  /*if(nodenum!=0)
				  {
					  //cout<<"Check"<<endl;
					  for(int i=0;i<xcolumns_val.getSize();i++)
					  {
						cout<<xcolumns_val[i]<<":"<<node[nodenum].xexclude[i]<<endl;
					  }
					  //system("pause");
				  }*/
				  for(int j=0;j<node[nodenum].numxexclude;j++)	//Check feasible xcolumn with the x columns in the exclude list
				  {
					dummy=0;
					if(xcol_check==-1)			//No xcolumn added
					{
						dummy=1;
						break;
					}
					for(int i=0;i<xcolumns_val.getSize();i++)
					{
						//cout<<xcolumns_val[i]<<": "<<xcol[node[nodenum].xexclude[j]][i]<<endl;
						if(xcolumns_val[i]!=xcol[node[nodenum].xexclude[j]][i])		//If a value of xcolumn is not same as the corresponding one of exclude column, check next
						{	
							dummy=1;
							break;
						}
					}
					//system("pause");
					if(dummy==0)
						break;
				  }
				  //If dummy=1 add new xcolumn
				  if(dummy==1)
				  {
					RC_FLAG = 1;
					cout<<"Storing xcolumns!"<<endl;
					//cout<<xcol_ctr<<":"<<xcolumns_val.getSize()<<endl;
					for(int j=0;j<xcolumns_val.getSize();j++)
					{
						xcol[xcol_ctr][j]=xcolumns_val[j];		//Keeps track of number of xcolumns added					
					}
					xcol_ctr++;			//Number of xcolumns added
					Col_select.add(IloNumVar(Col_cost(0) + Constraint_master(xcolumns_val)));//Obj fn coefficient of x column always 0
					isY[Col_Counter]=0;		//Column added is xcolumn
					Col_Counter++;
					//cout<<"Column Index : "<<j<<"\t"<<newCol_x_val<<endl;
					Column_index_r.push_back(R);//Column_index_r for y column varies from 0 t0 Num_destination-1
					Column_index_j.push_back(N);//Column_index_j for y column varies from 0 t0 Num_odoors-1
					//Column_index_plan.push_back(R+N);//I don't know
					cout<<'\n'<<'\n'<<"RC_FLAG is 1 because there is a new x-column"<<'\n'<<'\n';
					//if(nodenum!=0)
						//return;
				   }//system("pause");
				  else			//If xcolumn in exclude list then check next best xcolumn
				  {
					  node_number+=1;
					  cout<<koptvals[node_number]<<endl;
					  //Next kth best
					  if(node_number<=3)			//Since there are 6 possible xcolumns(hard-coded)
					  {
						  cout<<"Generating"<<node_number<<" best"<<endl;
						  //system("pause");
						  goto k_best;			//Goto the label k_best
					  }
					  else
					  {
						  res=0;
						  break;											  
					  }
				  }
			  }
			  /*else
			  {
				  if(nodenum!=0)		//Temporary
				  {
						  RC_FLAG=1;
						  node_number+=1.0;
				  }
				  cout<<'\n'<<'\n'<<"There is no new x-column"<<'\n'<<'\n';
				  
				  //system("pause");
			  }*/
		  }//if(node[nodenum].xincluded==0)
		  //cout<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
		  //system("pause");
	  }//while (RC_FLAG)
		if(cplex_master.getObjValue()>=feasibility)		//Check infeasibility if the solution is a large number meaning a column could not be made 0
		{
			res=0;
			//break;
			cout<<"Big value!"<<endl;
			system("pause");
		}
}

void branch(queue<IloModel> mod, IloNum &incumbent, IloCplex& cplex_master, IloNumVarArray var)
{		
	/*Load current model(branch) with non-Integer solutions*/
	/*string left = "left", right="right";
	ctr++;
	left=left+to_string(static_cast<long long>(ctr))+".lp";
	right=right+to_string(static_cast<long long>(ctr))+".lp";*/
	//cout<<left<<endl
	//cplex_master.setParam(IloCplex::EpOpt, 1e-2);
	//cplex_master.setParam(IloCplex::EpAGap, 1e-2);
	//cplex_master.setParam(IloCplex::EpGap, 1e-2);
	//cplex_master.setParam(IloCplex::EpInt, 1e-2);
	//cplex_master.setParam(IloCplex::EpRelax, 1e-2);
	/*Copy includes and excludes of the previous node*/
	int xexclude[100],yinclude[100],yexclude[100],numxexclude,numxinclude,numyexclude,numyinclude;
	IloNum xinclude[100];
	numxexclude=node[nodenum].numxexclude;
	numxinclude=node[nodenum].numxinclude;
	numyexclude=node[nodenum].numyexclude;
	numyinclude=node[nodenum].numyinclude;

	for(int j=0;j<node[nodenum].numxexclude;j++)
		xexclude[j]=node[nodenum].xexclude[j];
	for(int j=0;j<node[nodenum].numxinclude;j++)
		xinclude[j]=node[nodenum].xinclude[j];
	for(int j=0;j<node[nodenum].numyexclude;j++)
		yexclude[j]=node[nodenum].yexclude[j];
	for(int j=0;j<node[nodenum].numyinclude;j++)
		yinclude[j]=node[nodenum].yinclude[j];

	
	iter=0;
	if(mod.size()==0)
		return;
	IloModel model(mod.front());
	mod.pop();
	//cout<<"Inside branch"<<endl;
	//cout<<"Incumbent: "<<incumbent<<endl;
	//cplex_master.exportModel("master.lp");
	//system("pause");

	cplex_master.getValues(vals, Col_select);
	
	/*Locate position of non-integrality*/
	int i;
	IloNum value;
	int xbranch=0;
	//Branch on x-column first
	for( i=0;i<vals.getSize();i++)
	{
		if(isY[i]==0)						//Keeps track of which xcolumn in an array of xcolumn vectors
			xcol_check++;
		if(abs(vals[i]-floor(vals[i]))>=eps && isY[i]==0)		
		{
			//xcol_check++;
			xbranch=1;
			break;
		}
	}
	//Branch on y-column if xcolumn not found to branch on
	if(xbranch==0)
	{
		for( i=0;i<vals.getSize();i++)
		{
			if(isY[i]==0)
				xcol_check++;
			if(abs(vals[i]-floor(vals[i]))>=eps)// && isY[i]==0)		//Temporary...remove 2nd condition of always being x to break later
				break;
		}
	}
	int xcol_curr_node = xcol_check; 
	value=vals[i];	
	frac = i;
	IloNum objVal;
	//cout<<"Branching on col: "<<frac<<" with value: "<<value<<endl;
	//system("pause");

	/*Branch Left->lower bound*/
	//env.out()<<"\nleft BRANCHED PROBLEM"<<endl;
	//env.out()<<"-----------------------------------------------------------"<<endl;

	IloModel lb(model);	//LOWER BOUND PROBLEM
	nodenum++;			//NEXT NODE
	//cout<<"Node: "<<nodenum<<endl;
	//VALUES COPIED FROM PARENT NODE
	node[nodenum].numxexclude=numxexclude;
	node[nodenum].numxinclude=numxinclude;
	node[nodenum].numyexclude=numyexclude;
	node[nodenum].numyinclude=numyinclude;
	//cout<<"hi"<<endl;
	for(int j=0;j<node[nodenum].numxexclude;j++)
		node[nodenum].xexclude[j]=xexclude[j];
	for(int j=0;j<node[nodenum].numxinclude;j++)
		node[nodenum].xinclude[j]=xinclude[j];
	for(int j=0;j<node[nodenum].numyexclude;j++)
		node[nodenum].yexclude[j]=yexclude[j];
	for(int j=0;j<node[nodenum].numyinclude;j++)
		node[nodenum].yinclude[j]=yinclude[j];
	
	IloNum coeff=0;			//Original Coefficient of the objective function variable that is changed to M
	/*IMPLEMENTING M COEFFICIENT IN OBJECTIVE FUNCTION*/
	IloObjective obj = cplex_master.getObjective();
	IloExpr objexpr = (IloExpr) obj.getExpr();
	//objexpr.ge;
	IloExpr::LinearIterator it(objexpr.getLinearIterator());
	int j=0;
	//cout<<"hi"<<endl;
	for(;it.ok();j++)
	{
		if(j==frac)
		{
			if(isY[j]==1)
				coeff = it.getCoef();
			objexpr.setLinearCoef(Col_select[j],LN);
			break;
		}
		if(isY[j]==1)
			++it;
	}
	for(;j<Col_select.getSize();j++)
	{
		//cout<<"bye"<<endl;
		if(j==frac)
			objexpr.setLinearCoef(Col_select[j],LN);
	}
	cplex_master.getObjective().setExpr(objexpr);
	
	int r1=Column_index_r[frac];		//r of fractional column
	int j1=Column_index_j[frac];		//j of fractional column
	//CHECK IF COLUMN BEING BRANCHED ON IS X-COLUMN OR Y-COLUMN
	//cout<<"Ycolumn? "<<isY[frac]<<endl;
	if(isY[frac]==1)
	{				
		node[nodenum].yexclude[node[nodenum].numyexclude]=frac;
		node[nodenum].numyexclude+=1;		
	}
	else
	{
		cout<<xcol_check<<endl;
		node[nodenum].xexclude[node[nodenum].numxexclude]=xcol_check;
		node[nodenum].numxexclude++;
	}
	solve(cplex_master);

	cout<<"Excluded x List:-"<<endl;
	int c=-1;
	for(int j=0;j<node[nodenum].numxexclude;j++)
	{
		cout<<"Index: "<<node[nodenum].xexclude[j]<<" ";
		c=j;
	}
	if(c!=-1)
	{
		cout<<"Last excluded Column:"<<endl;
		for(int j=0;j<xcolumns_val.getSize();j++)
		{
			cout<<xcol[node[nodenum].xexclude[c]][j]<<": ";					
		}
	}
	cout<<endl;
	cout<<"Included x List:-"<<endl;
	if(node[nodenum].is_xincluded==1)
	{
		for(int j=0;j<R+S*M+1;j++)
		{
			cout<<node[nodenum].xinclude[j]<<" ";
		}
	}
	cout<<endl;
	cout<<"Excluded y List:-"<<endl;
	for(int j=0;j<node[nodenum].numyexclude;j++)
	{
		cout<<node[nodenum].yexclude[j]<<" ";
	}
	cout<<endl;
	cout<<"Included y List:-"<<endl;
	for(int j=0;j<node[nodenum].numyinclude;j++)
	{
		cout<<node[nodenum].yinclude[j]<<" ";
	}
	cout<<endl;
	//system("pause");
	//return;
	
	//Debug
	/*cout<<"Y-column? : "<<isY[i]<<endl;
	cout<<"Non-integral value: "<<value<<endl;
	cout<<"Index: "<<i<<endl;
	cout<<"r: "<<Column_index_r[i]<<endl;
	cout<<"j: "<<Column_index_j[i]<<endl;	
	cout<<"xcol_ctr: "<<xcol_ctr<<endl;
	for(int i=0;i<xcol_ctr;i++)
	{
		for(int j=0;j<R+(S*M)+1;j++)
		{
			cout<<xcolumns[i][j]<<" ";
		}
		cout<<endl;
	}*/
	//cout<<"planid: "<<Column_index_plan[i]<<endl;
	
	cplex_master.exportModel("lb.lp");
	//system("pause");
	cout<<"---------------------SUMMARY OF NODE "<<nodenum<<"(Left)---------------------"<<endl;
	cout<<"Incumbent solution: "<<incumbent<<endl;
	cout<<"Branching on col: "<<frac<<" with value: "<<value<<endl;
	if(res)	//IF FEASIBLE SOLUTION FOUND
	{
		mod.push(lb);		//PUSH THE MODEL INTO QUEUE
		//report1 (cplex_master, Col_select, Constraint_master, Column_index_r, Column_index_j, S, R, M, N);
		cplex_master.getValues(vals,var);
		objVal = cplex_master.getObjValue();
		//env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << objVal << endl;
		env.out() << "Solution vector = " << vals << endl;
		cout<<"Y/X column:- "<<endl;
		for(int j=0;j<vals.getSize();j++)
			cout<<isY[j]<<" ";
		cout<<endl;
		//system("pause");
		//return ;	
		if(vals.areElementsInteger())	//ARE ALL VALUES INTEGER
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))	//UPDATE INCUMBENT
			{
				incumbent = objVal;
				cplex_master.getValues(intSol,var);
				//env.out()<<vals;
			}
		}
		else if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))	//IS IT REQUIRED TO BRANCH
		{
			cout<<"Further branching"<<endl;
			system("pause");
			//return;
			branch(mod,incumbent,cplex_master,var);	//FURTHER BRANCHING
		}
	}	
	else env.out()<<"Lower Bound Infeasible!"<<endl;
	//return;
	//cout<<"Result of solution: "<<res<<endl;
	//cout<<"Node: "<<nodenum<<endl;
	system("pause");
	/*Branch right->upper bound*/
	//env.out()<<"\nright BRANCHED PROBLEM"<<endl;
	//env.out()<<"-----------------------------------------------------------"<<endl;

	IloModel ub(model);	//UPPER BOUND PROBLEM
	nodenum++;		//NEXT NODE
	//cout<<"Node: "<<nodenum<<endl;
	
	xcol_check = xcol_curr_node;		//Since xcol_check is a global variable and it may hold value from different node.
	//VALUES COPIED FROM PARENT NODE
	node[nodenum].numxexclude=numxexclude;
	node[nodenum].numxinclude=numxinclude;
	node[nodenum].numyexclude=numyexclude;
	node[nodenum].numyinclude=numyinclude;

	for(int j=0;j<node[nodenum].numxexclude;j++)
		node[nodenum].xexclude[j]=xexclude[j];
	for(int j=0;j<node[nodenum].numxinclude;j++)
		node[nodenum].xinclude[j]=xinclude[j];
	for(int j=0;j<node[nodenum].numyexclude;j++)
		node[nodenum].yexclude[j]=yexclude[j];
	for(int j=0;j<node[nodenum].numyinclude;j++)
		node[nodenum].yinclude[j]=yinclude[j];

	//CHANGE OBJECTIVE COEFFICIENT FOR UPPER BOUND USING M
	obj = cplex_master.getObjective();
	objexpr = (IloExpr) obj.getExpr();
	objexpr.setLinearCoef(Col_select[frac],coeff);

	for(int j=0;j<Col_select.getSize();j++)
	{
		if((j!=frac && isY[frac]==0 && isY[j]==0) || (j!=frac && isY[frac]==1 && isY[j]==1 && r1==Column_index_r[j])) //All other x cols and all other y cols with same 'r'
			objexpr.setLinearCoef(Col_select[j],LN);
	}
	cplex_master.getObjective().setExpr(objexpr);
	//cplex_master.exportModel("rit.lp");
	//cout<<"Is Y col? "<<isY[frac]<<endl;
	if(isY[frac]==1)
	{
		node[nodenum].yinclude[node[nodenum].numyinclude]=frac;
		node[nodenum].numyinclude+=1;
	}
	else
	{
		for(int j=0;j<R+S*M+1;j++)		
		{
			node[nodenum].xinclude[j]=xcol[xcol_check][j];
		}
		node[nodenum].numxinclude++;
		
		node[nodenum].is_xincluded=1;
		iter=1;
	}
	solve(cplex_master);
	cout<<"Excluded x List:-"<<endl;
	c=-1;
	for(int j=0;j<node[nodenum].numxexclude;j++)
	{
		cout<<node[nodenum].xexclude[j]<<" ";
		c=j;
	}
	cout<<"Last excluded Column:"<<endl;
	if(c!=-1)
	{
		for(int j=0;j<xcolumns_val.getSize();j++)
		{
			cout<<xcol[node[nodenum].xexclude[c]][j]<<": ";					
		}
	}
	cout<<endl;
	cout<<"Included x List:-"<<endl;
	if(node[nodenum].is_xincluded==1)
	{
		for(int j=0;j<R+S*M+1;j++)
		{
			cout<<node[nodenum].xinclude[j]<<" ";
		}
	}
	cout<<endl;
	cout<<"Excluded y List:-"<<endl;
	for(int j=0;j<node[nodenum].numyexclude;j++)
	{
		cout<<node[nodenum].yexclude[j]<<" ";
	}
	cout<<endl;
	cout<<"Included y List:-"<<endl;
	for(int j=0;j<node[nodenum].numyinclude;j++)
	{
		cout<<node[nodenum].yinclude[j]<<" ";
	}
	cout<<endl;
	//system("pause");
	
	//solve(cplex_master);

	cplex_master.exportModel("ub.lp");
	//system("pause");
	cout<<"---------------------SUMMARY OF NODE "<<nodenum<<"(Right)---------------------"<<endl;
	cout<<"Incumbent solution: "<<incumbent<<endl;
	cout<<"Branching on col: "<<frac<<" with value: "<<value<<endl;
	if(res)		//IS FEASIBLE SOLUTION FOUND
	{
		mod.push(ub);
		cplex_master.getValues(vals,var);
		objVal = cplex_master.getObjValue();;
		//env.out() << "Solution status = " << objVal << endl;
		env.out() << "Solution value  = " << cplex_master.getObjValue() << endl;
		env.out() << "Solution vector = " << vals << endl;
		cout<<"Y/X column:- "<<endl;
		for(int j=0;j<vals.getSize();j++)
			cout<<isY[j]<<" ";
		cout<<endl;
		//system("pause");
		if(vals.areElementsInteger())
		{
			if((objVal<incumbent && min_obj ==1) || (objVal>incumbent && min_obj ==0))
			{
				incumbent = objVal;
				cplex_master.getValues(intSol,var);
				//env.out()<<vals;
			}
		}
		else if((objVal<=incumbent && min_obj ==1) || (objVal>=incumbent && min_obj ==0))
		{
			cout<<"Further branching"<<endl;
			system("pause");
			branch(mod,incumbent,cplex_master,var);
		}
	}
	else env.out()<<"Upper Bound Infeasible!"<<endl;
	//cout<<"Result of solution: "<<res<<endl;
	//cout<<"Node: "<<nodenum<<endl;
	system("pause");
	//ub.remove(con2);
	//env.out()<<"Incum: "<<incumbent<<endl;
	return;
}


/*--------------------------Functions for CrossDock column generation problem---------------------------------*/

static void readData (const char* filename, Num2DMatrix& flow, Num2DMatrix& dist)
{
   ifstream in(filename);
   if (in) 
   {
      in >> flow;
      in >> dist;
   }
   else {
      cerr << "No such file: " << filename << endl;
      throw(1);
   }
}

static void report1 (IloCplex& cplex_master, IloNumVarArray Col_select, IloRangeArray Constraint_master, vector <int> Column_index_r, 
	                 vector <int> Column_index_j, int S, int R, int M, int N)
{
   cout << endl;
   cout << "MP Obj Value = " <<cplex_master.getObjValue()<<endl;
   cout << endl;
   for (int j = 0; j < Col_select.getSize(); j++) 
   {
      cout <<"Col" << j << " = " << cplex_master.getValue(Col_select[j]) << endl;
	  if(cplex_master.getValue(Col_select[j]) > 0)
	  {
		  if(Column_index_r[j] < R)//Print only y Columns for destination to Outdoor Assignment
		  {
			  cout<<"Column for Destination "<<Column_index_r[j]+1<<"\t Outbound Door "<<Column_index_j[j]+1<<endl;
			/*  for(int s=0; s<Num_source; s++)
			  {
				  for(int i=0; i<Num_idoors; i++)
				  {

				  }
			  }*/
		  }
	  }
   }

   for (int i = 0; i < Constraint_master.getSize(); i++) 
   {
      cout << "  Shadow Price" << i << " = " << cplex_master.getDual(Constraint_master[i]) << endl;
   }
   cout << endl;
}

static void report2 (IloAlgorithm& cplex_sub, IloNumArray newCol_val, IloObjective Objective_sub)
{
   cout << endl;
   cout << "Reduced cost is " << cplex_sub.getValue(Objective_sub) << endl;
   cout << endl;
   if (cplex_sub.getValue(Objective_sub) <= -RC_EPS) 
   {
	   cout<<"New Column: "<<newCol_val<<endl;
   }
}

static void report3 (IloCplex& cplex_master, IloNumVarArray Col_select)
{
   cout << endl;
   cout << "Best integer solution has Obj = " 
        << cplex_master.getObjValue() <<endl;
   cout << endl;
}

static void sort(IloEnv env, Num2DMatrix Unsorted_matrix, Num2DMatrix& Sorted_array, IloBool order)//Bubble Sort
{
	//Num2DMatrix Sorted_array(env);
	//Unsorted_matrix order: mXn, Sorted array arder: (m*n)X3
	for(int s=0; s<Unsorted_matrix.getSize(); s++)
	{
		for(int r=0; r<Unsorted_matrix[0].getSize(); r++)
		{
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][0] = Unsorted_matrix[s][r];
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][1] = s;
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][2] = r;
		}
	}
	cout<<"Unsorted_matrix Array UNSORTED "<<Sorted_array<<endl;
	IloNumArray swap(env, 3);
	for(int s=0; s<Sorted_array.getSize()-1; s++)
	{
		for(int d=0; d<Sorted_array.getSize()-s-1; d++)
		{
			if (order == 0)
			{
				if (Sorted_array[d][0] < Sorted_array[d+1][0]) /* For descending order use < */
				{
					for(int i=0; i<Sorted_array[0].getSize(); i++)
					{
						swap[i] = Sorted_array[d][i];
						Sorted_array[d][i]   = Sorted_array[d+1][i];
						Sorted_array[d+1][i] = swap[i];
					}
				}
			}//if (order == 0)
			else //(order == 1)
			{
				if (Sorted_array[d][0] > Sorted_array[d+1][0]) /* For ascending order use > */
				{
					for(int i=0; i<Sorted_array[0].getSize(); i++)
					{
						swap[i] = Sorted_array[d][i];
						Sorted_array[d][i]   = Sorted_array[d+1][i];
						Sorted_array[d+1][i] = swap[i];
					}
				}
			}//if (order == 1)
		}//for(int d=0; d<Sorted_array.getSize()-s-1; d++)
	}//for(int s=0; s<Sorted_array.getSize()-1; s++)
}

	//cout<<"Flow Array SORTED "<<flow_des_sort<<endl;
	
	//for (c = 0 ; c < ( n - 1 ); c++)
	//{
	//	for (d = 0 ; d < n - c - 1; d++)
	//	{
	//		if (array[d] > array[d+1]) /* For decreasing order use < */
	//		{
	//			swap       = array[d];
	//			array[d]   = array[d+1];
	//			array[d+1] = swap;
	//		}
	//	}
	//}

	//printf("Sorted list in ascending order:\n");

	//for ( c = 0 ; c < n ; c++ )
	//	printf("%d\n", array[c]);
//}

/* Example Input file:
110
{20, 45, 50, 55, 75}
{48, 35, 24, 10, 8}
*/
static void gen_cols(Num2DMatrix comb, int S, int R, Int2DMatrix SI, Int2DMatrix RJ, Int2DMatrix IJ,
	Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix)
{
	/*First Entry according to heuristic we choose the s,r,m,n corresponding to maximum flow and minimum distance*/
	comb[0][0] = flow_des_sort[0][1]; //Assign s
	comb[0][1] = flow_des_sort[0][2]; //Assign r
	comb[0][2] = distance_des_sort[0][1]; //Assign m
	comb[0][3] = distance_des_sort[0][2]; //Assign n
	comb[0][4] = flow_des_sort[0][0]; //Assign flow
	comb[0][5] = distance_des_sort[0][0]; //Assign distance
	/*Update s,i and r,j  and i,j in these matrices*/
	SI[0][0] = (IloInt)flow_des_sort[0][1]; //Assign s
	SI[0][1] = (IloInt)distance_des_sort[0][1]; //Assign m
	RJ[0][0] = (IloInt)flow_des_sort[0][2]; //Assign r
	RJ[0][1] = (IloInt)distance_des_sort[0][2]; //Assign n
	IJ[0][0] = (IloInt)distance_des_sort[0][1]; //Assign m
	IJ[0][1] = (IloInt)distance_des_sort[0][2]; //Assign n
	int num_rj=1, num_si=1, num_ij=1, num_comb=1;
	//cout<<"\ncomb[0]"<<endl;
	cout<<comb[0]<<endl;
	/*----------------------------COMB MATRIX UPDATE---------------------------------*/
	for(int i = 1;i < flow_des_sort.getSize(); i++)
	{
		comb[i][0] = flow_des_sort[i][1]; //s
		comb[i][1] = flow_des_sort[i][2]; //r
		comb[i][4] = flow_des_sort[i][0]; //flow
		
		/*-----------------j assignment-----------------*/
		int j_assigned=0;
		for(int rj=0;rj<num_rj;rj++)
		{
			if(comb[i][1]==RJ[rj][0])	//Is that r assigned to a j
			{
				j_assigned=1;
				//cout<<"1>Selected j: "<<RJ[rj][1]<<endl;
				comb[i][3] = (IloNum)RJ[rj][1]; //Select that j
			}
			if(j_assigned==1)
				break;
		}
		/*-----------------i assignment-----------------*/
		int i_assigned=0;
		for(int si=0;si<num_si;si++)
		{
			if(comb[i][0]==SI[si][0])	//Is that s assigned to a i
			{
				i_assigned=1;
				//cout<<"1>Selected i: "<<SI[si][1]<<endl;
				comb[i][2] = (IloNum)SI[si][1]; //Select that i
			}
			if(i_assigned==1)
				break;
		}
		/*Consider 4 cases of i,j assignment*/

		if(j_assigned==0 && i_assigned==0)		//Case 1
		{
			for(int dis=0;dis<distance_des_sort.getSize();dis++)
			{
				int unassigned_ij=1;
				for(int rj=0;rj<num_rj;rj++)
				{
					if(distance_des_sort[dis][2]==RJ[rj][1])
					{
						unassigned_ij=0;
						break;
					}
				}
				if(unassigned_ij==1)
				{
					for(int si=0;si<num_si;si++)
					{
						if(distance_des_sort[dis][1]==SI[si][1])
						{
							unassigned_ij=0;
							break;
						}
					}
				}
				if(unassigned_ij==1)
				{
					//Assign the new i,j pair
					comb[i][2] = distance_des_sort[dis][1]; 
					comb[i][3] = distance_des_sort[dis][2];	
					//cout<<"2>Selected i,j: "<<comb[i][2]<<":"<<comb[i][3]<<endl;
					break;
				}
			}
			RJ[num_rj][0] = (IloInt)comb[i][1];
			RJ[num_rj][1] = (IloInt)comb[i][3];
			SI[num_si][0] = (IloInt)comb[i][0];
			SI[num_si][1] = (IloInt)comb[i][2];
			//IJ[num_ij][0] = comb[i][2];
			//IJ[num_ij][1] = comb[i][3];
			num_rj++;
			num_si++;
			//num_ij++;
		}
		else if(j_assigned==1 && i_assigned==0)		//Case 2
		{
			int ij_found=0;
			for(int dis=0;dis<distance_des_sort.getSize();dis++)
			{
				//cout<<"dis: "<<distance_des_sort[dis]<<endl;
				if(distance_des_sort[dis][2]==comb[i][3])
				{
					ij_found=1;
					//cout<<"for j: "<<comb[i][3]<<endl;
					for(int si=0; si<num_si;si++)
					{
						//cout<<"i: "<<distance_des_sort[dis][1]<<":"<<SI[si][1]<<endl;
						if(distance_des_sort[dis][1]==SI[si][1])
						{							
							ij_found=0;
							break;
						}
					}
				}
				if(ij_found==1)
				{
					comb[i][2] = distance_des_sort[dis][1]; 
					comb[i][3] = distance_des_sort[dis][2];
					//cout<<"Chosen i,j: "<<comb[i][2]<<comb[i][3]<<endl;
					break;
				}
			}
			//RJ[num_rj][0] = comb[i][1];
			//RJ[num_rj][1] = comb[i][3];
			SI[num_si][0] = (IloInt)comb[i][0];
			SI[num_si][1] = (IloInt)comb[i][2];
			//IJ[num_ij][0] = comb[i][2];
			//IJ[num_ij][1] = comb[i][3];
			//num_rj++;
			num_si++;
			//num_ij++;
		}
		else if(j_assigned==0 && i_assigned==1)	//Case 3
		{
			int ij_found=0;
			for(int dis=0;dis<distance_des_sort.getSize();dis++)
			{
				if(distance_des_sort[dis][1]==comb[i][2])
				{
					ij_found=1;
					for(int rj=0; rj<num_rj;rj++)
					{
						if(distance_des_sort[dis][2]==RJ[rj][1])
						{
							ij_found=0;
							break;
						}
					}
				}
				if(ij_found==1)
				{
					comb[i][2] = distance_des_sort[dis][1]; 
					comb[i][3] = distance_des_sort[dis][2];
					//cout<<"Chosen i,j: "<<comb[i][2]<<comb[i][3]<<endl;
					break;
				}
			}
			RJ[num_rj][0] = (IloInt)comb[i][1];
			RJ[num_rj][1] = (IloInt)comb[i][3];
			//SI[num_si][0] = comb[i][0];
			//SI[num_si][1] = comb[i][2];
			//IJ[num_ij][0] = comb[i][2];
			//IJ[num_ij][1] = comb[i][3];
			num_rj++;
			//num_si++;
			//num_ij++;
		}
		//cout<<"\ni,j assgn over!"<<endl;
		
		/*-----------------------------distance assignment------------------------*/
		int flag=0;
		for(int c=0;c<distance_des_sort.getSize();c++)
		{
			for(int d=0;d<3;d++)
			{
				if(distance_des_sort[c][1]==comb[i][2] && distance_des_sort[c][2]==comb[i][3])
				{
					comb[i][5]=distance_des_sort[c][0];
					flag=1;
					break;
				}
			}
			if(flag==1)
				break;
		}
		//cout<<"Distance assign over!"<<endl;
		//cout<<"comb update"<<endl;
		cout<<comb[i]<<endl;
		num_comb++;		
	}
	//cout<<"\nComb"<<endl;
	//cout<<comb<<endl;
	
	/*----------------------------------COLUMN UPDATE------------------------------------*/
	//First R rows
	for(int rj=0;rj<RJ.getSize();rj++)
	{
		for(int r=0;r<R;r++)
		{
			if(r==RJ[rj][0])
				Col_Matrix[r][rj]=1;
			else
				Col_Matrix[r][rj]=0;
		}
	}
	
	//For each s,i corresponding to each r,j
	for(int rj=0;rj<RJ.getSize();rj++)
	{
		for(int c=0;c<comb.getSize();c++)
		{
			if(comb[c][1]==RJ[rj][0] && comb[c][3]==RJ[rj][1])
			{
				Col_Matrix[R+(IloInt)comb[c][0]*M+(IloInt)comb[c][2]][rj] = comb[c][4];
			}
		}
	}
	
	//For x column
	for(int r=0;r<R;r++)
		Col_Matrix[r][RJ.getSize()]=0;

	for(int row=R;row<Col_Matrix.getSize()-1;row++)
	{
		for(int col=0;col<Col_Matrix[0].getSize()-1;col++)
		{
			if(Col_Matrix[row][col]!=0)
			{
				Col_Matrix[row][Col_Matrix[0].getSize()-1]=-1;
				break;
			}
			else
				Col_Matrix[row][Col_Matrix[0].getSize()-1]=0;
		}
	}

	//For last row if x or not
	for(int r=0;r<R;r++)
		Col_Matrix[Col_Matrix.getSize()-1][r]=0;
	
	Col_Matrix[Col_Matrix.getSize()-1][R]=1;

	return;
}
/*
static void gen_col_matrix(Num2DMatrix comb, int S, int R, Int2DMatrix SI, Int2DMatrix RJ, Int2DMatrix IJ, 
	Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix)//Generate Initial Columns Matrix
{
		// (ii) Fill S, R & Flow values from flow_des_sort into comb matrix
			for(int row=0; row<S*R; row++)
			{
				for(int col=0; col<3; col++)
				{
					if(col==0)
					{
						comb[row][col+4] = flow_des_sort[row][col];
					}
					else
					{
						comb[row][col-1] = flow_des_sort[row][col];
					}
				}
			}

		// (iii) Fill remaining first row of comb matrix with i, j and distance from distance_des_sort matrix
			comb[0][2] = distance_des_sort[0][1];
			comb[0][3] = distance_des_sort[0][2];
			comb[0][5] = distance_des_sort[0][0];
			
		// (iv) Update SI and RJ matrix with first row of comb matrix
			SI[0][0] = comb[0][0];
			SI[0][1] = comb[0][2];

			RJ[0][0] = comb[0][1];
			RJ[0][1] = comb[0][3];
			
			cout<<"Debug"<<endl;
			cout<<comb<<endl;
			cout<<SI<<endl;
			cout<<RJ<<endl;
			system("pause");

		// (v) Start looping to complete other records of comb matrix
			char S_Flag;
			int s_temp;
			int i_temp;
			char R_Flag;
			int r_temp;
			int j_temp;
			int j_temp1;
			char j_temp1_flag;
			int i_temp1;
			char i_temp1_flag;

			for(int row=1; row<S*R; row++)
			{ // Begin Loop: 1

				// (A) Find whether S from new row of comb matrix already exists in SI matrix
				S_Flag = 'N';
				s_temp = 1000;
				i_temp = 1000;
				for(int s=0; s<S; s++)
				{
					if(comb[row][0] == SI[s][0])
					{
						S_Flag = 'Y';
						s_temp = comb[row][0];
						i_temp = SI[s][1];
					}
				}

				// (B) Find whether R from new row of comb matrix already exists in RJ matrix
				R_Flag = 'N';
				r_temp = 1000;
				j_temp = 1000;
				for(int r=0; r<R; r++){
					if(comb[row][1] == RJ[r][0]){
						R_Flag = 'Y';
						r_temp = comb[row][1];
						j_temp = RJ[r][1];
					}
				}

				// (C) For S_Flag = 'Y' and R_FLag = 'N' find j_temp, while i_temp if already
				//     obtained from S_Flag block (A)
				if(S_Flag == 'Y' && R_Flag == 'N'){
					j_temp1 = 1000;
					for(int ij=0; ij<S*R; ij++){
						if(i_temp == distance_des_sort[ij][1]){
							j_temp1 = distance_des_sort[ij][2];
							j_temp1_flag = 'N';
							for(int r=0; r<R; r++){
								if(j_temp1 == RJ[r][1]){
									j_temp1_flag = 'P';
								}
							}
							if(j_temp1_flag == 'N'){
								j_temp = j_temp1;
								goto update_comb;
							}
						}
					}
				}

				// (D) For S_Flag = 'N' and R_FLag = 'Y' find i_temp, while j_temp if already
				//     obtained from R_Flag block (B)
				if(S_Flag == 'N' && R_Flag == 'Y'){
					i_temp1 = 1000;
					for(int ij=0; ij<S*R; ij++){
						if(j_temp == distance_des_sort[ij][2]){
							i_temp1 = distance_des_sort[ij][1];
							i_temp1_flag = 'N';
							for(int s=0; s<S; s++){
								if(i_temp1 == SI[s][1]){
									i_temp1_flag = 'P';
								}
							}
							if(i_temp1_flag == 'N'){
								i_temp = i_temp1;
								goto update_comb;
							}
						}
					}
				}
				
				// (E) For S_Flag = 'Y' and R_FLag = 'Y' we already have S & R's respective
				//     input door as i_temp and output door as j_temp
				if(S_Flag == 'Y' && R_Flag == 'Y'){
					goto update_comb;
				}

				// (F) For S_Flag = 'N' and R_FLag = 'N': means when both source and destination
				//     from a new row of comb matrix are-not assigned to any doors I & J, respectively
				if(S_Flag == 'N' && R_Flag == 'N'){
					for(int ij=0; ij<S*R; ij++){
						i_temp1_flag = 'N';
						j_temp1_flag = 'N';
						for(int s=0; s<S; s++){
							if(distance_des_sort[ij][1] == SI[s][1]){
								i_temp1_flag = 'P';
							}
						}
						for(int r=0; r<R; r++){
							if(distance_des_sort[ij][2] == RJ[r][1]){
								j_temp1_flag = 'P';
							}
						}
						if(i_temp1_flag == 'N' && j_temp1_flag == 'N'){
							i_temp = distance_des_sort[ij][1];
							j_temp = distance_des_sort[ij][2];
							goto update_comb;
						}
					}
				}

				// (G) Update comb matrix with I, J and distance
				update_comb:
				// (G.1) Update I & J 
				comb[row][2] = i_temp;
				comb[row][3] = j_temp;
				// (G.2) Update distance
				for(int ij=0; ij<S*R; ij++){
					if(comb[row][2] == distance_des_sort[ij][1] && comb[row][3] == distance_des_sort[ij][2]){
						comb[row][5] = distance_des_sort[ij][0];
					}
				}

				// (H) Update RJ Matrix only when S_Flag = 'Y'
				if(S_Flag == 'Y' && R_Flag == 'N'){
					for(int r=0; r<R; r++){
						if(RJ[r][0] == -100){
							RJ[r][0] = comb[row][1];
							RJ[r][1] = comb[row][3];
							break;
						}
					}
				}
				
				// (I) Update SI Matrix only when R_Flag = 'Y'
				if(S_Flag == 'N' && R_Flag == 'Y'){
					for(int s=0; s<S; s++){
						if(SI[s][0] == -100){
							SI[s][0] = comb[row][0];
							SI[s][1] = comb[row][2];
							break;
						}
					}
				}
				
				// (J) Update both SI & RJ Matrix 
				if(S_Flag == 'N' && R_Flag == 'N'){
					for(int s=0; s<S; s++){
						if(SI[s][0] == -100){
							SI[s][0] = comb[row][0];
							SI[s][1] = comb[row][2];
							break;
						}
					}
					for(int r=0; r<R; r++){
						if(RJ[r][0] == -100){
							RJ[r][0] = comb[row][1];
							RJ[r][1] = comb[row][3];
							break;
						}
					}
				}
							

			} // End Loop: 1

			// (K) Update column matrix
				// (K.1) Update R-J in Col_matrix
				int j_temp2;
				for(int r=0; r<R; r++){
					j_temp2 = 1000;
					for(int r1=0; r1<R; r1++){
						if(r == RJ[r1][0]){
							j_temp2 = RJ[r1][1];
						}
					}
					for(int j=0; j<R; j++){
						if(j == j_temp2){
							Col_Matrix[r][j] = 1;
						}
					}
				}
				
				// (K.2) Update S-I-J in Col_matrix
				int si;
				for(int row=0; row<S*R; row++){
					si = 0;
					for(int s=0; s<S; s++){
						for(int i=0; i<S; i++){
							for(int j=0; j<R; j++){
								if( comb[row][0] == s && 
									comb[row][2] == i &&
									comb[row][3] == j){
										Col_Matrix[si+R][j] = comb[row][4];
								}
							}
							si = si + 1;
						}
					}
				}

				// (K.3) Update S-I-J in Col_matrix
				float sum_temp;
				for(int row = R; row < R + S*S; row++){
					sum_temp = 0;
					for(int j=0; j<R; j++){
						sum_temp += Col_Matrix[row][j];
					}
					if(sum_temp>0){
						Col_Matrix[row][R] = -1;
					}
				}
				Col_Matrix[R + S*S][R] = 1;
}
*/
void generateycolumn(int nodenum, int S, int M, int N, int R, int r, int j, Num2DMatrix dist, IloNumArray reduced_cost_y_column) 
{
	cout<<"\nInside y column generation function!"<<endl;
	int m, n, i;
	int nr,nc,k;
	int include[50]={0}; //rth entry will contain 1 or 0 to indicate include y (y=1) or not.
	int exclude[50][50]={0}; //(r,j) entry will contain global column number if to be exclude or 0.
	//node[0].numyinclude=0; //for root node;		
	//cout<<"YO: "<<node[nodenum].numyinclude<<endl;
	int found=0;   //indicates if a profitable column is found
	r++;	
	n=j+1;
	//n=j;
	//cout<<"gen::r,j: "<<r<<","<<n<<endl;
	//cout<<"Include"<<endl;	
	//include[r]=0;
	for (i=0; i < node[nodenum].numyinclude; i++)
	{
		//cout<<node[nodenum].yinclude[i]<<endl;
		//cout<<"r: "<<yheader[node[nodenum].yinclude[i]].r<<endl;
		include[yheader[node[nodenum].yinclude[i]].r+1] = 1; //make rth entry 1  to indicate include y(r)
		printf("\n INCLUDE, r= %d\n", yheader[node[nodenum].yinclude[i]].r);
	}
	
	
	//Debug
	/*cout<<"Exclude"<<endl;
	for(int i=1;i<=node[nodenum].numyexclude;i++)
	{
		cout<<node[nodenum].yexclude[i]<<endl;
		cout<<"r: "<<yheader[node[nodenum].yexclude[i]].r<<endl;
	}*/
	//system("pause");
	//cout<<"r: "<<r<<"Nodenum: "<<nodenum<<"\nnode[nodenum].numyinclude: "<<node[nodenum].numyinclude<<endl;
	//cout<<"yheader[node[nodenum].yinclude[i]].r: "<<yheader[node[nodenum].yinclude[i]].r<<endl;
	//system("pause");
	//yheader[nodenum].r=r;
	//yheader[nodenum].j=j;
	
	//for (r=1; r <= R && found==0; r++){ //for each retailer r
	//for (r=1; r <= R; r++) //for each retailer r  COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line
	//for (r=1; r <= 1; r++) //for each retailer r  COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line
	//cout<<include[r]<<": include[r]"<<endl;
	//system("pause");
	if(include[r]!=1)
	{      // generate column only if y is not included, e.g. y !=1
		//for each r,n combination select profitable plan which is not in the "excluded" list.
		double dualcost[20][100]={0.0}; //store evaluation of each (supplier-inbound dock) in this matrix
		int rank[20][100]={0};       //store rank of evaluation of each (supplier-inbound dock) in this matrix
		double costofcolumn;  //define variable for cost of column or objective value of subproblem
					
		//for (n=1; n <= N && found==0; n++){  //for each outbound dock (j or n) 
		// for (n=1; n <= N ; n++)  //for each outbound dock (j or n)COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line 
		//  for (n=1; n <= 1 ; n++)  //for each outbound dock (j or n)COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line 
		cleanspace();
		yassignment[r][n].nofplans=0;
		std::printf("\n \n FOR Retailer, r=%d,Outdock,j=%d \n",r,n);	
		std::printf("\n [supplier s, indock i, cost C= fsr*V(s,i)-fsr*dij]=\n");
		//system("pause");
		for(i=1; i <= flowretailer[r].numsuppliers; i++)
		{
			int suppierid=flowretailer[r].supplierids[i]; //supplier id
			//clean the space
			for(m=1; m <= M; m++)
			{ //each inbound dock (i or m) 
				int rowfordual=R+(suppierid-1)*S+m; //row reference for dual variable
				set[0].a[i][m]=flowretailer[r].flow[i]*node[nodenum].dual[rowfordual]-flowretailer[r].flow[i]*dist[m-1][n-1];
				//printf("\t[supplierid=%d, m=%d, set[0].a[i][m]=%f]",suppierid,m, set[0].a[i][m]);
				////system("pause");
			}
			std::printf("\n ");
		}

		nr = flowretailer[r].numsuppliers; // For each retailer r nr = number of suppliers for that retailer.
		nc = M;
		int p,q;
		p=getbestsol4semiassnprob(0,nr,nc);
		setsgenerated=1;
		q=0;
		k=1;	
		ksolutions[k].z=set[k].z;
					
		for(i=1;i <= nr; i++)
		{
			ksolutions[k].assignment[i]=set[q].assignment[i];
		}		
					
		double threshold= node[nodenum].dual[r];
		int kgenerated=1;
		int terminate =0;	
		//cout<<"check1: "<<ksolutions[k].z<<":"<<threshold<<":"<<TOL<<endl;

		while(ksolutions[k].z + threshold > TOL && terminate != 1)
		{
			ksolutions[k].z=set[q].z;
			for(i=1;i <= nr; i++)
			{
				ksolutions[k].assignment[i]=set[q].assignment[i];
			}
			//cout<<"check2"<<endl;
			//*******************add the column *****
			costofcolumn = ksolutions[k].z+node[nodenum].dual[r]; //add  Ur to cost of column or objective value
			//cout<<"WOOHOO: "<<costofcolumn<<": "<<TOL<<endl;
			//system("pause");
			if(costofcolumn > TOL)
			{
					//improving column is found
				int overlap=0, checkex;
				cout<<node[nodenum].numyexclude<<endl;
				//system("pause");
				for(checkex=1; checkex <= node[nodenum].numyexclude && overlap==0; checkex++)
				{
					int r1,j1,plan1,t2;
					//yheader[i].r =r;  yheader[i].j=n;yheader[i].planid=planid;
					r1=yheader[node[nodenum].yexclude[checkex]].r;
					j1=yheader[node[nodenum].yexclude[checkex]].j;
					plan1=yheader[node[nodenum].yexclude[checkex]].planid;
					//cout<<"r1: "<<r1<<", j1: "<<j1<<", p1: "<<plan1<<endl;
					//cout<<"r: "<<r<<", j: "<<j<<endl;
					//system("pause");
					if(r1==r && j1==n)
					{
						overlap=1; //assume there is an overlap
						for(t2=1; t2 <= flowretailer[r1].numsuppliers; t2++)
						{
							if(yassignment[r1][j1].assignment[plan1][t2]!= ksolutions[k].assignment[t2])
							{
								overlap=0; //no overlap
							}
						}
						printf("\n[overlap=%d for,r=%d, j=%d,plan=%d\n",overlap,r1,j1,plan1);
					}
				}
				//cout<<"overlap: "<<overlap<<endl;
				if(overlap==0)
				{
					//no overlap found with exclude list
					found=1;
					//generate column and put in RMP
					printf("\n Print Column, basissize=%d \n",basissize);	
					numycolumns++;
					reduced_cost_y_column[numycolumns] = costofcolumn;
								
					//system("pause");
					yassignment[r][n].nofplans++;
					yheader[numycolumns].r=r;
					yheader[numycolumns].j=n;
					yheader[numycolumns].planid=yassignment[r][n].nofplans;
					printf("\n[r=%d,j=%d, plan=%d,costAss=%f,cost=%f ", r,n,yassignment[r][n].nofplans,set[q].z,costofcolumn);
					double objvalue=0.0;
					objvalue=0.0;
					for(int t=1; t <= basissize ; t++)
					{
						ycolumns_semiassignment[t][numycolumns]=0.0;
					}
								
					ycolumns_semiassignment[r][numycolumns]=1.0;
					printf("\nassignment=[");
								
					for(int t1=1; t1 <= flowretailer[r].numsuppliers;t1++)
					{
						yassignment[r][n].assignment[yassignment[r][n].nofplans][t1]=ksolutions[k].assignment[t1];
						printf(" \n ksolutions[k].assignment[t1]%d ", ksolutions[k].assignment[t1]);
						int row=R+(flowretailer[r].supplierids[t1]-1)*S+ksolutions[k].assignment[t1];
						//cout<<'\n'<<'\n'<<"So far so good"<<'\n';
						////system("pause");
						objvalue=objvalue+flowretailer[r].flow[t1]*dist[ksolutions[k].assignment[t1]-1][n-1];
						ycolumns_semiassignment[row][numycolumns]=flowretailer[r].flow[t1];
					}
					//cout<<"Numycolumns: "<<numycolumns<<endl;			
					ycolumns_semiassignment[0][numycolumns]=objvalue;
					//printf("\n[row 0,%f for column=%d]", ycolumns[0][numycolumns],numycolumns);
					printf("\n%f", ycolumns_semiassignment[0][numycolumns]);
					cout<<"\nFinal ycolumn: "<<endl; 
					for(int t=1; t <= basissize ; t++)
					{
						//printf("\n[%d,%f]", t,ycolumns[t][numycolumns]);
						printf("\n%f", ycolumns_semiassignment[t][numycolumns]); // Final ycolumn
					}
					//system("pause");
             				
					//**********************************************************************************************	
				}
							
				q=setpartition(q,nr,nc);
				if (q==0 || setsremaining==0)
				{
					terminate=1;
					kgenerated--;
				}
				else if(set[q].z + threshold > TOL)
				{
					k++;
					kgenerated++;
					ksolutions[k].z=set[q].z;
								
					for(i=1;i <= nr; i++)
					{
						ksolutions[k].assignment[i]=set[q].assignment[i];
					}	
				}
	     		   
			} // improving column is foound when costofcolumn >TOL
		} // while loop ends here

	}
	//cout<<"gen::r,j: "<<r<<","<<n<<endl;
	//cout<<"\nGenY:S,M,N,R"<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
	//cout<<"GenY:S,M,N,R"<<S<<":"<<M<<":"<<N<<":"<<R<<endl;
	//system("pause");
	//cout<<yheader[numycolumns].r<<":"<<	yheader[numycolumns].j<<":"<< yheader[numycolumns].planid<<endl;
	//system("pause");
}

int getbestsol4semiassnprob(int setnum, int nr, int nc) //return 1 if a feasible solution is found else return 0
{
    int i,j; 
	// in each row pick the best column for maximization
    int bestindex[100];
    double bestvalue[100];
	set[setnum].z=0.0;
	for (i=1; i<=nr; i++){	
        bestvalue[i]= -LN;
		bestindex[i]=0;
		for (j=1; j<=nc; j++){
			if (set[setnum].a[i][j] > bestvalue[i]){
					bestvalue[i]=set[setnum].a[i][j];
		 			bestindex[i]=j;	
					//printf("\n insideloop %d,%f",bestindex[i],bestvalue[i]);
			}
		}
		    
		set[setnum].z=set[setnum].z+bestvalue[i];
		set[setnum].assignment[i]=bestindex[i];	
	//printf("[%d,%d,%f] ",i,bestindex[i],bestvalue[i]);
	//	printf("\n");
	}
/*
	for ( i=1; i <= nr; i++) 
	{	
	set[setnum].assignment[i]=bestindex[i];	
	//printf("[%d,%d] ",i,set[setnum].assignment[i]);
	//printf("\n");
	}
*/
	//printf("Total cost of assignment=%f\n ",set[setnum].z);	
	if(set[setnum].z > -LN)
  	return 1;
	else return 0;
}
int setpartition(int setnum,int nr, int nc)
  {
        int i,j,k,isetnumber,setlocation; 
	// create 'nr' number of new sets from set="setnum"
        for (isetnumber=1; isetnumber <= nr; isetnumber++)
           {
		setlocation=0;
		
		for ( i=1; i <= setsgenerated; i++) 
			{	
			if(set[i].status==0) 
                            {
				setlocation=i;
				i=setsgenerated+1;
                                break;
                            }	
			}
		if(setlocation==0)
                 {
		
			setsgenerated++; //increment number of sets generated by 1
			setlocation=setsgenerated;	
			if(setsgenerated >= SETSIZE)
                          {printf("\n not adequate space for set, increase size from %d", setsgenerated);
                           return 0;
                          }
                  }
		set[setlocation].status=1;
		if(isetnumber==nr) set[setnum].status=0; //when last partition of this set is done.

        	//copy parent set data to new set at "setlocation"
		set[setlocation].z=set[setnum].z;
	 	for ( i=1; i <= nr; i++) 
			{	
			set[setlocation].assignment[i]=set[setnum].assignment[i];	
			//printf("[%d,%d] ",i,set[setlocation].assignment[i]);
			//printf("\n");
			}
		for ( i=1; i <= nr; i++) 
			{	
 			for (j=1; j <= nc; j++)  
				{
 				set[setlocation].a[i][j]=set[setnum].a[i][j];	
 				//printf("[%d,%d,%f] ",i,bestindex[i],bestvalue[i]);
				//printf("\n");
				}
			}
			//impose "include" restrictions
                  if(isetnumber >1)
                   {
			for ( k=1; k < isetnumber; k++) 
			{	
			int include=set[setnum].assignment[k];	
 			for (j=1; j <= nc; j++)  
				{
                                 if(j!=include) set[setlocation].a[k][j]=-LN;
  				//printf("[%d,%d,%f] ",i,bestindex[i],bestvalue[i]);
				//printf("\n");
				}
			}
                  }
			//impose "exclude" restrictions for isetnumber
                       {
			int exclude=set[setnum].assignment[isetnumber];
			set[setlocation].a[isetnumber][exclude]=-LN;
                       }
			//print new set
			// printf("\n Set Number=%d location %d\n",isetnumber,setlocation);
		  for ( i=1; i <= nr; i++) 
			{	
			for (j=1; j <= nc; j++)  
				{
           			//printf("[%d,%d,%f] ",i,j,set[setlocation].a[i][j]);
				//printf("[%.3f ]",set[setlocation].a[i][j]);
				}
			//printf("\n");
			}
		 int p;
       		p=getbestsol4semiassnprob(setlocation,nr,nc);
		if(set[setlocation].z < -LN) set[setlocation].status=0;
		//printf("I have added [Cost of assignment=%.3f]\n ",set[setlocation].z);
                      
	    }
	    	
		//printf("\n Generate kth best solution, = %d",k);
		double bestvalue= -LN;
		int bestindex=0;
		setsremaining=0;
         	 for (i=1; i<=setsgenerated; i++)
            		 { 
                	  if(set[i].status==1)
                	 	 {
					//printf("\n[Active set no=%d, Z-Value=%f\n ",i,set[i].z);
			   		if(set[i].z > bestvalue)
                          	 	{
					bestvalue=set[i].z;
                               		bestindex=i;
                          	 	}
			  	setsremaining=setsremaining+1;
				}
	
       			}
		//printf("\n[Zbest set no=%d, Zbest Value=%f\n Assignment",bestindex,set[bestindex].z);
		for(i=1;i <=nr;i++)
                  {
		 	//printf("%d,",set[bestindex].assignment[i]);
                  }
		// printf("]\n ");
   return bestindex;
}

void cleanspace()
{
	int i,j, k;
	for (i=0; i < SETSIZE; i++)
	{
		set[i].z=0;
		set[i].status=0;
		for (j=0; j < 50; j ++)
        {
			set[i].assignment[j]=0;
			for(k=0; k<100 ; k++)
			{
				set[i].a[k][j]=0.0;
			}
		}
	}
}
