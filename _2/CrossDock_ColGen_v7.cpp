
// --------------------------------------------------------------------------
// Initial columns are generated using a heuristic
// x-columns sub-problem is solved using Miller's code
/* For the number of assignment solutions required from the Miller's Code, change the value of constant #define NUM_SOLUTIONS in both
   CrossDock_ColGen_v6.cpp and apqueue_test.cpp. */

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
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>//for vectors //ref: http://www.cprogramming.com/tutorial/stl/vector.html
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
//#include<CrossDock_CG_Header.h>

//#include <ilcplex/ilocplex.h>

//ILOSTLBEGIN

#define RC_EPS 1.0e-6//if shadow_price or Reduced Cost < RC_EPS, treat it as 0

#define INF (0x7FFFFFFF)
#define verbose (1)
#define NUM_SOLUTIONS 100
#define TOL 0.00001 /*define tolerance level to check improving obj value*/
#define LN 99999999.0 /* Large number */
#define SETSIZE 50
double ycolumns_semiassignment[30][30];
int basissize, numycolumns, numnodegen, numnodeactive;
double upperbound, lowerbound;

struct YHEADER
{
int planid;
int r;
int j;
} yheader[50];

struct YASSIGNMENT
{
int nofplans; 
int assignment[20][10]; //assignments for each plan
} yassignment[50][50]; 

struct FLOWDATA
{
int numsuppliers, supplierids[20];
double flow[20];
} flowretailer[50];

struct NODE_DATA  //template for data to be stored in each B & P node
{
    int status; //status=1 for active, 0 for pruned
    double lowerbound;
    double upperbound;
    double dual[100];
	int xcolumnindices [100],ycolumnindices[100]; //indices of variables in basis
	int numxinclude, xinclude[2],numxexclude,xexclude[100]; //indices of variables to be excluded
	int numyinclude,yinclude[100],numyexclude,yexclude[100]; //indices of variables to be included
} node[1000];

struct SET
{
int status;
double a[100][50],z;
int assignment[50];
} set[SETSIZE];

int getbestsol4semiassnprob(int setnum, int nr, int nc); // solve semiassignmentproblem
int setpartition(int setnum,int nr, int nc); //partition from node, return index of best set
void cleanspace();
int setsgenerated,setsremaining;

struct YSOLUTIONS
{
double z;
int assignment[50];
} ksolutions[100];

//typedef IloArray<IloNumArray> IloNumArray2;

//static void readData (const char* filename, IloNum& rollWidth, IloNumArray& size, IloNumArray& amount);
//static void report1 (IloCplex& cutSolver, IloNumVarArray Cut, IloRangeArray Fill);
//static void report2 (IloAlgorithm& cplex_sub, IloNumArray Allocation_val, IloObjective Objective_sub);
//static void report3 (IloCplex& cutSolver, IloNumVarArray Cut);

static void readData (const char* filename, Num2DMatrix& flow, Num2DMatrix& distance);
static void sort(IloEnv env, Num2DMatrix Unsorted_matrix, Num2DMatrix& Sorted_array, IloBool order);//Bubble Sort; 
// order = 1 for ascending; order = 0 for descending
static void report1 (IloCplex& cplex_master, IloNumVarArray Col_select, IloRangeArray Constraint_master, vector <int> Column_index_r, vector <int> Column_index_j, IloInt Num_source, IloInt Num_destination, IloInt Num_idoors, IloInt Num_odoors);
static void report2 (IloAlgorithm& cplex_sub, IloNumArray newCol_val, IloObjective Objective_sub);
static void report3 (IloCplex& cplex_master, IloNumVarArray Col_select);
static void gen_col_matrix(Num2DMatrix comb, IloInt Num_source, IloInt Num_destination, Num2DMatrix SI, Num2DMatrix RJ, Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix);
void Millers_Code(int S, int R, Num2DMatrix kbestsols, IloNumArray koptvals, IloNumArray shadow_price);
/// MAIN PROGRAM ///

int main(int argc, char **argv)
{
   IloEnv env;
   try 
   {
	  Num2DMatrix flow(env);
      Num2DMatrix dist(env);

      if ( argc > 1 )
         readData(argv[1], flow, dist);
      else
         readData("Data_CrossDock_CG_Ex_Illusttrative.txt", flow, dist);
	  
	  cout<<flow<<endl;
	  cout<<dist<<endl;

	  IloInt S = flow.getSize();
	  IloInt R = flow[0].getSize();
	  IloInt M = dist.getSize();
	  IloInt N = dist[0].getSize();

	  cout<<"No. of Sources = "<<S<<endl;
	  cout<<"No. of Destinations = "<<R<<endl;
	  cout<<"No. of Inbound Doors = "<<M<<endl;
	  cout<<"No. of Outbound Doors = "<<N<<endl;

	  Num2DMatrix flow_des_sort(env, S*R);//flow matrix as a sorted array in descending order of flows
	  for (IloInt s=0; s<S*R; s++)
	  {
		flow_des_sort[s] = IloNumArray(env, 3);//To store flow, source, destination
	  }
	  sort(env, flow, flow_des_sort, 0);
	  cout<<"Flow Array SORTED in Decending Order"<<flow_des_sort<<endl;

	  Num2DMatrix distance_des_sort(env, M*N);//flow matrix as a sorted array in descending order of flows
	  for (IloInt i=0; i<M*N; i++)
	  {
		distance_des_sort[i] = IloNumArray(env, 3);//To store flow, source, destination
	  }
	  sort(env, dist, distance_des_sort, 1);
	  cout<<"Distance Array SORTED in Ascending Order"<<distance_des_sort<<endl;
/* ------------------------------- START HEURISTIC TO GENERATE COLUMNS ------------------------------- */	
	  // (i) Declare required matrices:
			// a. combined matrix [S R M N Flow Distance]
			Num2DMatrix comb(env, S*R);
			for(int sr=0; sr<S*R; sr++){
				comb[sr] = IloNumArray(env, 6);
			}

			// b. SI matrix
			Num2DMatrix SI(env, S);
			for(int s=0; s<S; s++){
				SI[s] = IloNumArray(env,2);
			}
			// Declare each element of SI matrix with -100
			for(int s=0; s<S; s++){
				for(int col=0; col<2; col++){
					SI[s][col]=-100;
				}
			}
			
			// c. RJ martix
			Num2DMatrix RJ(env, R);
			for(int r=0; r<R; r++){
				RJ[r] = IloNumArray(env,2);
			}
			for(int r=0; r<R; r++){
				for(int col=0; col<2; col++){
					RJ[r][col]=-100;
				}
			}

			// d. Initial Column Matrix: Dimension = Col_Matrix[R + S*I +1]by[J+1]
				Num2DMatrix Col_Matrix(env, R + S*M + 1);
				for(int row=0; row<R + S*M + 1; row++){
					Col_Matrix[row] = IloNumArray(env, N + 1);
				}
      
		// (ii) Call required function
				gen_col_matrix(comb, S, R, SI, RJ, flow_des_sort, distance_des_sort, Col_Matrix);

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
				cout<<'\n'<<"|********* END: Columne matrix update *********|";
				cout<<'\n'<<"|--------------- END: HEURISTIC TO GENERATE INITIAL COLUMNS ---------------|";
				cout<<'\n';
				////system("pause");

/* ------------------------------- END: HEURISTIC TO GENERATE COLUMNS ------------------------------- */

	  /// MASTER PROBLEM: SELECTING AMONG GVEN COLUMNS ///
	  IloModel model_master(env);
	  IloObjective Col_cost = IloAdd(model_master, IloMinimize(env));
	  IloBoolArray UB_master_constraint(env, R+(S*M)+1);
	  IloBoolArray LB_master_constraint(env, R+(S*M)+1);
	  
	  for (IloInt r=0; r<R; r++){
		  UB_master_constraint[r] = 1;
		  LB_master_constraint[r] = 1;
	  }

	  for (IloInt s=0; s<S; s++){
		  for (IloInt i=0; i<M; i++){
			  UB_master_constraint[R+s*S+i] = 0;
			  LB_master_constraint[R+s*S+i] = -IloInfinity;
		  }
	  }

	  UB_master_constraint[R+(S*M)] = 1;
	  LB_master_constraint[R+(S*M)] = 1;
	  cout<<"UB_master_constraint: "<<UB_master_constraint<<endl;
	  cout<<"LB_master_constraint: "<<LB_master_constraint<<endl;
	  
	  IloRangeArray  Constraint_master = IloAdd(model_master, IloRangeArray(env, LB_master_constraint, UB_master_constraint));
	  IloNumVarArray Col_select(env);

	  //Initial Columns
	  /*for (IloInt r=0; r<Num_destination; r++)
	  {
		  Col_select.add(IloNumVar(Col_cost(27.5) + Constraint_master[r](1)));
	  }*/
	  /*Col_select.add(IloNumVar(Col_cost(27.5) + Constraint_master[0](1) + Constraint_master[Num_destination](flow[0][0]) + Constraint_master[Num_destination+Num_source*Num_idoors-1](flow[Num_source-1][0])));
	  Col_select.add(IloNumVar(Col_cost(22.5) + Constraint_master[1](1) + Constraint_master[Num_destination+Num_source+2-1](flow[1][1]) + Constraint_master[Num_destination+Num_source*Num_idoors-1](flow[Num_source-1][1])));
	  Col_select.add(IloNumVar(Col_cost(0) + Constraint_master[2](-1) + Constraint_master[2+3+1](-1) + Constraint_master[2+3+3+2](-1) + Constraint_master[2+3+3+3](1)));
   */

	  Col_select.add(IloNumVar(Col_cost(27.5) + Constraint_master[0](1) + Constraint_master[2](flow[0][0]) + Constraint_master[10](flow[S-1][0])));
	  Col_select.add(IloNumVar(Col_cost(22.5) + Constraint_master[1](1) + Constraint_master[6](flow[1][1]) + Constraint_master[10](flow[S-1][1])));
	  Col_select.add(IloNumVar(Col_cost(0) + Constraint_master[2](-1) + Constraint_master[6](-1) + Constraint_master[10](-1) + Constraint_master[11](1)));
   

	  IloCplex cplex_master(model_master);
	  
	  /// COLUMN-GENERATION PROBLEM- Column_y ///
	  IloModel model_sub_y(env);
	  IloBoolVarArray ycolumns(env, S*M);
	  //IloObjective ReducedCost_y = IloAdd(model_sub_y, IloMinimize(env, 1));
	  IloObjective Objective_sub_y = IloMinimize(env);
	  model_sub_y.add(Objective_sub_y);
	  //Add constraints to the y column subproblem
	  //IloRangeArray sub_y_constraint(env, Num_source);
	  IloRangeArray sub_y_constraint(env);//Constraints added later according to r
	  IloCplex cplex_sub_y(model_sub_y);
	  
	  /// COLUMN-GENERATION PROBLEM- Column_x ///
	  IloModel model_sub_x(env);
	  //IloObjective ReducedCost_x = IloAdd(model_sub_x, IloMinimize(env, 1));
	  IloBoolVarArray xcolumns(env, S*M);
	  IloObjective Objective_sub_x = IloMinimize(env);
	  model_sub_x.add(Objective_sub_x);
	  //Add constraints to the y column subproblem
	  IloRangeArray sub_x_constraint(env);

	  for(IloInt s=0; s<S; s++)
	  {
		  IloExpr LHS_x(env);
		  for(IloInt i=0; i<M; i++)
		  {
			  LHS_x+= xcolumns[s*M + i];
		  }
		  //sub_x_constraint[s] = IloRange(env, 1, LHS_x, 1); // \sum_i b(s, i) =  \forall s
		  sub_x_constraint.add(IloRange(env, 1, LHS_x, 1));
		  model_sub_x.add(sub_x_constraint[sub_x_constraint.getSize()-1]);
		  //IloRange sub_x_constraint(env, 1, LHS_x, 1); // \sum_i b(s, i) =  \forall s
		  //model_sub_x.add(sub_x_constraint);
		  LHS_x.end();
	  }

	  for(IloInt i=0; i<M; i++)
	  {
		  IloExpr LHS_x(env);
		  for(IloInt s=0; s<S; s++)
		  {
			  LHS_x+= xcolumns[s*M + i];
		  }
		  sub_x_constraint.add(IloRange(env, 1, LHS_x, 1));
		  model_sub_x.add(sub_x_constraint[sub_x_constraint.getSize()-1]);
		  //sub_x_constraint[Num_destination + Num_source + i] = IloRange(env, 1, LHS_x, 1); // \sum_s b(s, i) =  \forall i
		  //IloRange sub_x_constraint(env, 1, LHS_x, 1); // \sum_s b(s, i) =  \forall i
		  //model_sub_x.add(sub_x_constraint);
		  LHS_x.end();
	  }
	  //sub_x_constraint[Num_source + Num_idoors] = IloRange(env, 1, newCol_x[Num_source*Num_idoors], 1);//Last element in y column is always set to 1
	  //IloRange sub_x_constraint(env, 1, newCol_x[Num_destination + Num_source*Num_idoors + 0], 1); //Last element in y column is always set to 1
	  model_sub_x.add(sub_x_constraint);
	  IloCplex cplex_sub_x(model_sub_x);

	  IloNumArray shadow_price(env, R+(S*M)+1);//To store shadow prices
	  IloNumArray ycolumns_val(env, R+(S*M)+1);//To store the new column (y or x) obtained
	  IloNumArray xcolumns_val(env, R+(S*M)+1);//To store the new column (y or x) obtained
	  for(IloInt r=0; r<R; r++)
	  {
		ycolumns_val[r] = 0;//Later for column corresponding to r, newCol_val[r] set to 1 for y column
		xcolumns_val[r] = 0;
	  }
	  ycolumns_val[ycolumns_val.getSize()-1] = 0;//Last element in y column is always 0
	  xcolumns_val[xcolumns_val.getSize()-1] = 1;//Last element in x column is always 1
	  cout<<"newCol_y_val "<<ycolumns_val<<endl;
	  cout<<"newCol_x_val "<<xcolumns_val<<endl;
	  vector <int> Column_index_r;//This is to map the retailer (r) associated with the column added; r=/Num_destination for x column
	  vector <int> Column_index_j;//This is to map the outbound door (j) associated with the the column added; j=Num_odoors for x column

	  Column_index_r.push_back(0);//Initial Columns Y
	  Column_index_j.push_back(0);//Initial Columns Y
	  Column_index_r.push_back(1);//Initial Columns Y
	  Column_index_j.push_back(1);//Initial Columns Y
	  Column_index_r.push_back(R);//Initial Columns (X)
	  Column_index_j.push_back(N);//Initial Columns (X)

	  /// COLUMN-GENERATION PROCEDURE ///
	  IloNum RC_FLAG = 1;//Flag indicating negative reduced cost for at least one y or x column
	  IloInt Col_Counter=5;
	  
	  /******* SOME INITIAL DECLARATIONS FOR YCOLUMN - BEGIN *******/
	  // Below values are also hard coded; however these should be extracted from the flow matrix. 
	  // Change this once the code is finalized
	  // 1. r=1;
	  flowretailer[1].numsuppliers = 2;
		// 1.1 s=1;
		flowretailer[1].supplierids[1] = 1;
		flowretailer[1].flow[1] = 1;
		
		// 1.2 s=2;
		flowretailer[1].supplierids[2] = 3;
		flowretailer[1].flow[2] = 0.5;

	  // 2. r=2;
      flowretailer[2].numsuppliers = 2;
		// 2.1 s=1;
	    flowretailer[2].supplierids[1] = 2;
		flowretailer[2].flow[1] = 1;
		
		// 2.2 s=2; 
		flowretailer[2].supplierids[2] = 3;
		flowretailer[2].flow[2] = 0.5;

	  basissize = R + S*M + 1;

	  IloNumArray reduced_cost_y_column(env, 100);
	  /******* SOME INITIAL DECLARATIONS FOR YCOLUMN - END *******/

	  while (RC_FLAG)
	  {
		  RC_FLAG = 0;
		  /// OPTIMIZE OVER CURRENT COLUMNS ///
		  cplex_master.solve();
		  report1 (cplex_master, Col_select, Constraint_master, Column_index_r, Column_index_j, S, R, M, N);
		  /// FIND AND ADD y/x COLUMNS///
		  for(IloInt n=0; n<shadow_price.getSize(); n++){
			  shadow_price[n] = cplex_master.getDual(Constraint_master[n]);
			  if (IloAbs(shadow_price[n]) < RC_EPS){
				  shadow_price[n] = 0.0;
			  }
		  }
		  
		  /// FIND AND ADD A NEW Y COLUMN WITH NEGATIVE REDUCED COST ///
		  // 1. Read shadow prices
		 for(IloInt n=0; n<shadow_price.getSize(); n++){
			  node[0].dual[n+1]= shadow_price[n];
			  printf("\t[%d, %f]\n",n,node[0].dual[n+1]);
		}
		 cout<<'\n';
		 ////system("pause");

		  for (IloInt r=0; r<R; r++){
			  for(IloInt j=0; j<N; j++){//start y

				  for(int i=0; i<100; i++){ // Before begining every r-j pair, makes every entry of reduced cost = 0;
					  reduced_cost_y_column[i] = 0;
				  }

				  numycolumns = 0; // For now hard code this value; otherwise this should come from the code depending upon how many columns 
					               // have previously been added, and increment this value as and when a ycolumn is added

				  //generate y column by solving assignment subproblem 1. This is only a definition.
				  void generateycolumn(int nodenum, int S, int M, int N, int R, int r, int j, Num2DMatrix dist, IloNumArray reduced_cost_y_column);
				  cout<<"Solving y column generation problem for Destination "<<r+1<<", Outbound Door "<<j+1<<endl;
				 // //system("pause");

				  generateycolumn(0, S, M, N, R, r, j, dist, reduced_cost_y_column);

				  // Reading values of ycolumn[][] into ycolumns_val[] to add to restricted master problem
				  int y_node_output = 0; // r-j semi-assignment solution number
				  
				  for(int r=0; r<shadow_price.getSize(); r++){
					  ycolumns_val[r] = ycolumns_semiassignment[r+1][y_node_output+1];
				  }
				  
				  
				  /* I don't think we need this loop, since the ycolumn[][] output from Dr. Adil's code is already multiplied with the flow
				  for(IloInt s=0; s<S; s++){
					  for(IloInt i=0; i<M; i++){
						  if (cplex_sub_y.getValue(ycolumns[s*M+i]) > 0.5){
							  ycolumns_val[R+s*M+i] = 1.0*flow[s][r];
						  }
						  else{
							  ycolumns_val[R+s*M+i] = 0.0*flow[s][r];
						  }
					  }
				  }
				  */
				  //system("pause");
				  //report2(cplex_sub_y, ycolumns_val,  Objective_sub_y);
				  ////system("pause");
				  IloNum current_sub_red_cost = -reduced_cost_y_column[y_node_output+1];;
				  cout<<'\n'<<'\n'<<"y_sub_problem_objective value: "<<current_sub_red_cost<<'\n';
				  //system("pause");
				  if (current_sub_red_cost< -RC_EPS){
					  /*IloNum Col_cost_val = 0; 
					  for(IloInt s=0; s<S; s++){
						  for(IloInt i=0; i<M; i++){
							  Col_cost_val+= ycolumns_val[R + s*M + i]*dist[i][j];//Column elements already multiplied by flows
							  cout<<"Col_cost_val "<<Col_cost_val<<endl;
						  }
					  }*/
					  IloNum Col_cost_val = ycolumns_semiassignment[0][y_node_output+1];
		  			  cout<<"new_y_Col_val "<<ycolumns_val<<'\n';
					  cout<<"Cost of the new Column: "<<Col_cost_val<<endl;
					  RC_FLAG = 1;
	 				  cout<<'\n'<<"RC_FLAG is 1 because there is a new y-column for r="<<r+1<<" and j="<<j+1<<'\n'<<'\n';
					  system("pause");
					  Col_select.add(IloNumVar(Col_cost(Col_cost_val) + Constraint_master(ycolumns_val)));
					  Column_index_r.push_back(r);
					  Column_index_j.push_back(j);
				  }//if (current_sub_red_cost< -RC_EPS)
				  else{
					  cout<<'\n'<<"There is no new y-column for r="<<r+1<<" and j="<<j+1<<'\n'<<'\n';
					  system("pause");
				  }
					  
			  }//for(IloInt j=0; j<Num_odoors; j++) - endy
			  ycolumns_val[r] = 0;//Change it back to 0 by default. Will be changed to 1 for r inside the loop
		  }//for (IloInt r=0; r<Num_destination; r++)
		  
		  cout<<"Solving X column generation problem "<<endl;
		  /// FIND AND ADD A NEW X COLUMN WITH NEGATIVE REDUCED COST ///
		  // 1. Declare 2D Array, 2D Result Matrix and reduced cost place holder that will go as an input to hungarian
		  Num2DMatrix kbestsols(env, NUM_SOLUTIONS*S*S);
		  for (IloInt s=0; s<NUM_SOLUTIONS*S*S; s++){
			  kbestsols[s] = IloNumArray(env, 2);
		  }
		  	  
		  IloNumArray koptvals(env, NUM_SOLUTIONS);
		  
		  // 2. Call Miller function
		  double node_number = 0.0;

		  Millers_Code(S, R, kbestsols, koptvals, shadow_price);
		  cout<<'\n'<<"Objective function from assignment problem: "<<koptvals[node_number]<<'\n';
		  koptvals[node_number] += -shadow_price[shadow_price.getSize()-1];//shadow price of last constraint (w)
		  cout<<'\n'<<"Final reduced cost of x sub-problem: "<<koptvals[node_number]<<'\n';
		  
		  /* Consider root node only*/
		  for(IloInt s=0; s<S; s++){
			  for(IloInt i=0; i<M; i++){
				  if(kbestsols[node_number*S*S+s*S+i][0] == node_number && kbestsols[node_number*S*S+s*S+i][1] == 1.0){
					  xcolumns_val[R + s*M+i] = -1.0;
				  }	
				  else{
					  xcolumns_val[R + s*M+i] = 0.0;
				  }
			  }
		  }
		  
		  xcolumns_val[xcolumns_val.getSize()-1] = 1;//Last element in y column is always 0 for y column (1 for x column)
		  cout<<"newCol_x_val "<<xcolumns_val<<endl;
		  //cout<<"newCol_x_val "<<newCol_x_val<<endl;
		  //report2 (cplex_sub_x, xcolumns_val,  Objective_sub_x);
		  ////system("pause");
		  
		  //system("pause");
		  if (koptvals[0]< -RC_EPS){
			  RC_FLAG = 1;
			  Col_select.add(IloNumVar(Col_cost(0) + Constraint_master(xcolumns_val)));//Obj fn coefficient of x column always 0
			  //cout<<"Column Index : "<<j<<"\t"<<newCol_x_val<<endl;
			  Column_index_r.push_back(R);//Column_index_r for y column varies from 0 t0 Num_destination-1
			  Column_index_j.push_back(N);//Column_index_j for y column varies from 0 t0 Num_odoors-1
			  cout<<'\n'<<'\n'<<"RC_FLAG is 1 because there is a new x-column"<<'\n'<<'\n';
			  system("pause");
		  }
		  else{
			  cout<<'\n'<<'\n'<<"There is no new x-column"<<'\n'<<'\n';
			  system("pause");
		  }
	  }//while (RC_FLAG)
	  
	  
      /// COLUMN-GENERATION PROCEDURE ///
	 
      model_master.add(IloConversion(env, Col_select, ILOBOOL));

      cplex_master.solve();
      report3 (cplex_master, Col_select);
	  
   }
   catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

   env.end();

   return 0;
}


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
	                 vector <int> Column_index_j, IloInt S, IloInt R, IloInt M, IloInt N)
{
   cout << endl;
   cout << "MP Obj Value = " <<cplex_master.getObjValue()<<endl;
   cout << endl;
   for (IloInt j = 0; j < Col_select.getSize(); j++) 
   {
      cout <<"Col" << j << " = " << cplex_master.getValue(Col_select[j]) << endl;
	  if(cplex_master.getValue(Col_select[j]) > 0)
	  {
		  if(Column_index_r[j] < R)//Print only y Columns for destination to Outdoor Assignment
		  {
			  cout<<"Column for Destination "<<Column_index_r[j]+1<<"\t Outbound Door "<<Column_index_j[j]+1<<endl;
			/*  for(IloInt s=0; s<Num_source; s++)
			  {
				  for(IloInt i=0; i<Num_idoors; i++)
				  {

				  }
			  }*/
		  }
	  }
   }

   for (IloInt i = 0; i < Constraint_master.getSize(); i++) 
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
	for(IloInt s=0; s<Unsorted_matrix.getSize(); s++)
	{
		for(IloInt r=0; r<Unsorted_matrix[0].getSize(); r++)
		{
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][0] = Unsorted_matrix[s][r];
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][1] = s;
			Sorted_array[s*Unsorted_matrix[0].getSize()+r][2] = r;
		}
	}
	cout<<"Unsorted_matrix Array UNSORTED "<<Sorted_array<<endl;
	IloNumArray swap(env, 3);
	for(IloInt s=0; s<Sorted_array.getSize()-1; s++)
	{
		for(IloInt d=0; d<Sorted_array.getSize()-s-1; d++)
		{
			if (order == 0)
			{
				if (Sorted_array[d][0] < Sorted_array[d+1][0]) /* For decreasing order use < */
				{
					for(IloInt i=0; i<Sorted_array[0].getSize(); i++)
					{
						swap[i] = Sorted_array[d][i];
						Sorted_array[d][i]   = Sorted_array[d+1][i];
						Sorted_array[d+1][i] = swap[i];
					}
				}
			}//if (order == 0)
			else //(order == 1)
			{
				if (Sorted_array[d][0] > Sorted_array[d+1][0]) /* For decreasing order use < */
				{
					for(IloInt i=0; i<Sorted_array[0].getSize(); i++)
					{
						swap[i] = Sorted_array[d][i];
						Sorted_array[d][i]   = Sorted_array[d+1][i];
						Sorted_array[d+1][i] = swap[i];
					}
				}
			}//if (order == 0)
		}//for(IloInt d=0; d<Sorted_array.getSize()-s-1; d++)
	}//for(IloInt s=0; s<Sorted_array.getSize()-1; s++)
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
static void gen_col_matrix(Num2DMatrix comb, IloInt S, IloInt R, Num2DMatrix SI, Num2DMatrix RJ, 
	Num2DMatrix flow_des_sort, Num2DMatrix distance_des_sort, Num2DMatrix Col_Matrix)//Generate Initial Columns Matrix
{
		// (ii) Fill S, R & Flow values from flow_des_sort into comb matrix
			for(int row=0; row<S*R; row++){
				for(int col=0; col<3; col++){
					if(col==0){
						comb[row][col+4] = flow_des_sort[row][col];
					}
					else{
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

			for(int row=1; row<S*R; row++){ // Begin Loop: 1

				// (A) Find whether S from new row of comb matrix already exists in SI matrix
				S_Flag = 'N';
				s_temp = 1000;
				i_temp = 1000;
				for(int s=0; s<S; s++){
					if(comb[row][0] == SI[s][0]){
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

void generateycolumn(int nodenum, int S, int M, int N, int R, int r, int j, Num2DMatrix dist, IloNumArray reduced_cost_y_column) 
{
	int s, m, n, i;
	int nr,nc,k;
	int include[50]={0}; //rth entry will contain 1 or 0 to indicate include y (y=1) or not.
	int exclude[50][50]={0}; //(r,j) entry will contain global column number if to be exclude or 0.
		//node[0].numyinclude=0; //for root node;		
       for (i=1; i <= node[nodenum].numyinclude; i++){
			include[yheader[node[nodenum].yinclude[i]].r] = 1; //make rth entry 1  to indicate include y(r)
			printf("\n INCLUDE, r= %d\n", yheader[node[nodenum].yinclude[i]].r);
		}
	system("pause");
        int found=0;   //indicates if a profitable column is found
		 r++;	
		 n=j+1;
		 //for (r=1; r <= R && found==0; r++){ //for each retailer r
			//for (r=1; r <= R; r++) //for each retailer r  COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line
			//for (r=1; r <= 1; r++) //for each retailer r  COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line
			if(include[r]!=1){      // generate column only if y is not included, e.g. y !=1
				//for each r,n combination select profitable plan which is not in the "excluded" list.
				double dualcost[20][100]={0.0}; //store evaluation of each (supplier-inbound dock) in this matrix
				int rank[20][100]={0};       //store rank of evaluation of each (supplier-inbound dock) in this matrix
				double costofcolumn;  //define variable for cost of column or objective value of subproblem
					
				//for (n=1; n <= N && found==0; n++){  //for each outbound dock (j or n) 
					// for (n=1; n <= N ; n++)  //for each outbound dock (j or n)COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line 
					//  for (n=1; n <= 1 ; n++)  //for each outbound dock (j or n)COMMENT THIS LINE AFTER TESTING IS COMPLETE and UNCOMMENT previous line 
          	 		cleanspace();
					yassignment[r][n].nofplans=0;
					printf("\n \n FOR Retailer, r=%d,Outdock,j=%d \n",r,n);	
					printf("\n [supplier s, indock i, cost C= fsr*V(s,i)-fsr*dij]=\n");		
					for(i=1; i <= flowretailer[r].numsuppliers; i++){
						int suppierid=flowretailer[r].supplierids[i]; //supplier id
						//clean the space
						for(m=1; m <= M; m++){ //each inbound dock (i or m) 
							int rowfordual=R+(suppierid-1)*S+m; //row reference for dual variable
							set[0].a[i][m]=flowretailer[r].flow[i]*node[nodenum].dual[rowfordual]-flowretailer[r].flow[i]*dist[m-1][n-1];
							//printf("\t[supplierid=%d, m=%d, set[0].a[i][m]=%f]",suppierid,m, set[0].a[i][m]);
							////system("pause");
						}
						printf("\n ");
					}

					nr = flowretailer[r].numsuppliers; // For each retailer r nr = number of suppliers for that retailer.
					nc = M;
					int p,q;
       				p=getbestsol4semiassnprob(0,nr,nc);
					setsgenerated=1;
					q=0;
					k=1;	
					ksolutions[k].z=set[k].z;
					
					for(i=1;i <= nr; i++){
						ksolutions[k].assignment[i]=set[q].assignment[i];
					}		
					
					double threshold= node[nodenum].dual[r];
					int kgenerated=1;
					int terminate =0;	
					
					while(ksolutions[k].z + threshold > TOL && terminate != 1){
						ksolutions[k].z=set[q].z;
						for(i=1;i <= nr; i++){
							ksolutions[k].assignment[i]=set[q].assignment[i];
						}
						
						//*******************add the column *****
						costofcolumn = ksolutions[k].z+node[nodenum].dual[r]; //add  Ur to cost of column or objective value
						
						if(costofcolumn > TOL){  //improving column is found
							int overlap=0, checkex;

							for(checkex=1; checkex <= node[nodenum].numyexclude && overlap==0; checkex++){
								int r1,j1,plan1,t2;
								//yheader[i].r =r;  yheader[i].j=n;yheader[i].planid=planid;
								r1=yheader[node[nodenum].yexclude[checkex]].r;
								j1=yheader[node[nodenum].yexclude[checkex]].j;
								plan1=yheader[node[nodenum].yexclude[checkex]].planid;
								if(r1==r && j1==n){
									overlap=1; //assume there is an overlap
									for(t2=1; t2 <= flowretailer[r1].numsuppliers; t2++){
										if(yassignment[r1][j1].assignment[plan1][t2]!= ksolutions[k].assignment[t2]){
											overlap=0; //no overlap
										}
									}
									printf("\n[overlap=%d for,r=%d, j=%d,plan=%d\n",overlap,r1,j1,plan1);
								}
							}

							if(overlap==0){ //no overlap found with exclude list
								found=1;
								//generate column and put in RMP
								//printf("\n Print Column, basissize=%d \n",basissize);	
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
								for(int t=1; t <= basissize ; t++){
									ycolumns_semiassignment[t][numycolumns]=0.0;
								}
								
								ycolumns_semiassignment[r][numycolumns]=1.0;
								printf("\nassignment=[");
								
								for(int t1=1; t1 <= flowretailer[r].numsuppliers;t1++){
									yassignment[r][n].assignment[yassignment[r][n].nofplans][t1]=ksolutions[k].assignment[t1];
									printf(" %d ", ksolutions[k].assignment[t1]);
									int row=R+(flowretailer[r].supplierids[t1]-1)*S+ksolutions[k].assignment[t1];
									//cout<<'\n'<<'\n'<<"So far so good"<<'\n';
									////system("pause");
									objvalue=objvalue+flowretailer[r].flow[t1]*dist[ksolutions[k].assignment[t1]-1][n-1];
									ycolumns_semiassignment[row][numycolumns]=flowretailer[r].flow[t1];
								}
								
								ycolumns_semiassignment[0][numycolumns]=objvalue;
								//printf("\n[row 0,%f for column=%d]", ycolumns[0][numycolumns],numycolumns);
								printf("\n%f", ycolumns_semiassignment[0][numycolumns]);
								
								for(int t=1; t <= basissize ; t++){
									//printf("\n[%d,%f]", t,ycolumns[t][numycolumns]);
									printf("\n%f", ycolumns_semiassignment[t][numycolumns]); // Final yolumn
								}
									//system("pause");
             				
				 //**********************************************************************************************	
							}
							
							q=setpartition(q,nr,nc);
							if (q==0 || setsremaining==0){
								terminate=1;
								kgenerated--;
							}
							else if(set[q].z + threshold > TOL){
								k++;
								kgenerated++;
								ksolutions[k].z=set[q].z;
								
								for(i=1;i <= nr; i++){
									ksolutions[k].assignment[i]=set[q].assignment[i];
								}	
							}
	     		   
						} // improving column is foound when costofcolumn >TOL
					} // while lop ends here

						//}
				}

	 //}

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
