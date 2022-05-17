#include <ilcplex/ilocplex.h>
#include <stdio.h>
using namespace std;

ILOSTLBEGIN
typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;
const int  N_BINS = 10;
const int M = 1000;
IloInt i;
IloInt j;

int bin_SIZE[N_BINS][3] =
{
    

{
    36 ,36, 36},
{
    59, 39, 20},
{
    54, 40, 21},
{
    58 ,37, 21},
{
    52, 33, 20},
{
    40, 31, 21},
{
    31, 31, 17},
{
    31, 17, 16},
{
    26, 23, 14},
{
    33, 21, 4},

};
int bin[2] = {
     80,58 };

int main(void*) {
    
	

	IloEnv env;
	IloModel model(env) ;
	

	IloNumVarArray x(env, N_BINS);
	IloNumVarArray y(env, N_BINS);
	IloNumVarArray z(env, N_BINS);


	for ( i = 0; i < N_BINS; i++) {
    
		x[i]= IloNumVar(env);
		y[i]=IloNumVar(env);
		z[i]=IloNumVar(env);
	}

	IloNumVarArray xl(env, N_BINS);
	IloNumVarArray zl(env, N_BINS);
	IloNumVarArray zh(env, N_BINS);
	IloNumVarArray yw(env, N_BINS);
	for (i = 0; i < N_BINS; i++) {
    
		xl[i]= IloNumVar(env, 0, 1, IloNumVar::Bool);
	
		zl[i]=IloNumVar(env, 0, 1, IloNumVar::Bool);
	
		yw[i]=IloNumVar(env, 0, 1, IloNumVar::Bool);
	
		zh[i]=IloNumVar(env, 0, 1, IloNumVar::Bool);
	}

	IloArray<IloNumVarArray> a(env, N_BINS);
	IloArray<IloNumVarArray> b(env, N_BINS);
	IloArray<IloNumVarArray> c(env, N_BINS);
	for (i = 0; i < N_BINS; i++)
	{
    
		a[i] = IloNumVarArray(env, N_BINS);
		b[i] = IloNumVarArray(env, N_BINS);
		c[i] = IloNumVarArray(env, N_BINS);
		for (j = 0; j < N_BINS; j++)
		{
			a[i][j] =IloNumVar(env, 0, 1, ILOBOOL);
			b[i][j] = IloNumVar(env, 0, 1, ILOBOOL);
			c[i][j] = IloNumVar(env, 0, 1, ILOBOOL);
		}
	}

	IloIntVar h=IloNumVar(env, 0, 100, IloNumVar::Int);

	for ( i = 0u; i < N_BINS; i++)
	{
    
		for ( j = 0u; j < N_BINS; j++) {
    
			if (i == j) continue;
			model.add(x[i] + bin_SIZE[i][0] * xl[i]
				+ bin_SIZE[i][1] * (zl[i] - yw[i] + zh[i])
				+ bin_SIZE[i][2] * (1 - xl[i] - zl[i] + yw[i] - zh[i]) <=x[j]+M*(1-a[i][j]));
			model.add(y[i] + bin_SIZE[i][1] * yw[i]
				+ bin_SIZE[i][0] * (1 - xl[i] - zl[i])
				+ bin_SIZE[i][2] * (xl[i] + zl[i] - yw[i]) <= y[j]+M * (1 - b[i][j]));
			model.add(z[i] + bin_SIZE[i][2] * zh[i]
				+ bin_SIZE[i][1] * (1 - zl[i] - zh[i])
				+ bin_SIZE[i][0] * zl[i] <= z[j]+M * (1 - c[i][j]));
			model.add(a[i][j] + b[i][j] + c[i][j] +a[j][i] + b[j][i] + c[j][i] >= 1);
		}
	}
	
	
	for ( i = 0u; i < N_BINS; i++)
	{
    
		model.add(x[i] + bin_SIZE[i][0] * xl[i]
			+ bin_SIZE[i][1] * (zl[i] - yw[i] + zh[i])
			+ bin_SIZE[i][2] * (1 - xl[i] - zl[i] + yw[i] - zh[i]) <= bin[0]);
		model.add(y[i] + bin_SIZE[i][1] * yw[i]
			+ bin_SIZE[i][0] * (1 - xl[i] - zl[i])
			+ bin_SIZE[i][2] * (xl[i] + zl[i] - yw[i]) <= bin[1]);
		model.add(z[i] + bin_SIZE[i][2] * zh[i]
			+ bin_SIZE[i][1] * (1 - zl[i] - zh[i])
			+ bin_SIZE[i][0] * zl[i] <= h);
		model.add(xl[i] + zl[i] <= 1);
		model.add(zl[i] + zh[i] <= 1);
		model.add(zl[i] -yw[i] + zh[i] <= 1);
		model.add(zl[i] - yw[i] + zh[i] >= 0);
		model.add(1-xl[i] - zl[i]+yw[i] - zh[i] >= 0);
		model.add(1 - xl[i] - zl[i] + yw[i] - zh[i] <=1);
		model.add(xl[i] + zl[i] - yw[i] <= 1);
		model.add(xl[i] + zl[i] - yw[i] >= 0);

		}

	
	
	IloObjective obj(env, h, IloObjective::Minimize);

	
	model.add(obj);

	IloCplex cplex(model);
	bool solved = false;
	cplex.setParam(cplex.TiLim, 60);
	try {
    
		
		solved = cplex.solve();
	}
	catch (const IloException& e) {
    
		std::cerr << "\n\nCPLEX Raised an exception:\n";
		std::cerr << e << "\n";
		env.end();
		throw;
	}

	if (solved) {
    
		
		std::cout << "\n\nCplex success!\n";
		std::cout << "\tStatus: " << cplex.getStatus() << "\n";
		std::cout << "\tObjective value: " << cplex.getObjValue() << "\n";
		for (i = 0u; i < N_BINS; i++)

			std::cout << "\tBINS value:" << cplex.getValue(x[i]) << " " << cplex.getValue(y[i])<<" "<< cplex.getValue(z[i])<<"\n"
		 << "\tBINS parallel:" << cplex.getValue(xl[i]) << " " << cplex.getValue(zl[i]) << " " << cplex.getValue(zh[i])<<" "<< cplex.getValue(yw[i])<<"\n";

		
	}
	else {
    
		std::cerr << "\n\nCplex error!\n";
		std::cerr << "\tStatus: " << cplex.getStatus() << "\n";
		std::cerr << "\tSolver status: " << cplex.getCplexStatus() << "\n";
	}

	env.end();
	return 0;
}