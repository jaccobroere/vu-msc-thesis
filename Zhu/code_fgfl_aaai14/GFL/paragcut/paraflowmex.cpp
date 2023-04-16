#include "mex.h"

#include "stdafx.h"
#include "Parser.h"
#include "Solver.h"
#include "Utils.h"
#include "Graph.cpp"
#include "Parser.cpp"
#include "Solver.cpp"
#include "Utils.cpp"

using namespace paraF;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Graph graph;
    InputParser parser;
    Problem problem;
    Solver solver;

    if( !freopen("input.txt","r",stdin) )
        mexErrMsgTxt("failed to open input file.");
    
	parser.Parse(graph,problem);
    problem.hasLambdaRange = true;
    problem.minLambda = mxGetScalar(prhs[0]);
    problem.maxLambda = mxGetScalar(prhs[1]);
	mexPrintf("c Nodes: %d\n", graph.nodes.size());
	mexPrintf("c Arcs: %d\n", graph.arcs.size());
    
	solver.Solve(graph, problem);

    plhs[0] = mxCreateDoubleMatrix(graph.lambdas.size(),1,mxREAL);
    double* lambdas = (double*)mxGetData(plhs[0]);
    for(int i=0; i<graph.lambdas.size(); i++)
        lambdas[i] = graph.lambdas[i];

    plhs[1] = mxCreateNumericMatrix(graph.nodes.size(),1,mxINT32_CLASS,mxREAL);
    int* labels = (int*)mxGetData(plhs[1]);
    for(int i=0; i<graph.nodes.size(); i++)
        labels[i] = graph.nodes[i].label;    
}