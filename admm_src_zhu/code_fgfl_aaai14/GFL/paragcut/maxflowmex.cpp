#include "mex.h"

#include "stdafx.h"
#include "Parser.h"
#include "Solver.h"
#include "Utils.h"
using namespace paraF;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Graph graph;
    InputParser parser;
    Problem problem;
    Solver solver;
    MinCut mc;

    if( !freopen("input.txt","r",stdin) )
        mexErrMsgTxt("failed to open input file.");
    
	parser.Parse(graph,problem);
	mexPrintf("c Nodes: %d\n", graph.nodes.size());
	mexPrintf("c Arcs: %d\n", graph.arcs.size());

	graph.InitPreflowPush(0);
	graph.FindMaxPreflow();
	graph.PrintStats("st");
	graph.ConvertPreflowToFlow();
    graph.FindMinCut(mc);
    
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* flow = mxGetPr(plhs[0]);
    *flow = graph.GetSinkNode()->excess;
    
    plhs[1] = mxCreateNumericMatrix(graph.nodes.size(),1,mxINT32_CLASS,mxREAL);
    int* labels = (int*)mxGetData(plhs[1]);
    for(int i=0; i<graph.nodes.size(); i++)
        labels[i] = mc[i];
}