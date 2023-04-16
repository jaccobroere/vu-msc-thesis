# include <mex.h>
# include <stdio.h>
# include <math.h>
# include "quicksort.c"
# include "paragcut/Graph.h"
# include "paragcut/Graph.cpp"
# include "paragcut/Solver.h"
# include "paragcut/Solver.cpp"
# include "time.h"

// this is Efficient fused lasso signal approximation
//  min_w { 0.5*|w-z|_2^2 + lambda1*|w|_1 + lambda2*\sum{|w_i-w_j|} }
// 2013 12 10 


#define EPS 1.0E-10  /** allowed computational error bound **/
using namespace std;
using namespace paraF;
// all 0_index

// global variables
int nV;
double lambda1;
double lambda2;
double *z;
double * Edge_weight;
double * Edge_in;
double * Edge_out;
int nEdge;

double flag;

Graph graph;

double * t;
double * w;
double fw;
double gamma_xb;   //gamma (using gamma_xb to avoid redefinition)





struct mySet
{
	int ns;
	int * S;
};

double f_c(mySet S)
{
	// e'Le
	double fc = 0;

	int * select_stats = new int[nV];
	for(int i = 0; i < nV; i++)
	{
		select_stats[i] = 0;
	}
	for(int i = 0; i < S.ns; i++)
	{
		int ind = S.S[i];
		select_stats[ind] = 1;
	}
	for (int i = 0; i < nEdge; i++)
	{
		int node1 = (int)Edge_in[i]-1;
		int node2 = (int)Edge_out[i]-1;
		if(select_stats[node1]!=select_stats[node2])
		{
			fc = fc + Edge_weight[i];
		}
	}
	delete [] select_stats; 
	fc  = fc/2;
	return fc;
}

double g(mySet S)
{
	// e'Le+gamma_xb*1+(-1/lambda2)*z
	double fz = 0;
	for (int i = 0; i < S.ns; i++)
	{
		int ind = S.S[i];
		fz = fz + (-1/lambda2)*z[ind];
	}

	double f1 = 0;
	f1 = gamma_xb * S.ns;

	double fc = f_c(S);
	
	return f1 + fz + fc;
}

int calc_gamma_xb()
{
	double maxz = 0;
	for (int i = 0; i < nV; i++)
	{
		double tz = z[i];
		if (tz>maxz)
		{
			maxz = tz;
		}
	}
	double maxdegree = 10; // assume max degree in the graph is 10
	gamma_xb = maxdegree + (1/lambda2) * (maxz);
	if(gamma_xb < 0)
	{
		gamma_xb = 0;
	}
	return 0;
}

int paracut()
{

	Problem problem;
	Solver solver;
    double large = 1;
	double top = 1000;
    if(lambda2 < top)
    {
       large = top/lambda2;
    }

	// set problem info 
	problem.type = Problem::PARAMETRIC_FLOW_PROBLEM;
	problem.hasLambdaRange = true;
	problem.minLambda = -10e16*large;
	problem.maxLambda = gamma_xb*large;

	// build graph
	int n = nV+2; // # of nodes
	int m = nEdge+nV*2; // # of edges
	graph.Reserve(m * 2 + n * 4, n);
	  // add notes
	for (int i = 0; i < n; i++)
		graph.AddNode();
	int sourceNodeIndex = 0;
	int sinkNodeIndex = n-1;
	Node* sourceNode = &graph.nodes[sourceNodeIndex];
	Node* sinkNode = &graph.nodes[sinkNodeIndex];

	int tm = 0;
	paraF::Arc* arc;
	paraF::Arc* revArc;
		// add source and sink cap
	double cap1 = 0;
	double cap2 = 0;
	double slope1 = 0;
	double slope2 = 0;
	for(int i = 1; i < n-1; i++)
	{
		Node* currentNode = &graph.nodes[i];
		cap1 = -(1/lambda2)*z[i-1]+gamma_xb*1;
		cap2 = gamma_xb*1;
		slope1 = 0;
		slope2 = -1;
        cap1 = cap1 * large;
        cap2 = cap2 * large;


		arc = graph.AddArcFromSource(sourceNode, currentNode, cap1, slope1);
		revArc = graph.AddArc(currentNode, sourceNode, 0);
		arc->revArc = revArc;
		revArc->revArc = arc;
		tm++;
		arc = graph.AddArcToSink(currentNode, sinkNode, cap2, slope2);
		revArc = graph.AddArc(sinkNode, currentNode, 0);
		arc->revArc = revArc;
		revArc->revArc = arc;
		tm++;
	}
		 // add edge cap
	Node* Node1 = NULL;
	Node* Node2 = NULL;
	
	for (int t = 0; t < nEdge; t ++)
	{
		int i = (int)Edge_in[t];
		int j = (int)Edge_out[t];
		double w = Edge_weight[t];
		if (fabs(w) > EPS)
		{
            w = w * large;
                    
			Node1 = &graph.nodes[i];
			Node2 = &graph.nodes[j];

			arc = graph.AddArc(Node1,Node2,w);
			revArc = graph.AddArc(Node2, Node1, 0);
			arc->revArc = revArc;
			revArc->revArc = arc;
			tm ++;
		}
	}
	
	graph.AddAuxArcs(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);
	graph.PrepareNodes(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);
	graph.PrepareArcs();

	// solve 
	solver.Solve(graph, problem);

	return 0;
}

int calc_t()
{
	int nG = graph.lambdas.size()+1;
	double * g1_value = new double [nG];
	int * n_value = new int[nG];
	mySet T2;
	T2.S = new int [nV];
	int nt2 = 0;
	for (int s = 0; s < nG; s ++)
	{		
		for (int i = 0; i < nV; i++)
		{
			if (graph.nodes[i+1].label==nG+1-s)
			{
				T2.S[nt2] = i;
				nt2 = nt2 + 1;
			}
		}
		T2.ns = nt2;
 		g1_value[s] = g(T2);
		n_value[s] = nt2;
	}
	for ( int i = 0; i < nV; i++)
	{
		int xlabel = graph.nodes[i+1].label;
		int ind = nG-xlabel;
		t[i] = (g1_value[ind+1]-g1_value[ind])/(n_value[ind+1]-n_value[ind]);
	}
	delete [] T2.S;
	delete [] g1_value;
	delete [] n_value;
	return 0;
}

int solve_mnp()
{
	// solve tv using parametric flow

	// compute gamma_xb  Lemma4.4
	calc_gamma_xb();

	// apply paracut to find all S's
	paracut();

	// calc t from S's
	calc_t();

	// post process t   Lemma4.3
	for (int i = 0; i < nV; i++)
	{
		t[i] = t[i] - gamma_xb;
	}

	return 0;
}

double calcOmega()
{
	mySet ts;
	ts.S = new int[nV];
	double fw1 = 0;
	
	for(int j = 0; j < nV; j++)
	{
		ts.S[j] = j;
		if(lambda1!=0)
		{
			fw1 = fw1 + fabs(w[j]);
		}		  
	}
	
	quicksort3(w,ts.S,0,nV-1);
	double sf_previous = 0.0,sf_new = 0.0;
	double fw2 = 0;
	for (int j = 0; j < nV ; j++){
		ts.ns = j+1;
		int ind = ts.S[j];
		sf_new = f_c(ts);
		fw2 +=  w[ind] * (sf_new - sf_previous);
		sf_previous = sf_new;
	}
	delete [] ts.S;
	return lambda1*fw1 + lambda2*fw2;
}

void calcw()
{
	// calc w from t
	for (int i = 0; i < nV; i++)
	{
		w[i] = -lambda2 * t[i];

		if(fabs(w[i])<EPS)
		{
			w[i] = 0;
		}
		
	}    
}

void soft_thresh()
{
	for (int i = 0; i < nV; i++)
	{
		double t = fabs(w[i])-lambda1;

		if(t<0)
		{
			t = 0;
		}

		if(w[i]>=0)
		{		
			w[i] = t;
		}
		else
		{
			w[i] = -t;
		}
	} 
}



void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) 
{
	// check input and output
	if (nrhs != 9)
		mexErrMsgTxt("Bad number of inputs arguments");

	if (nlhs != 2) 
		mexErrMsgTxt("Bad number of output arguments");

	
	// get data from input
	nV = (int)(*(double *)mxGetPr(prhs[0]));
	z = (double *)mxGetPr(prhs[1]);
	lambda1 = *(double *)mxGetPr(prhs[2]);
	lambda2 = *(double *)mxGetPr(prhs[3]);
	nEdge = (int)(*(double *)mxGetPr(prhs[4]));
	Edge_in = (double *)mxGetPr(prhs[5]);
	Edge_out = (double *)mxGetPr(prhs[6]);
	Edge_weight = (double *)mxGetPr(prhs[7]);
	flag = *(double*)mxGetPr(prhs[8]);

	// allocate t and w
	t = new double [nV];
	w = new double [nV];

	///////// sovle FLSA /////////////////

		// calc mnp: min_{t\in B(g)} { 0.5*|t|_2^2 }
		solve_mnp();

		// calc w from t;
		calcw();

		// soft-thresholding
		if (lambda1!=0)
		{
			soft_thresh();
		}
		
		// calc Omega(w)=\sum{ |w_i-w_j| }
		fw = 0;
		if (flag)
		{
			fw = calcOmega();
		}
	//////////////////////////////////////


	// output
	plhs[0] = mxCreateDoubleMatrix(nV,1,mxREAL); 
	double * pw = (double *)mxGetPr(plhs[0]); 
	for(int i = 0; i < nV; i++)
	{
		pw[i] = w[i]; 
	} 
	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
	double * pf = (double *)mxGetPr(plhs[1]);
	pf[0] = (double)fw;


	// release memo
	delete [] t;
	delete [] w;

}



