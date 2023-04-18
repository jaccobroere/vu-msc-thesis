//(C) Microsoft Corporation. All rights reserved.
#include "stdafx.h"
#include <iostream>
#include "Parser.h"
#include "Graph.h"
using namespace std;

namespace paraF
{

const int MAX_LINE_LENGTH = 256;

void InputParser::Parse(Graph& graph, Problem& problem)
{
	char line[MAX_LINE_LENGTH];
	if (!ReadLine(line) || line[0] != 'p')
		throw "problem line expected";

	char problemSingature[MAX_LINE_LENGTH];
	int n;
	int m;

	if (sscanf(line, "%*c %s %d %d", problemSingature, &n, &m) < 3)
		throw "error in problem line format";

	if (strcmp(problemSingature, "par") == 0)
		problem.type = Problem::PARAMETRIC_FLOW_PROBLEM;
	else if (strcmp(problemSingature, "max") == 0)
		problem.type = Problem::MAX_FLOW_PROBLEM;
	else
		throw "wrong problem signature";

	problem.hasLambdaRange = false;

	graph.Reserve(m * 2 + n * 4, n);
	for (int i = 0; i < n; i++)
		graph.AddNode();

	int nArcsRead = 0;

	int sourceNodeIndex = -1;
	int sinkNodeIndex = -1;

	while (ReadLine(line))
	{
		switch (line[0])
		{
		case 'n':
			char nodeType;
			int nodeIndex;
			if (sscanf(line, "%*c %d %c", &nodeIndex, &nodeType) < 2)
				throw "error in node line format";
			if (nodeIndex < 1 || nodeIndex > n)
				throw "node index is out of range";
			nodeIndex--;
			switch (nodeType)
			{
			case 's':
				sourceNodeIndex = nodeIndex;
				break;
			case 't':
				sinkNodeIndex = nodeIndex;
				break;
			default:
				throw "unrecognized node type";
			}
			break;

		case 'a':
			{
				if (sourceNodeIndex < 0 || sinkNodeIndex < 0)
					throw "definition of source and sink nodes should precede arc definitions";
				if (nArcsRead >= m)
					throw "too many arcs in the input";

				nArcsRead++;

				int tailNodeIndex;
				int headNodeIndex;

				if (sscanf(line, "%*c %d %d", &tailNodeIndex, &headNodeIndex) < 2)
					throw "error in arc line format";

				if (tailNodeIndex < 1 || tailNodeIndex > n)
					throw "tail node index is out of range";
				if (headNodeIndex < 1 || headNodeIndex > n)
					throw "head node index is out of range";

				headNodeIndex--;
				tailNodeIndex--;

				Node* headNode = &graph.nodes[headNodeIndex];
				Node* tailNode = &graph.nodes[tailNodeIndex];

				double cap = 0;
				double slope = 0;

				Arc* arc;
				if (problem.type == Problem::PARAMETRIC_FLOW_PROBLEM &&
					(tailNodeIndex == sourceNodeIndex || headNodeIndex == sinkNodeIndex))
				{
					if (sscanf(line, "%*c %*d %*d %lf %lf", &cap, &slope) < 2)
						throw "error in arc line format";
					if (slope < 0)
						throw "slope should be non-negative";
                    
					if (tailNodeIndex == sourceNodeIndex && headNodeIndex == sinkNodeIndex)
						continue;

					if (tailNodeIndex == sourceNodeIndex)
						arc = graph.AddArcFromSource(tailNode, headNode, cap, slope);
					else
						arc = graph.AddArcToSink(tailNode, headNode, cap, -slope);
				}
				else
				{
					if (headNodeIndex == sourceNodeIndex)
						continue;
					if (tailNodeIndex == sinkNodeIndex)
						continue;

					if (sscanf(line, "%*c %*d %*d %lf", &cap) < 1)
						throw "error in arc line format";
                    
					if (tailNodeIndex == sourceNodeIndex)
						arc = graph.AddArcFromSource(tailNode, headNode, cap, 0);
					else if (headNodeIndex == sinkNodeIndex)
						arc = graph.AddArcToSink(tailNode, headNode, cap, 0);
					else
						arc = graph.AddArc(tailNode, headNode, cap);
				}

				Arc* revArc = graph.AddArc(headNode, tailNode, 0);
				arc->revArc = revArc;
				revArc->revArc = arc;

				break;
			}


		case 'r':
			{
				if (problem.type != Problem::PARAMETRIC_FLOW_PROBLEM)
					throw "lambda range line is only applicable to parametric problem";

				if (problem.hasLambdaRange)
					throw "duplicate lambda range line";

				if (sscanf(line, "%*c %d %d", &problem.minLambda, &problem.maxLambda) < 2)
					throw "error in lambda range line";

				if (problem.minLambda >= problem.maxLambda)
					throw "invalid lambda range";

				problem.hasLambdaRange = true;
				break;
			}

		default:
			throw "unrecognized line type";
		}
	}

	if (nArcsRead < m)
		throw "too few arcs in the input";
	if (sourceNodeIndex < 0)
		throw "source node is not chosen";
	if (sinkNodeIndex < 0)
		throw "sink node is not chosen";

	if (problem.type == Problem::PARAMETRIC_FLOW_PROBLEM)
		graph.AddAuxArcs(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);

	graph.PrepareNodes(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);
	graph.PrepareArcs();
}

void OutputParser::Parse(Output& output)
{
	char line[MAX_LINE_LENGTH];

	output.breakpoints.clear();
	output.multiplier = -1;

	while (ReadLine(line))
	{
		switch (line[0])
		{
		case 'm':
			if (sscanf(line, "%*c %d", &output.multiplier) < 1)
				throw "error in multiplier line format";
			break;
			
		case 'z':
			break;

		case 'l':
			{
				std::pair<LongType, LongType> breakpoint;
				if (sscanf(line, "%*c %" LONGTYPE_SPEC " %" LONGTYPE_SPEC, &breakpoint.first, &breakpoint.second) < 2)
					throw "error in breakpoint line format";
				output.breakpoints.push_back(breakpoint);
				break;
			}

		default:
			throw "unrecognized line type";
		}
	}

	if (output.multiplier < 0)
		throw "multiplier is not specified";
}

bool ParserBase::ReadLine(char *line)
{
	do
	{
		if (!fgets(line, MAX_LINE_LENGTH, stdin))
			return false;
	} while (line[0] == 'c' || line[0] == '\n' || line[0] == 0);
	return true;
}

}
