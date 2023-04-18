//(C) Microsoft corporation. All rights reserved
#ifndef Parser_h_SEEN
#define Parser_h_SEEN

#include "Graph.h"

namespace paraF
{

struct Problem
{
	enum ProblemType
	{
		MAX_FLOW_PROBLEM,
		PARAMETRIC_FLOW_PROBLEM
	};

	ProblemType type;

	bool hasLambdaRange;
	double minLambda;
	double maxLambda;
};

class ParserBase
{
protected:
	bool ReadLine(char* line);
};

class InputParser : public ParserBase
{
public:
	void Parse(Graph& graph, Problem& problem);

};

class Output
{
public:
	double multiplier;
	std::vector<std::pair<LongType, LongType> > breakpoints;
};

class OutputParser : public ParserBase
{
public:
	void Parse(Output& output);

};

}
#endif

