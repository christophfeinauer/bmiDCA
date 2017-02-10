#ifndef ENV_H
#define ENV_H
#include "state.hpp"
#include <vector>
class Env{
	double beta;
	std::vector<std::vector<std::vector<double> > > coup;
	std::vector<std::vector<double> > fields;
	State state;

	//RNG generator
	std::mt19937_64 mtgen;
	//Becomes RNG for [0,1]
	std::uniform_real_distribution<double> * rdist;
	public:	
		Env(State);
		~Env();
};
#endif
