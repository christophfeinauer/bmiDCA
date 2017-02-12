#ifndef ENV_H
#define ENV_H
#include "state.hpp"
#include "efu.hpp"
#include <vector>
class Env{
	double beta;
	State * state;
	Efu  * efu;

	//RNG generator
	std::mt19937_64 mtgen;
	//Becomes RNG for [0,1]
	std::uniform_real_distribution<double> * rdist;
	public:	
		Env(State*, Efu*);
		~Env();
};
#endif
