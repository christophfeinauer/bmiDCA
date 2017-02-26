#ifndef ENV_H
#define ENV_H
#include "state.hpp"
#include "efu.hpp"
#include <vector>
class Env{
	double beta;
	Efu  * efu;

	//RNG generator
	std::mt19937_64 mtgen;
	//Becomes RNG for [0,1]
	std::uniform_real_distribution<double> * rdist;
	public:	
		State * state;
		Env(State*, Efu*, int);
		bool decide_move();
		bool step();
		double get_energy() {return efu->get_energy(*state);};
		~Env();
};
#endif
