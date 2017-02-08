#include <stdio.h>
#include <time.h>
#include "state.hpp"

State::State(int len, int q, int seed) : len(len), q(q){
	seq.resize(len);
	// Standard constructor: All positions are variable
	varpos.resize(len);
	for (int pos=0; pos<len; ++pos)
		varpos[pos]=pos;
	// Seed the generator
	mtgen.seed(seed);
	// RNG that proposes position
	pos_dist = new std::uniform_int_distribution<int>(0,varpos.size()-1);
	// RNG that proposes color
	color_dist = new std::uniform_int_distribution<int>(0,q-1);
}

State::~State(){
	delete pos_dist;
	delete color_dist;
}

void State::propose_move(){
	pos_prop = (*pos_dist)(mtgen);
	color_prop = (*color_dist)(mtgen);
}
