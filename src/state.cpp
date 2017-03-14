#include <stdio.h>
#include <time.h>
#include "state.hpp"

State::State(int len, int q, int seed) : len(len), q(q), pos_dist(0,len-1), color_dist(0,q-1){
	seq.resize(len,0);
	// Standard constructor: All positions are variable
	varpos.resize(len);
	for (int pos=0; pos<len; ++pos)
		varpos[pos]=pos;
	// Seed the generator
	mtgen.seed(seed);
}

double State::acc(){
	return (double) moves_made/ (double) moves_proposed;
}

void State::propose_move(){
	pos_prop = pos_dist(mtgen);
	color_prop = color_dist(mtgen);
	moves_proposed++;
	moves_proposed_total++;
}

void State::make_move(){
	seq[pos_prop] = color_prop;
	moves_made++;
	moves_made_total++;
}

void State::reset_moves(){
	moves_made = 0;
	moves_proposed = 0;
}

void State::reset_all_moves(){
	moves_made = 0;
	moves_proposed = 0;
	moves_made_total=0;
	moves_proposed_total=0;
}
