#include <vector>
#include <stdio.h>
#include <random>
#include <time.h>
#include "state.hpp"

class State{
	int len,q;
	// contains state sequence
	std::vector<int> seq;
	// contains variable positions
	std::vector<int> varpos;
	// RNG generator
	std::mt19937_64 mtgen;
	// RNG that chooses position in a move
	std::uniform_int_distribution<int> * pos_dist;
	// RNG that chooses color in a move
	std::uniform_int_distribution<int> * color_dist;
	// Keep track of steps
	unsigned long long eq_tries, eq_acc, batch_tries, batch_acc, all_tries, all_acc;
	public:
		State(int,int,int);
		~State();
		void propose_move();
		// variables for proposed moves
		int pos_prop,color_prop = 0;
};

State::State(int _len, int _q, int seed){
	len = _len;
	q = _q;
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

int main(void){
	State state(10,21,100);
	for (int i = 0; i<100; ++i){
		state.propose_move();
		printf("%d %d\n",state.pos_prop,state.color_prop);
	}
	return 0;
}
