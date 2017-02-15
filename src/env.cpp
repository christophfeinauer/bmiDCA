#include "env.hpp"
// Constructor that produces a random environment with beta=1
Env::Env(State * state_ptr , Efu * efu) : state(state_ptr) , efu(efu){
	beta=1.0;
	rdist = new std::uniform_real_distribution<double>;
}

Env::~Env(){
	delete rdist;
}

bool Env::decide_move(){
	double d = (*rdist)(mtgen);
	printf("%lf\n",d);
	return d > exp(efu->get_move_endiff(*state));
	//return (*rdist)(mtgen) > exp(efu->get_move_endiff(*state));
}

bool Env::step(){
	state->propose_move();
	if (decide_move())
		state->make_move();
}
