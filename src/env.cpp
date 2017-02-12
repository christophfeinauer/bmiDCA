#include "env.hpp"
// Constructor that produces a random environment with beta=1
Env::Env(State * state , Efu * efu) : state(state) , efu(efu){
	beta=1.0;
}

Env::~Env(){
	printf("Destroying Environment\n");
	delete rdist;
}

