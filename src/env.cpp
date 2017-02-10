#include "env.hpp"
// Constructor that produces a random environment with beta=1
Env::Env(State state) : state(state){
	beta=1.0;
	// Resize parameter vectors
	fields.resize(state.q);
	for (int c=0; c<state.q;++c)
		fields[c].resize(state.len);
	coup.resize(state.q);
	// Number of site pairs
	int nsp = state.len * (state.len-1)/2;
	for (int a=0; a<state.q; ++a){
		coup[a].resize(state.q);
			for (int b=0; b<state.q; ++b)
				coup[a][b].resize(nsp);
	}
	// Make RNG
	rdist = new std::uniform_real_distribution<double>(0.0,1.0);
	// Random fields
	for (int c=0; c<state.q; ++c)
		for (int n=0; n<state.len; ++n)
			fields[c][n] = 2*(*rdist)(mtgen)-1.0;
	// Random couplings
	for (int l=0; l<nsp; ++l)
		for (int a=0; a<state.q; ++a)
			for (int b=0; b<state.q; ++b)
				coup[a][b][l] = 2*(*rdist)(mtgen)-1.0;
}

Env::~Env(){
	printf("Destroying Environment\n");
	delete rdist;
}
