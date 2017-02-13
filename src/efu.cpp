#include "efu.hpp"

// constructor using random couplings
Pairwise::Pairwise(int len, int q) : len(len), q(q) {
	
	// Resize parameter vectors
	fields.resize(q);
	for (int c=0; c<q;++c)
		fields[c].resize(len);
	coup.resize(q);
	// Number of site pairs
	int nsp = len * (len-1)/2;
	for (int a=0; a<q; ++a){
		coup[a].resize(q);
			for (int b=0; b<q; ++b)
				coup[a][b].resize(nsp);
	}
	std::uniform_real_distribution<double> rdist(0.0,1.0);
	std::mt19937_64 mtgen;
	// Random fields
	for (int c=0; c<q; ++c)
		for (int n=0; n<len; ++n)
			fields[c][n] = 2*rdist(mtgen)-1.0;
	// Random couplings
	for (int l=0; l<nsp; ++l)
		for (int a=0; a<q; ++a)
			for (int b=0; b<q; ++b)
				coup[a][b][l] = 2*rdist(mtgen)-1.0;

}

double Pairwise::get_energy(State& state){
	double en=0;
	long int l=0;
	for (int i=0; i<len; ++i){
		en+=fields[state.seq[i]][i];
		for (int j=i+1; j<len; ++j)
			en+=coup[state.seq[i]][state.seq[j]][l];
		l++;
	}
	return en;
}

double Pairwise::get_move_energy_diff(State& state){
	double en=0;
	long int l=0;
	for (int i=0; i<len; ++i){
		en+=fields[state.seq[i]][i];
		for (int j=i+1; j<len; ++j)
			en+=coup[state.seq[i]][state.seq[j]][l];
		l++;
	}
	return en;
}


