#include "efu.hpp"

Pairwise::Pairwise(int len, int q, tens3 coup, tens2 fields ) : len(len), q(q), coup(coup), fields(fields) {
	
}

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
		en-=fields[state.seq[i]][i];
		for (int j=i+1; j<len; ++j)
			en-=coup[state.seq[i]][state.seq[j]][++l];
	}
	return en;
}

int ind2(int i, int j,int N){
	return j - 1 + i*(2*N-3-i)/2;
}

double Pairwise::get_move_endiff(State& state){
	int i = state.pos_prop;
	double endiff=fields[ state.seq[i] ][ i ] - fields[ state.color_prop ][ i ];
	// we start at the l-index for (0,i) 
	// we add state.len-2-state.pos_prop to go from (i,j) to (i+1,j) for any i<j 
	int l=ind2(0,i,state.len);
	for (int j=0; j<i; ++j){
		endiff+=coup[state.seq[j]][state.seq[i]][l];
		l+=state.len-2-i;
	}
	// we start at the l-index for (i,i+1)
	// add 1 to the l-index of (i,j) to (i,j+1)
	l=ind2(i,i+1,state.len);
	for (int j=i+1; j<state.len-1; ++j){
		endiff+=coup[state.seq[i]][state.seq[j]][l];
		l++;
	}
	return endiff;
}


