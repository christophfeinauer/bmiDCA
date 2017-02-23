#include "efu.hpp"

Pairwise::Pairwise(int len, int q, tens3 coup, tens2 fields ) : len(len), q(q), coup(coup), fields(fields) {
	
}

Pairwise::Pairwise(std::string fn, std::string coup_name, std::string fields_name){
	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);		
}

// constructor using random couplings
Pairwise::Pairwise(int len, int q, int seed) : len(len), q(q) {
	
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
	// Random fields
	for (int c=0; c<q; ++c)
		for (int n=0; n<len; ++n)
			fields[c][n] =  2*rand()-1.0;

	// Random couplings
	for (int l=0; l<nsp; ++l)
		for (int a=0; a<q; ++a)
			for (int b=0; b<q; ++b)
				coup[a][b][l] = 2*rand()-1.0;

}

double Pairwise::get_energy(State& state){
	double en=0.0;
	long int l=0;
	for (int i=0; i<len; ++i){
		en-=fields[state.seq[i]][i];
		for (int j=i+1; j<len; ++j){
			en-=coup[state.seq[i]][state.seq[j]][l];
			l++;
	}
	}
	return en;
}

int ind2(int i, int j,int N){
	return j - 1 + i*(2*N-3-i)/2;
}

double Pairwise::get_move_endiff(State& state){
	int i = state.pos_prop;
	//printf("%d\n",i);
	double endiff=fields[ state.seq[i] ][ i ] - fields[ state.color_prop ][ i ];
	// we start at the l-index for (0,i) 
	// we add state.len-2-state.pos_prop to go from (j,i) to (j+1,i) for any j<i
	int l=ind2(0,i,state.len);
	for (int j=0; j<i; ++j){
		endiff += coup[state.seq[j]][state.seq[i]][l] - coup[state.seq[j]][state.color_prop][l];
		l+=state.len-2-j;
	}
	// we start at the l-index for (i,i+1)
	// add 1 to the l-index of (i,j) to (i,j+1)
	l=ind2(i,i+1,state.len);
	for (int j=i+1; j<state.len; ++j){
		endiff += coup[state.seq[i]][state.seq[j]][l] - coup[state.color_prop][state.seq[j]][l];
		l++;
	}
	return endiff;
}


