#include "efu.hpp"
#include "hdfio.hpp"
#include <stdexcept>
#include <string>

typedef long long unsigned llu;

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

Efu::Efu(std::string storage_order) : storage_order(storage_order) {};

Pairwise::Pairwise(int len, int q,std::string storage_order) : Efu(storage_order), len(len), q(q) {

	llu lenbn2 = (len*(len-1))/2;
	if (storage_order == "fortran"){
		coup_ptr.reset(new tens3(boost::extents[q][q][lenbn2],boost::fortran_storage_order()));
		fields_ptr.reset(new tens2(boost::extents[q][len],boost::fortran_storage_order()));
	}
	else{
		coup_ptr.reset(new tens3(boost::extents[q][q][lenbn2]));
		fields_ptr.reset(new tens2(boost::extents[q][len]));
	}
	std::fill_n(coup_ptr->data(),coup_ptr->num_elements(), 0.0);
	std::fill_n(fields_ptr->data(),fields_ptr->num_elements(), 0.0);
	

}

Pairwise::Pairwise(std::string fn, std::string coup_name, std::string fields_name, std::string storage_order) : Efu(storage_order) {

	readtens3(fn,coup_ptr,coup_name,storage_order);
	readtens2(fn,fields_ptr,fields_name,storage_order);
	q = (int) fields_ptr->shape()[0];
	len = (int) fields_ptr->shape()[1];
	if (q!=coup_ptr->shape()[0] || q!=coup_ptr->shape()[1] || (len*(len-1))/2 != coup_ptr->shape()[2])
		throw ReadError("Coupling and field dimensions non consistent");
}

double Pairwise::get_energy(State const & state){
	double en=0.0;
	long int l=0;
	for (int i=0; i<len; ++i){
		en-=(*fields_ptr)[state.seq[i]][i];
		for (int j=i+1; j<len; ++j){
			en-=(*coup_ptr)[state.seq[i]][state.seq[j]][l];
			l++;
	}
	}
	return en;
}

inline int ind2(int i, int j,int N){
	return j - 1 + i*(2*N-3-i)/2;
}

double Pairwise::get_move_endiff(State const & state){
	tens3& coup = *coup_ptr;
	tens2& fields = *fields_ptr;
	int i = state.pos_prop;
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

void Pairwise::add(std::vector<double> & to_add){
	double* fields_flat = fields_ptr->data();
	double* coup_flat = coup_ptr->data();
	llu k = 0;
	for (llu i=0; i<fields_ptr->num_elements(); ++i)
		fields_flat[i]=fields_flat[i]+to_add[k++];
	for (llu i=0; i<coup_ptr->num_elements(); ++i)
		coup_flat[i]=coup_flat[i]+to_add[k++];
	
}

void Pairwise::subtract_with_factor(std::vector<double> & g,double lambda){
	double* fields_flat = fields_ptr->data();
	double* coup_flat = coup_ptr->data();
	llu k = 0;
	for (llu i=0; i<fields_ptr->num_elements(); ++i)
		g[k] = g[k] - lambda*fields_flat[i];
	for (llu i=0; i<coup_ptr->num_elements(); ++i)
		g[k] = g[k] - lambda*coup_flat[i];
	

}


