#include "efu.hpp"
#include "hdf5.h"
#include <stdexcept>
#include <string>

typedef long long unsigned llu;

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

Pairwise::Pairwise(int len, int q, tens3_ptr coup_ptr, tens2_ptr fields_ptr ) : len(len), q(q), coup_ptr(coup_ptr), fields_ptr(fields_ptr) {
	
}

//ATTENTION: This function expects a qxqxbinomial(N,2) column-major memory layout for couplings
Pairwise::Pairwise(std::string fn, std::string coup_name, std::string fields_name){
	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);		
	// read couplings
	hid_t dset = H5Dopen1(hfid,coup_name.c_str());
	hid_t dspace = H5Dget_space(dset);
	int ndims_coup = H5Sget_simple_extent_ndims(dspace);
	if (ndims_coup!=3)
		throw ReadError("HDF5Error: Dimensions of couplings ≠ 3");
	printf("Reading couplings from HDF5 data..."); 
	std::vector<hsize_t> dims_coup(ndims_coup,0);
	H5Sget_simple_extent_dims(dspace,&dims_coup[0],NULL);
	llu lenbn2 = (llu) dims_coup[0];
	q = (int) dims_coup[1];
	len = (int) (1+sqrt(1+8*lenbn2))/2;
	if (dims_coup[2] != q )
		throw ReadError("HDF5Error: Dimensions not consistent");
	// the init is necessary because tens3_ptr coup_ptr(...) shadows the actual member
	tens3_ptr coup_ptr_init(new tens3(boost::extents[q][q][lenbn2],boost::fortran_storage_order()));
	coup_ptr = coup_ptr_init;
	H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,coup_ptr->data());
	H5Sclose(dspace);
	H5Dclose(dset);
	printf("done\n");
	// read fields
	dset = H5Dopen1(hfid,fields_name.c_str());
	dspace = H5Dget_space(dset);
	int ndims_fields = H5Sget_simple_extent_ndims(dspace);
	if (ndims_fields!=2)
		throw ReadError("HDF5Error: Dimensions of fields ≠ 2");
	printf("Reading fields from HDF5 data..."); 
	std::vector<hsize_t> dims_fields(ndims_fields,0);
	H5Sget_simple_extent_dims(dspace,&dims_fields[0],NULL);
	if (dims_fields[0] != len || dims_fields[1] != q)
		throw ReadError("HDF5Eror: Coupling- and fields-dimensions not consistent");
	tens2_ptr fields_ptr_init(new tens2(boost::extents[q][len],boost::fortran_storage_order()));
	fields_ptr = fields_ptr_init;
	H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,fields_ptr->data());
	H5Sclose(dspace);
	H5Dclose(dset);
	H5Fclose(hfid);
	printf("done\n");
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


