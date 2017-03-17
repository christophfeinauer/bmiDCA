#include "hdfio.hpp"
#include <boost/filesystem.hpp>

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

struct WriteError : public std::runtime_error{
	WriteError(std::string const & message)
		: std::runtime_error(message)
	{}
};

void readtens3(std::string fn,std::shared_ptr<tens3>& tens3_ptr,std::string tens3_name, std::string storage_order){

	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);		
	hid_t dset = H5Dopen1(hfid,tens3_name.c_str());
	hid_t dspace = H5Dget_space(dset);
	int ndims_tens3 = H5Sget_simple_extent_ndims(dspace);
	if (ndims_tens3!=3)
		throw ReadError("HDF5Error: Dimensions of tens3 ≠ 3");
	printf("Reading %s from %s...",tens3_name.c_str(),fn.c_str()); 
	std::vector<hsize_t> dims_tens3(ndims_tens3,0);
	H5Sget_simple_extent_dims(dspace,&dims_tens3[0],NULL);
	// todo:: this actually expects fortran storage order
	llu lenbn2 = (llu) dims_tens3[0];
	int q = (int) dims_tens3[1];
	int len = (int) (1+sqrt(1+8*lenbn2))/2;
	if (dims_tens3[2] != q )
		throw ReadError("HDF5Error: Dimensions not consistent");
	if (storage_order=="fortran")
		tens3_ptr.reset(new tens3(boost::extents[q][q][lenbn2],boost::fortran_storage_order()));
	else
		tens3_ptr.reset(new tens3(boost::extents[q][q][lenbn2]));
	H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens3_ptr->data());
	H5Sclose(dspace);
	H5Dclose(dset);
	H5Fclose(hfid);
	printf("done\n");

}
void readtens2(std::string fn,std::shared_ptr<tens2>& tens2_ptr,std::string tens2_name, std::string storage_order){

	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);		
	hid_t dset = H5Dopen1(hfid,tens2_name.c_str());
	hid_t dspace = H5Dget_space(dset);
	int ndims_tens2 = H5Sget_simple_extent_ndims(dspace);
	if (ndims_tens2!=2)
		throw ReadError("HDF5Error: Dimensions of fields ≠ 2");
	printf("Reading %s from %s...",tens2_name.c_str(),fn.c_str()); 
	std::vector<hsize_t> dims_tens2(ndims_tens2,0);
	H5Sget_simple_extent_dims(dspace,&dims_tens2[0],NULL);
	// todo: this actually expects fortran storage order
	int q = (int) dims_tens2[1];
	int len = (int) dims_tens2[0];
	if (storage_order=="fortran")
		tens2_ptr.reset(new tens2(boost::extents[q][len],boost::fortran_storage_order()));
	else
		tens2_ptr.reset(new tens2(boost::extents[q][len]));
	H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens2_ptr->data());
	H5Sclose(dspace);
	H5Dclose(dset);
	H5Fclose(hfid);
	printf("done\n");
}

void writetens3(std::string fn,std::shared_ptr<tens3>& tens3_ptr,std::string tens3_name, std::string storage_order){
        
        if (tens3_ptr->num_dimensions()!=3)
                throw WriteError("dimensions of tens3 not 3");
        hid_t hfid = H5Fcreate(fn.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);              
        hsize_t dims[3]; 
        if (storage_order=="fortran")
                for (int i=0; i<3; ++i)
                        dims[i] = tens3_ptr->shape()[2-i];
        else
                for (int i=0;i<3;++i)
                        dims[i] = tens3_ptr->shape()[i];
        hid_t dspace = H5Screate_simple(3,dims,NULL);
        hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
        herr_t status = H5Pset_layout (dcpl, H5D_CONTIGUOUS);
        hid_t dset = H5Dcreate (hfid, tens3_name.c_str(), H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, dcpl,H5P_DEFAULT);
        printf("Writing %s to %s...",tens3_name.c_str(),fn.c_str()); 
        H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens3_ptr->data());
        H5Pclose(dcpl);
        H5Sclose(dspace);
        H5Dclose(dset);
        H5Fclose(hfid);
        printf("done\n");

}


void addtens3(std::string fn,std::shared_ptr<tens3>& tens3_ptr,std::string tens3_name, std::string storage_order){
	
	if (tens3_ptr->num_dimensions()!=3)
		throw WriteError("dimensions of tens3 not 3");
	
	if (!boost::filesystem::exists( fn )) throw WriteError("trying to add to non-existing file");

	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);		
	hsize_t dims[3]; 
	if (storage_order=="fortran")
		for (int i=0; i<3; ++i)
			dims[i] = tens3_ptr->shape()[2-i];
	else
		for (int i=0;i<3;++i)
			dims[i] = tens3_ptr->shape()[i];
	hid_t dspace = H5Screate_simple(3,dims,NULL);
	hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
	herr_t status = H5Pset_layout (dcpl, H5D_CONTIGUOUS);
	hid_t dset = H5Dcreate (hfid, tens3_name.c_str(), H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, dcpl,H5P_DEFAULT);
	printf("Writing %s to %s...",tens3_name.c_str(),fn.c_str()); 
	H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens3_ptr->data());
	H5Pclose(dcpl);
	H5Sclose(dspace);
	H5Dclose(dset);
	H5Fclose(hfid);
	printf("done\n");

}

void writetens2(std::string fn,std::shared_ptr<tens2>& tens2_ptr,std::string tens2_name, std::string storage_order){
        
        if (tens2_ptr->num_dimensions()!=2)
                throw WriteError("dimensions of tens2 not 2");
        hid_t hfid = H5Fcreate(fn.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);              
        hsize_t dims[2]; 
        if (storage_order=="fortran")
                for (int i=0; i<2; ++i)
                        dims[i] = tens2_ptr->shape()[1-i];
        else
                for (int i=0;i<2;++i)
                        dims[i] = tens2_ptr->shape()[i];
        hid_t dspace = H5Screate_simple(2,dims,NULL);
        hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
        herr_t status = H5Pset_layout (dcpl, H5D_CONTIGUOUS);
        hid_t dset = H5Dcreate (hfid, tens2_name.c_str(), H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, dcpl,H5P_DEFAULT);
        printf("Writing %s to %s...",tens2_name.c_str(),fn.c_str()); 
        H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens2_ptr->data());
        H5Pclose(dcpl);
        H5Sclose(dspace);
        H5Dclose(dset);
        H5Fclose(hfid);
        printf("done\n");

}


void addtens2(std::string fn,std::shared_ptr<tens2>& tens2_ptr,std::string tens2_name, std::string storage_order){
	
	if (tens2_ptr->num_dimensions()!=2)
		throw WriteError("dimensions of tens2 not 2");
	
	if (!boost::filesystem::exists( fn )) throw WriteError("trying to add to non-existing file");

	hid_t hfid = H5Fopen(fn.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);		
	hsize_t dims[2]; 
	if (storage_order=="fortran")
		for (int i=0; i<2; ++i)
			dims[i] = tens2_ptr->shape()[1-i];
	else
		for (int i=0;i<2;++i)
			dims[i] = tens2_ptr->shape()[i];
	hid_t dspace = H5Screate_simple(2,dims,NULL);
	hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
	herr_t status = H5Pset_layout (dcpl, H5D_CONTIGUOUS);
	hid_t dset = H5Dcreate (hfid, tens2_name.c_str(), H5T_NATIVE_DOUBLE, dspace, H5P_DEFAULT, dcpl,H5P_DEFAULT);
	printf("Writing %s to %s...",tens2_name.c_str(),fn.c_str()); 
	H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,tens2_ptr->data());
	H5Pclose(dcpl);
	H5Sclose(dspace);
	H5Dclose(dset);
	H5Fclose(hfid);
	printf("done\n");

}

