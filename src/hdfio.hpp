#ifndef HDFIO
#define HDFIO
#include "hdf5.h"
#include <string>
#include "boost/multi_array.hpp"
typedef long long unsigned llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;

void readtens3(std::string,std::shared_ptr<tens3>&,std::string, std::string="fortran");
void readtens2(std::string,std::shared_ptr<tens2>&,std::string, std::string="fortran");
void writetens3(std::string,std::shared_ptr<tens3>&,std::string,std::string);
void addtens3(std::string,std::shared_ptr<tens3>&,std::string,std::string);
void writetens2(std::string,std::shared_ptr<tens2>&,std::string,std::string);
void addtens2(std::string,std::shared_ptr<tens2>&,std::string,std::string);

#endif
