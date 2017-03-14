#ifndef EFU_H
#define EFU_H
#include "state.hpp"
#include <vector>
#include <string>
#include "boost/multi_array.hpp"
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;
typedef std::shared_ptr<double> double_sarray_ptr;
typedef long long unsigned llu;


// pure virtual base class for energy functions
class Efu{
	public:
		virtual double get_energy(State const & state) = 0;
		virtual double get_move_endiff(State const & state) = 0;
		llu nparameters;
		virtual void add(std::vector<double>&) = 0;
		virtual void subtract_with_factor(std::vector<double>&,double) = 0;
		std::string storage_order = "fortran";
		Efu(std::string);
};


class Pairwise : public Efu{
	tens3_ptr coup_ptr;
	tens2_ptr fields_ptr;
	public:
		int len;
		int q;
		double get_energy(State const & state);
		double get_move_endiff(State const & state);
		void add(std::vector<double>&);
		void subtract_with_factor(std::vector<double>&, double);
		//Pairwise(int, int,int);
		Pairwise(int,int,std::string="fortran");
		Pairwise(std::string,std::string,std::string, std::string="fortran");

};
#endif
