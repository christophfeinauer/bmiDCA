#ifndef EFU_H
#define EFU_H
#include "state.hpp"
#include <vector>
#include <string>
#include "boost/multi_array.hpp"
#include "boost/shared_ptr.hpp"
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef boost::shared_ptr<tens3> tens3_ptr;
typedef boost::shared_ptr<tens2> tens2_ptr;


// pure virtual base class for energy functions
class Efu{
	public:
		virtual double get_energy(State const & state) = 0;
		virtual double get_move_endiff(State const & state) = 0;
};


class Pairwise : public Efu{
	tens3_ptr coup_ptr;
	tens2_ptr fields_ptr;
	public:
		int len;
		int q;
		double get_energy(State const & state);
		double get_move_endiff(State const & state);
		//Pairwise(int, int,int);
		Pairwise(int,int,tens3_ptr,tens2_ptr);
		Pairwise(std::string,std::string,std::string);

};
#endif
