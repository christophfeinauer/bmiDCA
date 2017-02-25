#ifndef EFU_H
#define EFU_h
#include "state.hpp"
#include <vector>
#include <string>
typedef std::vector<std::vector<std::vector<double> > > tens3;
typedef std::vector<std::vector<double> > tens2;


// pure virtual base class for energy functions
class Efu{
	public:
		virtual double get_energy(State& state) = 0;
		virtual double get_move_endiff(State& state) = 0;
};


class Pairwise : public Efu{
	tens3 coup;
	tens2 fields;
	public:
		int len;
		int q;
		double get_energy(State& state);
		double get_move_endiff(State& state);
		Pairwise(int, int,int);
		Pairwise(int,int,tens3,tens2);
		Pairwise(std::string,std::string,std::string);

};
#endif
