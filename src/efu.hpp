#ifndef EFU_H
#define EFU_h
#include "state.hpp"
#include <vector>
typedef std::vector<std::vector<std::vector<double> > > tens3;
typedef std::vector<std::vector<double> > tens2;


// pure virtual base class for energy functions
class Efu{
	public:
		virtual double get_energy(State& state) = 0;
		virtual double get_move_energy_diff(State& state) = 0;
};


class Pairwise : public Efu{
	tens3 coup;
	tens2 fields;
	int len;
	int q;
	public:
		double get_energy(State& state);
		double get_move_endiff(State& state);
		Pairwise(int, int);
		Pairwise(int,int,tens3,tens2);

};
#endif
