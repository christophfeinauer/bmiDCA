#ifndef ENV_H
#define ENV_H
#include "state.hpp"
#include "efu.hpp"
#include <vector>
class Env{
	double beta;
	//RNG generator
	std::mt19937_64 mtgen;
	//Becomes RNG for [0,1] and an interger seed for the states
	std::uniform_real_distribution<double>  rdist_double;
	std::uniform_int_distribution<int>  rdist_int;
	public:	
		std::shared_ptr<Efu> efu;
		int nr_states;
		std::vector<State> state_vec;
		//Env(State, Efu*, int);
		Env(std::shared_ptr<Efu>, int, int, int, int);
		Env(int, int, int, std::string="fortran");
		//Functions operating on state_ptr_vec[0]
		bool decide_move();
		void reset_moves();
		void reset_moves(int);
		void reset_all_moves();
		void reset_all_moves(int);
		bool step();
		double get_energy() {return efu->get_energy(state_vec[0]);};
		//Functions operating on state_ptr_vec[i]
		bool decide_move(int);
		bool step(int);
		double get_energy(int i) {return efu->get_energy(state_vec[i]);};
		std::vector<double> get_energies();
		std::vector<State>* get_ptr_state_vec() {return &state_vec;};
		
};
#endif
