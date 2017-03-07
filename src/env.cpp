#include "env.hpp"
//simple constructor with all len/q/beta the same
Env::Env(std::shared_ptr<Efu> efu, int nr_states,int len, int q, int seed) : efu(efu), nr_states(nr_states), rdist_double(), rdist_int(){
	mtgen.seed(seed);
	for (int i=0; i<nr_states; ++i){
		state_vec.push_back(State(len,q,rdist_int(mtgen)));
	}
	beta=1.0;
}

std::vector<double> Env::get_energies(){
	std::vector<double> envec(nr_states,0.0);
	for (int i=0; i<nr_states; ++i)
		envec[i]=get_energy(i);
	return envec;
}

bool Env::decide_move(){
	double d = rdist_double(mtgen);
	return d < exp(-efu->get_move_endiff(state_vec[0]));
}

bool Env::step(){
	state_vec[0].propose_move();
	if (decide_move()){
		state_vec[0].make_move();
		return true;}
	return false;	
}

bool Env::decide_move(int i){
	double d = rdist_double(mtgen);
	return d < exp(-efu->get_move_endiff(state_vec[i]));
}

bool Env::step(int i){
	state_vec[i].propose_move();
	if (decide_move(i)){
		state_vec[i].make_move();
		return true;}
	return false;	
}

void Env::reset_moves(int i){
	state_vec[i].reset_moves();
}

void Env::reset_moves(){
	for (int i=0; i<state_vec.size(); ++i)
		state_vec[i].reset_moves();
}
