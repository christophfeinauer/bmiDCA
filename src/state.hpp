#ifndef STATE_H
#define STATE_H
#include <vector>
#include <random>
class State{
	// contains variable positions
	std::vector<int> varpos;
	// RNG generator
	std::mt19937_64 mtgen;
	// RNG that chooses position in a move
	std::uniform_int_distribution<int> * pos_dist;
	// RNG that chooses color in a move
	std::uniform_int_distribution<int> * color_dist;
	// Keep track of steps
	unsigned long long eq_tries, eq_acc, batch_tries, batch_acc, all_tries, all_acc;
	public:
		State(int,int,int);
		State(const State&);
		~State();
		void propose_move();
		// variables for proposed moves
		int pos_prop,color_prop = 0;
		int len,q;
		// contains state sequence
		std::vector<int> seq;

};
#endif
