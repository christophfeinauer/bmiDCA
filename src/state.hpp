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
	std::uniform_int_distribution<int> pos_dist;
	// RNG that chooses color in a move
	std::uniform_int_distribution<int> color_dist;
	// Keep track of steps
	unsigned long long moves_proposed = 0, moves_made = 0;
	public:
		State(int,int,int);
		void propose_move();
		void make_move();
		void reset_moves();
		void reset_all_moves();
		double acc();
		// variables for proposed moves
		int pos_prop,color_prop = 0;
		int len,q;
		unsigned long long moves_proposed_total = 0, moves_made_total = 0;
		// contains state sequence
		std::vector<int> seq;

};
#endif
