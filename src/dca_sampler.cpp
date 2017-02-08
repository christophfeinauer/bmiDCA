#include "environ.hpp"
#include "state.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	State state(10,21,rand());
	Env env = Env(state);
	for (int i = 0; i<100; ++i){
		state.propose_move();
		printf("%d %d\n",state.pos_prop,state.color_prop);
	}
	return 0;
}
