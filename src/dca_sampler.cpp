#include "env.hpp"
#include "state.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	int len=10;
	int q=21;
	State state(len,q,rand());
	Pairwise pw = Pairwise(len,q);
	Env env = Env(&state,&pw);
	//printf("%lf\n",pw.get_energy(state));
	
	return 0;
}
