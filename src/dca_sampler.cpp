#include "env.hpp"
#include "state.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	int len=10;
	int q=21;
	State state(len,q,rand());
	Pairwise pw = Pairwise(len,q,rand());
	Env env = Env(&state,&pw);
	double n,o;
	for (int i=0; i<10000; ++i){
		env.step();
		n = pw.get_energy(state);
		if (n!=o){
			printf("%lf\n",n);
			o=n;
		}
	}
	return 0;
}
