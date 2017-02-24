#include "env.hpp"
#include "state.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	int len=10;
	int q=21;
	State state(len,q,rand());
	std::string fn = "/home/christoph/upmc/varun_shachi/FGF/alignment.fa.numbered.plmdca.hdf5";
	Pairwise pw = Pairwise(fn,"J","h");
	//Env env = Env(&state,&pw, rand());
	//for (int i=0; i<1000000; ++i){
	//	env.step();
	//}
	//printf("%lf\n",state.acc());
	return 0;
}
