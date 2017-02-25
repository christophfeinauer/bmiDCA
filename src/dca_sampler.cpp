#include "env.hpp"
#include "state.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	std::string fn = "/home/christoph/.dropbox-alt/Dropbox (HuGeF)/varun_shachi/FGF/alignment.fa.numbered.plmdca.hdf5";
	Pairwise pw = Pairwise(fn,"J","h");
	State state(pw.len,pw.q,rand());
	Env env = Env(&state,&pw, rand());
	for (int k=0; k<10; ++k){
		for (int i=0; i<10000; ++i){
			env.step();
		}
		printf("%lf\n",pw.get_energy(state));
		for (int n=0; n<pw.len; ++n)
			printf("%d ",state.seq[n]);
		printf("\n");
	}
	printf("%lf\n",state.acc());
	return 0;
}
