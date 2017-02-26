#include "env.hpp"
#include "state.hpp"
#include "efu.hpp"
#include "mcrun.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	std::string fn = "/home/christoph/delete_me/PF00014_full.txt.plmdca.hdf5";
	Pairwise pw = Pairwise(fn,"J","h");
	State state(pw.len,pw.q,rand());
	Env env = Env(&state,&pw, rand());

	std::string ofile = "test";
	MCMCRun mcrun(&env,50000,1000000,100000,ofile);
	mcrun.run();
	return 0;
}
