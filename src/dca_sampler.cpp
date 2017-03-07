#include "env.hpp"
#include "state.hpp"
#include "efu.hpp"
#include "mcrun.hpp"
#include <stdio.h>
#include <stdlib.h>

int main(void){
	srand(time(NULL));
	std::string infile = "/home/christoph/delete_me/PF00014_full.txt.plmdca.hdf5";
	std::string ofile = "test";
	MCMCRun mcrun(infile,10000,1000000,100000,ofile);
	mcrun.run();
	return 0;
}
