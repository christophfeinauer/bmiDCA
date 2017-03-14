#include "bmirun.hpp"
#include "env.hpp"


int main(void){

	llu nsamples = 500;
	llu eq_steps = 100000;
	llu delta_steps = 10000;
	llu max_runs = 100000;
	double lambda = 0.01;
	double epsilon = 0.1;
	double alpha = 0.5;

	Simple bmirun("/home/christoph/Dropbox/delete_me/PF00014_full.txt.theta_0.2.f.hdf5",epsilon,nsamples,eq_steps,delta_steps,lambda,alpha,max_runs);
	bmirun.run();
}
