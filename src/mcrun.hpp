#ifndef MCRUN_H
#define MCRUN_H
#include "env.hpp"
typedef long long int llu;
typedef std::vector<std::vector<std::vector<double> > > tens3;
typedef std::vector<std::vector<double> > tens2;

class MCRun{
	public:
		Env * env;
		llu nsamples;
		void virtual run() = 0;
		MCRun(Env*,llu);
		MCRun(Env*,llu,tens3*);
		tens3 * f_tens3 = NULL;
		llu f_tens3_count = 0;

};

class MCMCRun : MCRun {
	public: 
		llu eq_steps;
		llu delta_steps;	
		std::string ofile;
		MCMCRun(Env*,llu, llu , llu, std::string);
		void run();
		void write_state_to_file();
		void write_state_to_f_tens3();
};
#endif
