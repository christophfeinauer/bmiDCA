#ifndef MCRUN_H
#define MCRUN_H
#include "env.hpp"
typedef long long int llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;

class MCRun{
	public:
		std::shared_ptr<Env> env;
		llu nsamples;
		void virtual run() = 0;
		MCRun(std::shared_ptr<Env>,llu,tens3*);
		tens3 * f_tens3 = NULL;
		llu f_tens3_count = 0;
};

class MCMCRun : MCRun {
	public: 
		llu eq_steps;
		llu delta_steps;	
		std::string ofile;
		MCMCRun(std::shared_ptr<Env>,llu, llu , llu, std::string);
		MCMCRun(std::string,llu, llu , llu, std::string);
		void run();
		void run_serial();
		void write_state_to_file(int);
		void write_state_f_tens3(int);
		void normalize_f_tens3();
		void reset_f_tens3();
		FILE* pfile;
		void open_file_handle() { pfile = fopen(ofile.c_str(),"w");};
		void close_file_handle() { fclose(pfile);}; 
		int max_threads;
		std::vector<llu> get_thread_loads();
};

#endif
