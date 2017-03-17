#ifndef BMIRUN_HPP
#define BMIRUN_HPP
#include "boost/multi_array.hpp" 
#include "env.hpp"
#include "mcrun.hpp"
typedef long long unsigned llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;

class BmiRun{
	double epsilon;
	int q;
	llu len;
	std::shared_ptr<Env> env_ptr;
	std::vector<double> gradient;
	public:
		BmiRun(std::string,std::string,std::string,std::string,double,llu,llu,llu,double,double,double,llu,std::string="fij",std::string="fi",std::string="fortran");
		void make_gradient();
		double get_squared_error();
		void add_pseudocount();
		void write_tensors(MCMCRun&);
		void set_ptrs_from_mcrun(MCMCRun&);
		llu lenbn2;
		// Tensors holding freqs for the objective function
		tens3_ptr f3tens_obj_ptr;
		tens2_ptr f2tens_obj_ptr;
		tens3_ptr f3tens_ptr;
		tens2_ptr f2tens_ptr;
		tens3_ptr coup_ptr;
		tens2_ptr fields_ptr;

		llu nsamples;
		llu eq_steps;
		llu delta_steps;
		llu max_runs;
		std::string output;
		std::string strgy;
		double lambda;
		double alpha;
		double pc;
		void run();
		void run_adaptive();
};
#endif
