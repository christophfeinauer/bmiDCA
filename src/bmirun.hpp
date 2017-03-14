#ifndef BMIRUN_HPP
#define BMIRUN_HPP
#include "boost/multi_array.hpp" 
#include "env.hpp"
typedef long long unsigned llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;

class Simple{
	double epsilon;
	int q;
	llu len;
	std::shared_ptr<Env> env_ptr;
	std::vector<double> gradient;
	public:
		Simple(std::string,double,llu,llu,llu,double,double,llu,std::string="fij",std::string="fi",std::string="fortran");
		void make_gradient();
		double get_squared_error();
		llu lenbn2;
		// Tensors holding freqs for the objective function
		tens3_ptr f3tens_obj_ptr;
		tens2_ptr f2tens_obj_ptr;
		tens3_ptr f3tens_ptr;
		tens2_ptr f2tens_ptr;
		llu nsamples;
		llu eq_steps;
		llu delta_steps;
		llu max_runs;
		double lambda;
		double alpha;

		void run();
};
#endif
