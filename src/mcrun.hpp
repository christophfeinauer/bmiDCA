#ifndef MCRUN_H
#define MCRUN_H
#include "env.hpp"
typedef long long unsigned llu;
typedef boost::multi_array<double,3> tens3;
typedef boost::multi_array<double,2> tens2;
typedef std::shared_ptr<tens3> tens3_ptr;
typedef std::shared_ptr<tens2> tens2_ptr;

class MCRun{
	public:
		std::shared_ptr<Env> env;
		llu nsamples;
		void virtual run() = 0;
		MCRun(std::shared_ptr<Env>,llu,tens3_ptr=nullptr);
		tens3_ptr f3tens_ptr = nullptr;
		tens2_ptr f2tens_ptr = nullptr;
		llu f3tens_count = 0;
		llu f2tens_count = 0;
};

class MCMCRun : public MCRun {
	public: 
		llu eq_steps;
		llu delta_steps;	
		std::string ofile;
		MCMCRun(std::shared_ptr<Env>,llu, llu , llu, std::string);
		MCMCRun(std::string,llu, llu , llu, std::string);
		MCMCRun(int,int,llu, llu , llu, std::string, std::string="fortran");
		void set_f3tens_ptr(tens3_ptr _f3) {f3tens_ptr=_f3;};
		void set_f2tens_ptr(tens2_ptr _f2) {f2tens_ptr=_f2;};
		void set_ftens_ptrs(tens3_ptr _f3, tens2_ptr _f2) {set_f3tens_ptr(_f3); set_f2tens_ptr(_f2);};
		void run();
		void run(std::string);
		void run_serial();
		void write_state_to_file(int);
		void write_state_f3tens(int);
		void normalize_f3tens();
		void reset_f3tens();
		void write_state_f2tens(int);
		void normalize_f2tens();
		void normalize(){normalize_f3tens(); normalize_f2tens();};
		void reset_f2tens();
		void reset() {reset_f3tens(); reset_f2tens(); env->reset_all_moves();};
		FILE* pfile;
		void open_file_handle() { pfile = (ofile == "") ? nullptr : fopen(ofile.c_str(),"w");};
		void close_file_handle() { if (pfile != nullptr) fclose(pfile);}; 
		int max_threads;
		void report_threads();
		std::vector<llu> get_thread_loads();
};

#endif
