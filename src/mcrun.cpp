#include "mcrun.hpp"
#include <stdio.h>
#include <omp.h>

struct AccessError : public std::runtime_error{
	AccessError(std::string const & message)
		: std::runtime_error(message)
	{}
};

void MCMCRun::report_threads() { printf("Number of threads: %d\n",max_threads); }

MCRun::MCRun(std::shared_ptr<Env> env, 
		llu nsamples, 
		tens3_ptr f3tens_ptr) 
		: 
		env(env), 
		nsamples(nsamples), 
		f3tens_ptr(f3tens_ptr) { }

MCMCRun::MCMCRun(std::shared_ptr<Env> env, 
		llu nsamples, 
		llu eq_steps, 
		llu delta_steps, 
		std::string ofile) 
		: 
		MCRun(env,nsamples), 
		eq_steps(eq_steps), 
		delta_steps(delta_steps), 
		ofile(ofile), 
		max_threads(omp_get_max_threads()) { report_threads(); }

MCMCRun::MCMCRun(std::string infile, 
		llu nsamples, 
		llu eq_steps, 
		llu delta_steps, 
		std::string ofile)
		: 
		MCRun(nullptr,nsamples), 
		eq_steps(eq_steps), 
		delta_steps(delta_steps), 
		ofile(ofile), 
		max_threads(omp_get_max_threads()){

	std::shared_ptr<Pairwise> pw(new Pairwise(infile,"J","h"));	
	env.reset(new Env(pw,max_threads,pw->len,pw->q,rand()));
	report_threads();
}

MCMCRun::MCMCRun(int len, 
		int q, 
		llu nsamples, 
		llu eq_steps, 
		llu delta_steps, 
		std::string ofile,
		std::string storage_order)
		:
	 	MCRun(nullptr,nsamples), 
		eq_steps(eq_steps), 
		delta_steps(delta_steps), 
		ofile(ofile), 
		max_threads(omp_get_max_threads()){

	std::shared_ptr<Pairwise> pw(new Pairwise(len,q));	
	env.reset(new Env(pw,max_threads,pw->len,pw->q,rand()));
	report_threads();
}

void MCMCRun::run(std::string _ofile){
	ofile = _ofile;
	run();
}

void MCMCRun::run(){
	open_file_handle();
	#pragma omp parallel for
	for (int state_nr = 0; state_nr<env->nr_states; ++state_nr)
		for (llu eq_step = 0; eq_step < eq_steps; ++eq_step)
			env->step(state_nr);
	
	std::vector<llu> thread_loads = get_thread_loads();
	llu samples_done=0;
	#pragma omp parallel for
	for (int state_nr = 0; state_nr<env->nr_states; ++state_nr)
		for (llu nsample = 0; nsample < thread_loads[state_nr]; ++nsample){
			env->reset_moves(state_nr);
			for (llu mc_step = 0; mc_step < delta_steps; ++mc_step)
				env->step(state_nr);
			#pragma omp critical
			{
				write_state_to_file(state_nr);
				if (f3tens_ptr!=NULL)
					write_state_f3tens(state_nr);
				if (f2tens_ptr!=NULL)
					write_state_f2tens(state_nr);
				samples_done++;
				printf("\rsample %llu/%llu",samples_done,nsamples);
				fflush(stdout);
			}
	}
	printf("\n");
	close_file_handle();
}


std::vector<llu> MCMCRun::get_thread_loads(){
	std::vector<llu> thread_loads(env->nr_states,floor((long double) nsamples / (long double) env->nr_states));
	for (int i=0; i<nsamples % env->nr_states; ++i)
		thread_loads[i]++;
	return thread_loads;
	
}

void MCMCRun::write_state_f3tens(int i){
	if (f3tens_ptr==NULL)
		throw AccessError("Trying to write to f3tens, but f3tens_ptr not set");
	State& state = env->state_vec[i];
	int len = state.len;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int j=i+1; j<len; ++j)
			(*f3tens_ptr)[state.seq[i]][state.seq[j]][l++]+=1.0f;
	f3tens_count++;
}

void MCMCRun::write_state_f2tens(int i){
	if (f2tens_ptr==NULL)
		throw AccessError("Trying to write to f2tens, but f2tens_ptr not set");
	State& state = env->state_vec[i];
	int len = state.len;
	for (int i=0; i<len; ++i)
		(*f2tens_ptr)[state.seq[i]][i]+=1.0f;
	f2tens_count++;
}

void MCMCRun::normalize_f3tens(){
	if (f3tens_ptr==NULL)
		throw AccessError("Trying to normalize ftens3, but ftens3 not set");
	State& state = env->state_vec[0];
	int len = state.len;
	int q = state.q;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int j=i+1; j<len; ++j){
			for (int a=0; a<q; ++a)
				for (int b=0; b<q; ++b)
					(*f3tens_ptr)[a][b][l]/=f3tens_count;
			l++;
		}

}

void MCMCRun::normalize_f2tens(){
	if (f2tens_ptr==NULL)
		throw AccessError("Trying to normalize ftens2, but ftens2_ptr not set");
	State& state = env->state_vec[0];
	int len = state.len;
	int q = state.q;
	for (int i=0; i<len; ++i)
		for (int a=0; a<q; ++a)
			(*f2tens_ptr)[a][i]/=f2tens_count;

}

void MCMCRun::reset_f3tens(){
	if (f3tens_ptr==NULL)
		throw AccessError("Trying to reset f3tens, but f3tens not set");
	State& state = env->state_vec[0];
	int len = state.len;
	int q = state.q;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int j=i+1; j<len; ++j){
			for (int a=0; a<q; ++a)
				for (int b=0; b<q; ++b)
					(*f3tens_ptr)[a][b][l]=0.0;
			l++;
		}
	f3tens_count = 0;
}

void MCMCRun::reset_f2tens(){
	if (f2tens_ptr==NULL)
		throw AccessError("Trying to reset f3tens, but f3tens not set");
	State& state = env->state_vec[0];
	int len = state.len;
	int q = state.q;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int a=0; a<q; ++a)
			(*f2tens_ptr)[a][i]=0.0;
	f2tens_count = 0;
}

void MCMCRun::write_state_to_file(int i){
	State& state = env->state_vec[i];
	double en = env->get_energy();
	double acc = state.acc();
	fprintf(pfile,">thread_%d_step_%llu_en_%lf_acc_%lf\n%d",i,state.moves_proposed_total,en,acc,state.seq[0]);
	for (int i=1; i<state.len; ++i){
		fprintf(pfile," %d",state.seq[i]);
	}
	fprintf(pfile,"\n");
}
