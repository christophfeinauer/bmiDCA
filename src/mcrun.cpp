#include "mcrun.hpp"
#include <stdio.h>
#include <omp.h>

MCRun::MCRun(Env* env, llu nsamples) : env(env), nsamples(nsamples) {}
MCRun::MCRun(Env* env, llu nsamples, tens3 * corr_tens3) : env(env), nsamples(nsamples), f_tens3(f_tens3) {}
MCMCRun::MCMCRun(Env * env, llu nsamples, llu eq_steps, llu delta_steps, std::string ofile) : MCRun(env,nsamples), eq_steps(eq_steps), delta_steps(delta_steps), ofile(ofile), max_threads(omp_get_max_threads()) {}


void MCMCRun::run(){
	open_file_handle();
	printf("Number of threads: %d\n",omp_get_max_threads());
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
				if (f_tens3!=NULL)
					write_state_f_tens3(state_nr);
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

void MCMCRun::write_state_f_tens3(int i){
	State& state = env->state_vec[i];

	int len = state.len;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int j=i+1; j<len; ++j)
			(*f_tens3)[state.seq[i]][state.seq[j]][l++]+=1.0f;
	f_tens3_count++;
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
