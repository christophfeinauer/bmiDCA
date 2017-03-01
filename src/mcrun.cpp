#include "mcrun.hpp"
#include <stdio.h>

MCRun::MCRun(Env* env, llu nsamples) : env(env), nsamples(nsamples) {}
MCRun::MCRun(Env* env, llu nsamples, tens3 * corr_tens3) : env(env), nsamples(nsamples), f_tens3(f_tens3) {}

MCMCRun::MCMCRun(Env * env, llu nsamples, llu eq_steps, llu delta_steps, std::string ofile) : MCRun(env,nsamples), eq_steps(eq_steps), delta_steps(delta_steps), ofile(ofile) {}

void MCMCRun::run(){
	for (llu eq_step = 0; eq_step < eq_steps; ++eq_step)
		env->step();
	for (llu nsample = 0; nsample < nsamples; ++nsample){
		printf("\rsample %llu/%llu",nsample+1,nsamples);
		fflush(stdout);
		env->state->reset_moves();
		for (llu mc_step = 0; mc_step < delta_steps; ++mc_step)
			env->step();
		write_state_to_file();
		if (f_tens3!=NULL)
			write_state_f_tens3();
	}
	printf("\n");
			
}

void MCMCRun::write_state_f_tens3(){
	int len = env->state->len;
	llu l=0;
	for (int i=0; i<len; ++i)
		for (int j=i+1; j<len; ++j)
			(*f_tens3)[env->state->seq[i]][env->state->seq[j]][l++]+=1.0f;
}

void MCMCRun::write_state_to_file(){
	FILE * pfile;
	pfile = fopen(ofile.c_str(),"a");
	double en = env->get_energy();
	double acc = env->state->acc();
	fprintf(pfile,">step_%llu_en_%lf_acc_%lf\n%d ",env->state->moves_proposed_total,en,acc,env->state->seq[0]);
	for (int i=1; i<env->state->len; ++i){
		fprintf(pfile," %d",env->state->seq[i]);
	}
	fprintf(pfile,"\n");
	fclose(pfile);
}
