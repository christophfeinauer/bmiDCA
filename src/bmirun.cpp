#include "bmirun.hpp"
#include "mcrun.hpp"
#include "hdfio.hpp"
#include "msastruct.hpp"

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

BmiRun::BmiRun(std::string fn, 
		std::string input_type,
		std::string output,
		std::string strgy,
		double epsilon, 
		llu nsamples, 
		llu eq_steps, 
		llu delta_steps, 
		double lambda,
		double alpha,
		double pc,
		double theta,
		llu max_runs,
		std::string f3tens_name, 
		std::string f2tens_name, 
		std::string charmap,
		char na_char,
		std::string storage_order) 
		: 
		output(output),
		strgy(strgy),
		epsilon(epsilon), 
		nsamples(nsamples), 
		eq_steps(eq_steps), 
		delta_steps(delta_steps),
		lambda(lambda),
		alpha(alpha),
		pc(pc),
		max_runs(max_runs){

	if (input_type=="hdf5") read_hdf5(fn, f3tens_name, f2tens_name, storage_order);
	if (input_type=="fasta") read_fasta(fn,theta, charmap, na_char, storage_order);

	q = (int) f2tens_obj_ptr->shape()[0];
	len = (int) f2tens_obj_ptr->shape()[1];
	
	if (q!=f3tens_obj_ptr->shape()[0] 
		|| q!=f3tens_obj_ptr->shape()[1] 
		|| (len*(len-1))/2 != f3tens_obj_ptr->shape()[2])
		throw ReadError("f3tens and f2tens dimensions non consistent");
	
	llu lenbn2 = (len*(len-1))/2;
	if (storage_order == "fortran"){
		f3tens_ptr.reset(new tens3(boost::extents[q][q][lenbn2],boost::fortran_storage_order()));
		f2tens_ptr.reset(new tens2(boost::extents[q][len],boost::fortran_storage_order()));
	}
	else{
		f3tens_ptr.reset(new tens3(boost::extents[q][q][lenbn2]));
		f2tens_ptr.reset(new tens2(boost::extents[q][len]));
	}

	llu gradient_nelements = f3tens_ptr->num_elements()+f2tens_ptr->num_elements();

	gradient.resize(gradient_nelements,0.0);	
}

void BmiRun::read_fasta(std::string& fn, double theta, std::string charmap, char na_char, std::string storage_order){

	MSAStruct msa =  MSAStruct(fn,theta,charmap,na_char);
	msa.get_frequencies(f3tens_obj_ptr, f2tens_obj_ptr, theta,storage_order);
	
}

void BmiRun::read_hdf5(std::string& fn, std::string& f3tens_name, std::string& f2tens_name, std::string storage_order){

	readtens3(fn,f3tens_obj_ptr,f3tens_name,storage_order);
	readtens2(fn,f2tens_obj_ptr,f2tens_name,storage_order);
}

void BmiRun::make_gradient(){
	
	double* f3_obj_data = f3tens_obj_ptr->data();
	double* f2_obj_data = f2tens_obj_ptr->data();
	double* f3_data = f3tens_ptr->data();
	double* f2_data = f2tens_ptr->data();
	
	llu k=0;	
	for (llu i=0; i<f2tens_ptr->num_elements(); ++i)
		gradient[k++] = alpha*gradient[k] + (1-alpha)*epsilon*(f2_obj_data[i] - f2_data[i]);
	for (llu i=0; i<f3tens_ptr->num_elements(); ++i)
		gradient[k++] = alpha*gradient[k] + (1-alpha)*epsilon*(f3_obj_data[i] - f3_data[i]); 

}

double BmiRun::get_squared_error(){

	double* f3_obj_data = f3tens_obj_ptr->data();
	double* f3_data = f3tens_ptr->data();
	double err=0.0;
	for (llu i=0; i<f3tens_obj_ptr->num_elements(); ++i)
			err += pow((f3_obj_data[i] - f3_data[i]),2);
	return err;
}

void BmiRun::add_pseudocount(){

	double* f3_flat = f3tens_obj_ptr->data();
	double* f2_flat = f2tens_obj_ptr->data();
	double dq = (double) q;
	double invq = 1/dq;
	double invq2 = invq*invq;
	
	for (llu i=0; i<f3tens_obj_ptr->num_elements();++i)
		f3_flat[i] = (1-pc)*f3_flat[i] + pc*invq2;
	for (llu i=0; i<f2tens_obj_ptr->num_elements();++i)
		f2_flat[i] = (1-pc)*f2_flat[i] + pc*invq;
}

void BmiRun::run(){
	
	if (pc!=0.0) add_pseudocount();
	if (strgy=="adaptive") run_adaptive();

}

// adpative run - makes gradient steps and checks if the squared error
// diminishes; if the squared error doesn't reach a new minimum for more than
// 10 steps, it either reduces epsilon or increases sample size
void BmiRun::run_adaptive(){

	MCMCRun mcrun(len,q,nsamples,eq_steps,delta_steps,"");
	mcrun.set_ftens_ptrs(f3tens_ptr,f2tens_ptr);
	set_ptrs_from_mcrun(mcrun);
	double err=0.0;
	double minerr = std::numeric_limits<double>::max();
	int c=0;
	bool epsilon_flag=true;
	bool samples_flag=false;
	for (llu r=0; r<max_runs; r++){
		mcrun.run("");
		mcrun.normalize();
		write_tensors(mcrun);
		make_gradient();
		if (lambda != 0.0) mcrun.env->efu->subtract_with_factor(gradient,2*lambda);
		mcrun.env->efu->add(gradient);
		double new_err = get_squared_error();
		if (new_err<minerr){
			minerr=new_err;	
			c=0;
		}
		else 
			c++;
		if (c>10){
			if (epsilon_flag){
				epsilon=0.5*epsilon;
				epsilon_flag=false;
				samples_flag=true;
				printf("NEW EPS: %lf\n",epsilon);
			}
			else {
				mcrun.nsamples= (llu) floor(1.5 * (double) mcrun.nsamples);
				epsilon_flag=true;
				samples_flag=false;
				printf("NEW SAMPLE SIZE: %llu\n",mcrun.nsamples);
			}
			minerr=new_err;
			c=0;
		}
			
		printf("ERR: %lf EPS: %lf\n",get_squared_error(),epsilon);
		mcrun.reset();
		r++;
	}
}

void BmiRun::set_ptrs_from_mcrun(MCMCRun& mcrun){
	env_ptr = mcrun.env;
	coup_ptr = mcrun.env->efu->get_coup_ptr();
	fields_ptr =  mcrun.env->efu->get_fields_ptr();
}

void BmiRun::write_tensors(MCMCRun& mcrun){
	writetens3(output,f3tens_ptr, "fij",mcrun.env->efu->storage_order);
	addtens2(output,f2tens_ptr,"fi",mcrun.env->efu->storage_order);
	addtens3(output,coup_ptr,"coup",mcrun.env->efu->storage_order);
	addtens2(output,fields_ptr,"fields",mcrun.env->efu->storage_order);
}



