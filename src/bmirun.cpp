#include "bmirun.hpp"
#include "mcrun.hpp"
#include "hdfio.hpp"

struct ReadError : public std::runtime_error{
	ReadError(std::string const & message)
		: std::runtime_error(message)
	{}
};

Simple::Simple(std::string fn, 
		double epsilon, 
		llu nsamples, 
		llu eq_steps, 
		llu delta_steps, 
		double lambda,
		double alpha,
		llu max_runs,
		std::string f3tens_name, 
		std::string f2tens_name, 
		std::string storage_order) 
		: 
		epsilon(epsilon), 
		nsamples(nsamples), 
		eq_steps(eq_steps), 
		delta_steps(delta_steps),
		lambda(lambda),
		alpha(alpha),
		max_runs(max_runs) {

	readtens3(fn,f3tens_obj_ptr,f3tens_name,storage_order);
	readtens2(fn,f2tens_obj_ptr,f2tens_name,storage_order);
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
		f2tens_ptr.reset(new tens2(boost::extents[q][len]));}
	llu gradient_nelements = f3tens_ptr->num_elements()+f2tens_ptr->num_elements();
	gradient.resize(gradient_nelements,0.0);	
	
}

void Simple::make_gradient(){
	
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

double Simple::get_squared_error(){
	double* f3_obj_data = f3tens_obj_ptr->data();
	double* f3_data = f3tens_ptr->data();
	double err=0.0;
	for (llu i=0; i<f3tens_obj_ptr->num_elements(); ++i)
			err += pow((f3_obj_data[i] - f3_data[i]),2);
	return err;
}

void Simple::run(){

	MCMCRun mcrun(len,q,nsamples,eq_steps,delta_steps,"test");
	env_ptr = mcrun.env;
	mcrun.set_ftens_ptrs(f3tens_ptr,f2tens_ptr);

	llu r=0;
	double err=0.0;
	double minerr = 0.0;
	int c=0;
	bool epsilon_flag=true;
	bool samples_flag=false;
	while(r<max_runs){
		std::stringstream ofile;
		ofile << "test_" << r;
		mcrun.run(ofile.str().c_str());
		mcrun.normalize();
		make_gradient();
		mcrun.env->efu->subtract_with_factor(gradient,2*lambda);
		writetens3("test.hdf",f3tens_ptr, ofile.str(),mcrun.env->efu->storage_order);
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
				printf("new epsilon: %lf\n",epsilon);
			}
			else {
				mcrun.nsamples=2*mcrun.nsamples;
				epsilon_flag=true;
				samples_flag=false;
				printf("new nsamples: %llu\n",mcrun.nsamples);
			}
			minerr=new_err;
			c=0;
		}
			
		printf("ERR: %lf\n",get_squared_error());
		mcrun.reset();
		r++;
	}
}



