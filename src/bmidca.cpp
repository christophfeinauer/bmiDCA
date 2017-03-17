#include "bmirun.hpp"
#include "env.hpp"
#include <iostream>
#include "boost/program_options.hpp" 
namespace po = boost::program_options;

int main(int argc, char** argv){

	std::string strgy, input, input_type, output, storage_order, f3tens_name, f2tens_name;
	llu nsamples, eq_steps, delta_steps,  max_runs;
	double lambda, epsilon, alpha, pc;

	po::options_description desc("Options");
	desc.add_options()
		("help","Show help")
		("input",po::value<std::string>(&input),"Input file (fasta or hdf5)")
		("input_type",po::value<std::string>(&input_type)->default_value(std::string("hdf5")),"Input file type ('fasta' or 'hdf5')")
		("output",po::value<std::string>(&output),"Output file (HDF5 format, snapshots of parameters and frequencies are saved after every gradient step")
		("storage_order",po::value<std::string>(&storage_order)->default_value("fortran"),"Storage order for tensor if input is HDF5 (tensor should be q x q x N(N-1)/2) in fortran order). Default value should work with Julia/Matlab output") 
		("f3tens_name",po::value<std::string>(&f3tens_name)->default_value("fij"),"dataset name of two-site frequency tensor if input is hdf5")
		("f2tens_name",po::value<std::string>(&f2tens_name)->default_value("fi"),"dataset name of one-site frequency tensor if input is hdf5")
		("strategy",po::value<std::string>(&strgy)->default_value(std::string("adaptive")),"minimization strategy - only 'adaptive' for now")
		("nsamples",po::value<llu>(&nsamples)->default_value(500),"Initial number of MC samples used in gradient calculation")
		("eq_steps",po::value<llu>(&eq_steps)->default_value(100000),"Number of equilibration MC-steps per chain")
		("delta_steps",po::value<llu>(&delta_steps)->default_value(10000),"Number of steps between MC-samples")
		("max_runs",po::value<llu>(&max_runs)->default_value(1000),"Maximum number of gradient descent steps")
		("epsilon",po::value<double>(&epsilon)->default_value(0.1),"Initial gradient descent step size")
		("pc",po::value<double>(&pc)->default_value(0.001),"Pseudocount to be added to objective function")
		("lambda",po::value<double>(&lambda)->default_value(0.0),"L2 regularization constant");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    std::cout << desc << "\n";
	    return 1;
	}

	if (vm.count("output") != 1 || vm.count("input")!=1) {
	    std::cout << "Please provide an input and output file" << "\n";
	    return 1;
	}

	BmiRun bmirun(input,input_type,output,strgy,epsilon,nsamples,eq_steps,delta_steps,lambda,alpha,pc,max_runs, f3tens_name, f2tens_name);
	bmirun.run();
	return 0;
}
