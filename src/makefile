CXX=g++
CPPFLAGS= -L/usr/lib/x86_64-linux-gnu -std=c++11 -fopenmp -g -O3
CFLAGS=-I/usr/include/hdf5/serial 
LIB = -lhdf5_serial -lboost_program_options -lboost_filesystem -lboost_system


all: bmidca dca_sampler

bmidca: bmidca.o bmirun.o state.o env.o efu.o mcrun.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bmidca bmidca.o bmirun.o state.o env.o efu.o mcrun.o hdfio.o $(LIB)
dca_sampler: dca_sampler.o state.o env.o efu.o mcrun.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o dca_sampler dca_sampler.o state.o env.o efu.o mcrun.o hdfio.o $(LIB)
bmidca.o: bmidca.cpp mcrun.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c bmidca.cpp 

bmirun.o: bmirun.cpp mcrun.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c bmirun.cpp

dca_sampler.o: dca_sampler.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c dca_sampler.cpp

mcrun.o: mcrun.cpp env.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c mcrun.cpp
env.o: env.cpp state.o efu.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c env.cpp
state.o: state.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c state.cpp
efu.o: efu.cpp hdfio.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c efu.cpp 
hdfio.o: hdfio.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c hdfio.cpp
