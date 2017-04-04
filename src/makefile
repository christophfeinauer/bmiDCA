CXX=g++
CPPFLAGS= -L/usr/lib/x86_64-linux-gnu -std=c++11 -fopenmp -O3
CFLAGS=-I/usr/include/hdf5/serial 
LIB = -lhdf5_serial -lboost_program_options -lboost_filesystem -lboost_system


all: bmidca dca_sampler

bmidca: bmidca.o bmirun.o state.o env.o efu.o mcrun.o msastruct.o hdfio.o
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bmidca bmidca.o bmirun.o state.o env.o efu.o mcrun.o hdfio.o msastruct.o $(LIB)
dca_sampler: dca_sampler.o state.o env.o efu.o mcrun.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o dca_sampler dca_sampler.o state.o env.o efu.o mcrun.o hdfio.o $(LIB)
bmidca.o: bmidca.cpp env.hpp bmirun.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c bmidca.cpp 

bmirun.o: bmirun.cpp env.hpp mcrun.hpp hdfio.hpp msastruct.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c bmirun.cpp

dca_sampler.o: dca_sampler.cpp env.hpp state.hpp efu.hpp mcrun.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c dca_sampler.cpp

mcrun.o: mcrun.cpp mcrun.hpp env.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c mcrun.cpp
env.o: env.cpp env.hpp state.hpp efu.hpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c env.cpp
state.o: state.cpp state.hpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c state.cpp
efu.o: efu.cpp efu.hpp hdfio.hpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c efu.cpp 
hdfio.o: hdfio.cpp hdfio.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c hdfio.cpp
msastruct.o: msastruct.cpp msastruct.hpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c msastruct.cpp
