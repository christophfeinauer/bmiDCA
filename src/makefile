CXX=g++
CPPFLAGS= -L/usr/lib/x86_64-linux-gnu -std=c++11 -g
CFLAGS=-I/usr/include/hdf5/serial 
LIB = -lhdf5_serial

dca_sampler: dca_sampler.o state.o env.o efu.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o dca_sampler dca_sampler.o state.o env.o efu.o $(LIB)

dca_sampler.o: dca_sampler.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c dca_sampler.cpp
env.o: env.cpp state.o efu.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c env.cpp
state.o: state.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c state.cpp
efu.o: efu.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c efu.cpp
