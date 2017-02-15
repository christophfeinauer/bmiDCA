CXX=g++
CPPFLAGS=-std=c++11 -g

dca_sampler: dca_sampler.o state.o env.o efu.o
	$(CXX) $(CPPFLAGS) -o dca_sampler dca_sampler.o state.o env.o efu.o

dca_sampler.o: dca_sampler.cpp
	$(CXX) $(CPPFLAGS) -c dca_sampler.cpp
env.o: env.cpp state.o efu.o
	$(CXX) $(CPPFLAGS) -c env.cpp
state.o: state.cpp
	$(CXX) $(CPPFLAGS) -c state.cpp
efu.o: efu.cpp
	$(CXX) $(CPPFLAGS) -c efu.cpp
