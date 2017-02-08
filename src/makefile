CXX=g++
CPPFLAGS=-std=c++11

dca_sampler: dca_sampler.o state.o environ.o
	$(CXX) $(CPPFLAGS) -o dca_sampler dca_sampler.o state.o environ.o

dca_sampler.o: dca_sampler.cpp
	$(CXX) $(CPPFLAGS) -c dca_sampler.cpp
environ.o: environ.cpp state.o
	$(CXX) $(CPPFLAGS) -c environ.cpp
state.o: state.cpp
	$(CXX) $(CPPFLAGS) -c state.cpp
