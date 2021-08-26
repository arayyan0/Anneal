#file   Makefile
#author Ahmed Rayyan
#date   December 2, 2019
#brief  classical Monte Carlo make file
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -c
SIM = multipole

all: sim

sim: sim_${SIM}.o MC.o lattice.o hamiltonian.o common.o
	${CXX} sim_${SIM}.o MC.o lattice.o hamiltonian.o common.o -o sim

sim_${SIM}.o: sim_${SIM}.cpp MC.hpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} sim_${SIM}.cpp

MC.o: MC.cpp MC.hpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} MC.cpp
	
lattice.o: lattice.cpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} lattice.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} hamiltonian.cpp

common.o: common.cpp common.hpp
	${CXX} ${CXXFLAGS} common.cpp

clean:
	rm *.o sim
