#file   Makefile
#author Ahmed Rayyan
#date   December 2, 2019
#brief  classical Monte Carlo make file
CXX = g++
CXXFLAGS = -std=c++17 -O3 -c
SIM = sim_kga

all: sim

sim: ${SIM}.o lattice.o hamiltonian.o spin.o common.o
	${CXX} ${SIM}.o lattice.o hamiltonian.o spin.o common.o -o sim

${SIM}.o: ${SIM}.cpp lattice.hpp hamiltonian.hpp
	${CXX} ${CXXFLAGS} ${SIM}.cpp

lattice.o: lattice.cpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} lattice.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp spin.hpp common.hpp
	${CXX} ${CXXFLAGS} hamiltonian.cpp

spin.o: spin.cpp spin.hpp common.hpp
	${CXX} ${CXXFLAGS} spin.cpp

common.o: common.cpp common.hpp
	${CXX} ${CXXFLAGS} common.cpp

clean:
	rm *.o sim
