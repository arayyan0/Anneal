#file   Makefile
#author Ahmed Rayyan
#date   December 2, 2019
#brief  classical Monte Carlo make file
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -c
#DEBUGFLAGS = -pg -g
SIM = multipole

LATTICE = lattice

all: sim

sim: sim_${SIM}.o MC.o ${LATTICE}.o hamiltonian.o common.o
	${CXX} sim_${SIM}.o MC.o ${LATTICE}.o hamiltonian.o common.o -o sim

sim_${SIM}.o: sim_${SIM}.cpp MC.hpp ${LATTICE}.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} sim_${SIM}.cpp

MC.o: MC.cpp MC.hpp ${LATTICE}.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} MC.cpp

${LATTICE}.o: ${LATTICE}.cpp ${LATTICE}.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} ${LATTICE}.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} hamiltonian.cpp

common.o: common.cpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} common.cpp

clean:
	rm *.o sim
