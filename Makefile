#file   Makefile
#author Ahmed Rayyan
#date   December 2, 2019
#brief  classical Monte Carlo make file
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -c
#DEBUGFLAGS = -pg -g
SIM = multipole
<<<<<<< HEAD

LATTICE = lattice
=======
>>>>>>> 569b9a13e5a38ae3866330adf20f7d4a8729bcdc

all: sim

sim: sim_${SIM}.o MC.o lattice.o hamiltonian.o common.o
	${CXX} sim_${SIM}.o MC.o lattice.o hamiltonian.o common.o -o sim

sim_${SIM}.o: sim_${SIM}.cpp MC.hpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} sim_${SIM}.cpp

MC.o: MC.cpp MC.hpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} MC.cpp

lattice.o: lattice.cpp lattice.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} lattice.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} hamiltonian.cpp

common.o: common.cpp common.hpp
	${CXX} ${CXXFLAGS} ${DEBUGFLAGS} common.cpp

clean:
	rm *.o sim
