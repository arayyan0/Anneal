#file   Makefile
#author Ahmed Rayyan
#date   December 2, 2019
#brief  classical Monte Carlo make file
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -c
SIM = multipole_PT
LATTICE = triangular

all: sim

sim: sim_${SIM}.o ${LATTICE}.o hamiltonian.o spin.o common.o
	${CXX} sim_${SIM}.o ${LATTICE}.o hamiltonian.o spin.o common.o -o sim

sim_${SIM}.o: sim_${SIM}.cpp ${LATTICE}.hpp hamiltonian.hpp
	${CXX} ${CXXFLAGS} sim_${SIM}.cpp

${LATTICE}.o: ${LATTICE}.cpp ${LATTICE}.hpp hamiltonian.hpp common.hpp
	${CXX} ${CXXFLAGS} ${LATTICE}.cpp

hamiltonian.o: hamiltonian.cpp hamiltonian.hpp spin.hpp common.hpp
	${CXX} ${CXXFLAGS} hamiltonian.cpp

spin.o: spin.cpp spin.hpp common.hpp
	${CXX} ${CXXFLAGS} spin.cpp

common.o: common.cpp common.hpp
	${CXX} ${CXXFLAGS} common.cpp

clean:
	rm *.o sim
