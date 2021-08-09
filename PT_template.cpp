#include <chrono>
#include "Eigen/Dense"
#include <iostream>
#include <mpi.h>
#include <random>
#include <vector>

int main(int argc, char *argv[]){
  int rank, size, tag=1;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Datatype MPI_VECTOR3LD;
  MPI_Type_contiguous(3, MPI_LONG_DOUBLE, &MPI_VECTOR3LD);
  MPI_Type_commit(&MPI_VECTOR3LD);

  // data buffers
  long double e, other_e;

  typedef Eigen::Matrix<long double, 3, 1> Vector3LD;
  std::vector<Vector3LD> vv, other_vv;

  //filling data in each process
  e=rank + 0.5;
  Vector3LD v1, v2;
  v1 << 1+6*rank, 2+6*rank, 3+6*rank;
  v2 << 4+6*rank, 5+6*rank, 6+6*rank;
  vv.push_back(v1);
  vv.push_back(v2);

  //figure out other process to exchange with
  int other_rank = rank%2==0? rank+1 : rank-1;// partners are  0<->1, 2<->3, ... size-2<->size-1
  // int other_rank = rank%2==0? rank-1 : rank+1;// partners are -1<->0, 1<->2, ... size-1<->size
  // std::cout << rank << " " << other_rank << std::endl;

  //if statement excludes the edge cases in the second statement
  if ((other_rank >= 0) and (other_rank < size)){
    MPI_Sendrecv(      &e, 1, MPI_LONG_DOUBLE, other_rank, tag,
                 &other_e, 1, MPI_LONG_DOUBLE, other_rank, tag,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //figure out whether to swap based on random number generation
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 RNG(seed);
    std::bernoulli_distribution boo(0.5);
    bool acceptswap = false;
    if (rank%2 == 0){
       acceptswap = boo(RNG);
      MPI_Ssend(&acceptswap, 1, MPI_C_BOOL, other_rank, tag, MPI_COMM_WORLD);
    } else {
      MPI_Recv( &acceptswap, 1, MPI_C_BOOL, other_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //exchange e and vector
    if (acceptswap){
      e = other_e;

      other_vv.resize(vv.size());
      MPI_Sendrecv(       vv.data(),       vv.size(), MPI_VECTOR3LD, other_rank, tag,
                    other_vv.data(), other_vv.size(), MPI_VECTOR3LD, other_rank, tag,
                                                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      vv = other_vv;
    }

    //output resutls
    for (int j=0;j<size;j++){
      if (rank == j){
        std::cout << rank << " " << acceptswap << " " << vv[0].transpose() << " " << vv[1].transpose() << std::endl;
      }
    }
  }

  MPI_Finalize();
  return 0;
}
