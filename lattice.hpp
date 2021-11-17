///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the types of Lattices
#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "hamiltonian.hpp"

class Reciprocal
{
  public:
    Reciprocal(Matrix3LD& primitive_cell, const uint& which_lattice);
    Matrix3XLD QPoints;
    vector<std::string> QPointsLabel;
  private:
    Matrix3LD CrystalReciprocalVectors;
    void CreateReciprocalVectors(const Matrix3LD& translation_vectors,
                                       Matrix3LD& recip_vectors);
    void CreateSymmetryPoints(const uint& which_lattice);
};

struct SiteInfo
{
  //Hamiltonian, type, index of NN, which NNvector
  vector<std::tuple<Matrix3LD, uint, uint, uint>> NNInfo;
};

class Lattice
{
  public:
    uint NumSites;
    Matrix3XLD Cluster;
    long double ClusterEnergy;
    Vector3LD hField;
    vector<SiteInfo> NNSiteInfo;
    const uint L1, L2, L3;
    uint NumTotalSublattices;
    long double EBar, E2Bar;
    // constructor
    Lattice(const uint& which_lattice,
               const uint& shape,
               const uint& l1, const uint& l2, const uint& l3,
               Hamiltonia& haminfo
              );
    void PrintLatticeParameters(std::ostream &out);
    void PrintHamiltonianParameters(std::ostream &out);
    void PrintThermalObservables(std::ostream &out);
    void CalculateOP();
    void PrintOP(std::ostream &out);
    void CreateDefectPositions();
    void AddDefectHamiltonia(const long double& defect_quad,
                             const long double& defect_octo,
                             const long double& defect_lengthscale);
    void CheckHamiltonians();
  private:
    const uint WhichLattice;
    const uint NumUnitCells;
    Hamiltonia HamInfo;
    uint LatticeDimensions;
    uint CoordinationNumber;
    uint NumSublattices;
    vector<int> NNBondType;
    Matrix3LD PrimitiveCell;
    Matrix3LD ConventionalCell;
    Matrix3LD LocalOctahedra;
    Matrix3LD SpinBasis;
    Matrix3XLD SublatticeVectors;
    Matrix3XLD NNVectors;
    Matrix3XLD Positions;
    MatrixXXLDc Phases;
    Matrix3XLD QPoints;
    vector<std::string> QPointsLabel;
    Matrix3XLD SSFComponents;

    vector<uint> Defects;
    long double DefectQuad       ;
    long double DefectOcto       ;
    long double DefectLengthScale;

    ABC_DIR;

    void SpecifyLattice(const uint& shape);
    void CreateLatticePositions();
    void FindNearestNeighbours();
    void CalculatePhases();

    void CheckPositions();
    void CheckNNBonds();

};

#endif
