///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "hamiltonian.hpp"

struct SiteInfo2D
{
  //nn_1, nn_2, sub, bond-dep Hamiltonian, real_nn_1, real_nn_2 (real_sub = sub)
  vector<std::tuple<uint, uint, uint, Matrix3LD, int, int>> NearestNeighbours;
  Vector2LD Position;
};

class Honeycomb
{
  public:
    const uint NumSites;
    vector<Vector3LD> Cluster;
    vector<SiteInfo2D> ClusterInfo;
    const Vector3LD hField;
    const uint L1, L2, NumSublattices;
    long double ClusterEnergy;

    Honeycomb(const uint& hc_or_kek, const uint& type,
              const uint& num_sublattices, const uint& l1, const uint& l2,
              Hamiltonia& haminfo);
    void PrintLatticeParameters(std::ostream &out);
    void PrintHamiltonianParameters(std::ostream &out);
    void CalculateClusterOP();
    void PrintOP(std::ostream &out);

  private:
    Hamiltonia HamInfo;

    const uint HcOrKekule, ClusterType, NumUnitCells;

    vector<Vector2LD> SublatticeOffset;
    Vector2LD Translation1, Translation2;
    LATTICE_DIR;

    void CreateRhombicCluster();
    void CreateRectangularCluster();
    void CreateC3Cluster();
    void CreateKekuleCluster();

};

class Triangular
{
  public:
    const uint NumSites;
    vector<Vector3LD> Cluster;
    vector<SiteInfo2D> ClusterInfo;
    const Vector3LD hField;
    const uint L1, L2, NumSublattices;
    long double ClusterEnergy;

    Vector3LD ClusterFMOP,ClusterCombinedOP;
    Vector2LD ClusterStripyOP;
    Eigen::Matrix<long double, 3, 2> ClusterStripyOPMatrix;

    Triangular(const uint& l1, const uint& l2, const uint& num_sublattices,
                     const uint& num_defects,
                     Hamiltonia& haminfo,
                     const long double& defect_strength,
                     const long double& defect_lengthscale);
    void PrintLatticeParameters(std::ostream &out);
    void PrintHamiltonianParameters(std::ostream &out);
    void CalculateClusterOP();
    void PrintOP(std::ostream &out);

  private:
    Hamiltonia HamInfo;

    Vector2LD Translation1, Translation2;
    vector<Vector2LD> SublatticeOffset;
    LATTICE_DIR;

    //nn_1, nn_2, sub
    vector<std::tuple<uint,uint,uint>> Defects;
    const long double DefectStrength, DefectLengthScale;
    uint NumDefects;

    ArrayXLD StripySignsX, StripySignsY, StripySignsZ;

    void CreateRhombicCluster();

    void CreateDefectPositions();
    void AddDefectHamiltonia();

    void CreateStripySignMatrices();
    void SelectStripyOP();
    void AverageStripyOP();

};

#endif
