///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the types of Lattices
#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "hamiltonian.hpp"

struct SiteInfo2D
{
  //nn_1, nn_2, sub, bond-dep Hamiltonian, real_nn_1, real_nn_2 (real_sub = sub)
  vector<std::tuple<uint, uint, uint, Matrix3LD, int, int, uint>> NearestNeighbours;
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
    vector<Vector3LD> ClusterSSf;

    long double EBar, E2Bar;
    vector<Vector3LD> SSfBar;

    vector<Vector2LD> SSFPoints;

    Honeycomb(const uint& hc_or_kek, const uint& type,
              const uint& num_sublattices, const uint& l1, const uint& l2,
              Hamiltonia& haminfo);
    void PrintLatticeParameters(std::ostream &out);
    void PrintHamiltonianParameters(std::ostream &out);
    void CalculateClusterOP();
    void PrintOP(std::ostream &out);
    void PrintThermalObservables(std::ostream &out);

  private:
    Hamiltonia HamInfo;

    const uint HcOrKekule, ClusterType, NumUnitCells;

    vector<Vector2LD> SublatticeOffset;
    Vector2LD Translation1, Translation2;
    LATTICE_DIR;
    RECIP_DIR;


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
    vector<Vector3LD> ClusterSSf;

    long double EBar, E2Bar;
    vector<Vector3LD> SSfBar;

    vector<Vector2LD> SSFPoints;

    Triangular(const uint& l1, const uint& l2, const uint& num_sublattices,
                     const uint& num_defects,
                     Hamiltonia& haminfo,
                     const long double& defect_quad,
                     const long double& defect_octo,
                     const long double& defect_lengthscale);
    void PrintLatticeParameters(std::ostream &out);
    void PrintHamiltonianParameters(std::ostream &out);
    void CalculateClusterOP();
    void PrintOP(std::ostream &out);
    void PrintThermalObservables(std::ostream &out);

  private:
    Hamiltonia HamInfo;

    Vector2LD Translation1, Translation2;
    vector<Vector2LD> SublatticeOffset;
    LATTICE_DIR;
    RECIP_DIR;
    SYM_DIR;

    vector<std::string> SSFPointsLabel;

    //nn_1, nn_2, sub
    vector<std::tuple<uint,uint,uint>> Defects;
    const long double DefectQuad, DefectOcto, DefectLengthScale;
    uint NumDefects;

    void CreateRhombicCluster();

    void CreateDefectPositions();
    void AddDefectHamiltonia();

    void DebugHamiltonians();

};

#endif
