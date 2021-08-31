///  @file     lattice.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the types of Lattices
#include "lattice.hpp"

Honeycomb::Honeycomb(const uint& hc_or_kek, const uint& type,
                     const uint& num_sublattices, const uint& l1, const uint& l2,
                     Hamiltonia& haminfo
                   ):
 HcOrKekule(hc_or_kek), ClusterType(type), NumSublattices(num_sublattices), L1(l1), L2(l2),
 NumUnitCells(l1*l2), NumSites(l1*l2*num_sublattices),
 HamInfo(haminfo), hField(haminfo.hField)
{
  SublatticeOffset.resize(NumSublattices);

  Cluster.resize(NumSites);
  ClusterInfo.resize(NumSites);

  if (HcOrKekule == 0){
    if      (NumSublattices == 2) CreateRhombicCluster();
    else if (NumSublattices == 4) CreateRectangularCluster();
    else if (NumSublattices == 6) CreateC3Cluster();
  }
  else if (HcOrKekule == 1){
    if (NumSublattices == 6) CreateKekuleCluster();
  }

}

void Honeycomb::CreateRhombicCluster()
{
  Matrix3LD h1, h2, h3;
  if (ClusterType == 1){ //rhombus oriented around z-bond
    Translation1 = a1;
    Translation2 = a2;
    h1 = HamInfo.Hx;
    h2 = HamInfo.Hy;
    h3 = HamInfo.Hz;
  } else if (ClusterType == 2){ //rhombus oriented around x-bond
    Translation1 = a1-a2;
    Translation2 = a1;
    h1 = HamInfo.Hy;
    h2 = HamInfo.Hz;
    h3 = HamInfo.Hx;
  } else if (ClusterType == 3){ //rhombus oriented around y-bond
    Translation1 = a2;
    Translation2 = a2-a1;
    h1 = HamInfo.Hz;
    h2 = HamInfo.Hx;
    h3 = HamInfo.Hy;
  }

  SublatticeOffset[0] = (1.0*Translation1+1.0*Translation2)/3.0;
  SublatticeOffset[1] = (2.0*Translation1+2.0*Translation2)/3.0;

  Vector3LD some_vec;
  int higher_x, higher_y, lower_x, lower_y;
  uint flat_index;
  Vector2LD position;
  for (int y=0; y<L2; ++y){
    for (int x=0; x<L1; ++x){ //both x,y should NOT be ints as they may become negative
      for (uint sub =0; sub<NumSublattices; ++sub){
        // cout << "unit cell: " << x << " " << y << endl;
        higher_x = (x+1)%L1;
        higher_y = (y+1)%L2;
        //if index becomes -1, wrap it around to l-1
        lower_x = (x-1);
        lower_y = (y-1);
        if (lower_x >= 0){lower_x = lower_x%L1;}
        else if (lower_x < 0){lower_x = L1-1;}
        if (lower_y >= 0){lower_y = lower_y%L2;}
        else if (lower_y < 0){lower_y = L2-1;}
        // cout << lower_x << " " << lower_y << endl;
        position = x*Translation1+y*Translation2+SublatticeOffset[sub];
        BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);

        Cluster[flat_index] = some_vec;
        if (flat_index%NumSublattices==0){
          ClusterInfo[flat_index+0] = {
                                      {std::make_tuple(lower_x,       y, 1, h1, x-1,   y),
                                       std::make_tuple(      x, lower_y, 1, h2,   x, y-1),
                                       std::make_tuple(      x,       y, 1, h3,   x,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+1] = {
                                      {std::make_tuple(higher_x,        y, 0, h1, x+1,   y),
                                       std::make_tuple(       x, higher_y, 0, h2,   x, y+1),
                                       std::make_tuple(       x,        y, 0, h3,   x,   y)},
                                       position
                                     };
        }
      }
    }
  }
}

void Honeycomb::CreateRectangularCluster()
{
  Matrix3LD h1, h2, h3;
  if (ClusterType == 1){        //rectangle oriented around y-bond
    Translation1 = a1;
    Translation2 = 2*a2-a1;
    h1 = HamInfo.Hx;
    h2 = HamInfo.Hy;
    h3 = HamInfo.Hz;
  } else if (ClusterType == 2){ //rectangle oriented around z-bond
    Translation1 = a1-a2;
    Translation2 = a1+a2;
    h1 = HamInfo.Hy;
    h2 = HamInfo.Hz;
    h3 = HamInfo.Hx;
  } else if (ClusterType == 3){ //rectangle oriented around x-bond
    Translation1 = a2;
    Translation2 = a2-2*a1;
    h1 = HamInfo.Hz;
    h2 = HamInfo.Hx;
    h3 = HamInfo.Hy;
  }

  SublatticeOffset[0] = (3.0*Translation1+1.0*Translation2)/6.0;
  SublatticeOffset[1] = (6.0*Translation1+2.0*Translation2)/6.0;
  SublatticeOffset[2] = (6.0*Translation1+4.0*Translation2)/6.0;
  SublatticeOffset[3] = (3.0*Translation1+5.0*Translation2)/6.0;

  Vector3LD some_vec;
  int higher_x, higher_y, lower_x, lower_y;
  uint flat_index;
  Vector2LD position;
  for (int y=0; y<L2; ++y){
    for (int x=0; x<L1; ++x){
      for (uint sub =0; sub<NumSublattices; ++sub){
        // cout << "unit cell: " << x << " " << y << endl;
        higher_x = (x+1)%L1;
        higher_y = (y+1)%L2;
        //if index becomes -1, wrap it around to l-1
        lower_x = (x-1);
        lower_y = (y-1);
        if (lower_x >= 0){lower_x = lower_x%L1;}
        else if (lower_x < 0){lower_x = L1-1;}
        if (lower_y >= 0){lower_y = lower_y%L2;}
        else if (lower_y < 0){lower_y = L2-1;}
        // cout << lower_x << " " << lower_y << endl;
        position = x*Translation1+y*Translation2+SublatticeOffset[sub];
        BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);

        Cluster[flat_index] = some_vec;
        if (flat_index%NumSublattices==0){
          ClusterInfo[flat_index+0] = {
                                      {std::make_tuple(lower_x,       y, 1, h1, x-1,   y),
                                       std::make_tuple(      x, lower_y, 3, h2,   x, y-1),
                                       std::make_tuple(      x,       y, 1, h3,   x,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+1] = {
                                      {std::make_tuple(higher_x,        y, 0, h1, x+1,   y),
                                       std::make_tuple(       x,        y, 2, h2,   x,   y),
                                       std::make_tuple(       x,        y, 0, h3,   x,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+2] = {
                                      {std::make_tuple(       x,       y, 3, h1,   x,   y),
                                       std::make_tuple(       x,       y, 1, h2,   x,   y),
                                       std::make_tuple(higher_x,       y, 3, h3, x+1,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+3] = {
                                      {std::make_tuple(       x,        y, 2, h1,   x,   y),
                                       std::make_tuple(       x, higher_y, 0, h2,   x, y+1),
                                       std::make_tuple( lower_x,        y, 2, h3, x-1,   y)},
                                       position
                                     };
        }
      }
    }
  }
}

void Honeycomb::CreateC3Cluster()
{
  //C3 cluster with z-bonds arranged vertically
  Translation1 = 2*a1-a2;
  Translation2 = 2*a2-a1;

  SublatticeOffset[0] = (1.0*Translation1+1.0*Translation2)/6.0;
  SublatticeOffset[1] = (3.0*Translation1+1.0*Translation2)/6.0;
  SublatticeOffset[2] = (5.0*Translation1+3.0*Translation2)/6.0;
  SublatticeOffset[3] = (5.0*Translation1+5.0*Translation2)/6.0;
  SublatticeOffset[4] = (3.0*Translation1+5.0*Translation2)/6.0;
  SublatticeOffset[5] = (1.0*Translation1+3.0*Translation2)/6.0;

  Vector3LD some_vec;
  int higher_x, higher_y, lower_x, lower_y;
  uint flat_index;
  Vector2LD position;
  for (int y=0; y<L2; ++y){
    for (int x=0; x<L1; ++x){
      for (uint sub =0; sub<NumSublattices; ++sub){
        // cout << "unit cell: " << x << " " << y << endl;
        higher_x = (x+1)%L1;
        higher_y = (y+1)%L2;
        //if index becomes -1, wrap it around to l-1
        lower_x = (x-1);
        lower_y = (y-1);
        if (lower_x >= 0){lower_x = lower_x%L1;}
        else if (lower_x < 0){lower_x = L1-1;}
        if (lower_y >= 0){lower_y = lower_y%L2;}
        else if (lower_y < 0){lower_y = L2-1;}
        // cout << lower_x << " " << lower_y << endl;
        position = x*Translation1+y*Translation2+SublatticeOffset[sub];
        BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);

        Cluster[flat_index] = some_vec;
        if (flat_index%NumSublattices==0){
          ClusterInfo[flat_index+0] = {
                                      {std::make_tuple(      x,       y, 1, HamInfo.Hx,   x,   y),
                                       std::make_tuple(      x,       y, 5, HamInfo.Hy,   x,   y),
                                       std::make_tuple(lower_x, lower_y, 3, HamInfo.Hz, x-1, y-1)},
                                       position
                                     };
          ClusterInfo[flat_index+1] = {
                                      {std::make_tuple(       x,        y, 0, HamInfo.Hx,   x,   y),
                                       std::make_tuple(       x,  lower_y, 4, HamInfo.Hy,   x, y-1),
                                       std::make_tuple(       x,        y, 2, HamInfo.Hz,   x,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+2] = {
                                      {std::make_tuple(higher_x,       y, 5, HamInfo.Hx, x+1,   y),
                                       std::make_tuple(       x,       y, 3, HamInfo.Hy,   x,   y),
                                       std::make_tuple(       x,       y, 1, HamInfo.Hz,   x,   y)},
                                       position
                                     };
          ClusterInfo[flat_index+3] = {
                                      {std::make_tuple(       x,        y, 4, HamInfo.Hx,   x,   y),
                                       std::make_tuple(       x,        y, 2, HamInfo.Hy,   x,   y),
                                       std::make_tuple(higher_x, higher_y, 0, HamInfo.Hz, x+1, y+1)},
                                       position
                                     };
          ClusterInfo[flat_index+4] = {
                                     {std::make_tuple(       x,       y, 3, HamInfo.Hx,   x,   y),
                                      std::make_tuple(       x,higher_y, 1, HamInfo.Hy,   x, y+1),
                                      std::make_tuple(       x,       y, 5, HamInfo.Hz,   x,   y)},
                                      position
                                    };
          ClusterInfo[flat_index+5] = {
                                     {std::make_tuple( lower_x,        y, 2, HamInfo.Hx, x-1,   y),
                                      std::make_tuple(       x,        y, 0, HamInfo.Hy,   x,   y),
                                      std::make_tuple(       x,        y, 4, HamInfo.Hz,   x,   y)},
                                      position
                                    };
        }
      }
    }
  }
}

void Honeycomb::CreateKekuleCluster(){}      //to be implemented later
void Honeycomb::CalculateClusterOP(){}       //to be implemented later
void Honeycomb::PrintOP(std::ostream &out){} //to be implemented later

void Honeycomb::PrintLatticeParameters(std::ostream &out){
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Honeycomb/Kekule? Which cluster type?\n";
  out << HcOrKekule << " " << ClusterType << "\n";
  out << "Number of Sublattices/Unit Cells (s, l1, l2)\n";
  out << NumSublattices << " " << L1 << " " << L2 << " " << "\n";
}

void Honeycomb::PrintHamiltonianParameters(std::ostream &out){
  out << HamInfo.ParameterOutput;
}

Triangular::Triangular(const uint& l1, const uint& l2, const uint& num_sublattices,
                       const uint& num_defects, Hamiltonia& haminfo,
                       const long double& defect_quad,
                       const long double& defect_octo,
                       const long double& defect_lengthscale):
L1(l1), L2(l2), NumSublattices(num_sublattices), NumSites(l1*l2*num_sublattices),
NumDefects(num_defects),
HamInfo(haminfo),
hField(haminfo.hField),
DefectQuad(defect_quad),
DefectOcto(defect_octo),
DefectLengthScale(defect_lengthscale)
{
  if (not ( (num_defects == 1) || (num_defects == 3) || (num_defects == 9) )){
    std::cerr << "Only 1, 3, or 9 defects are allowed currently." << endl;
    std::cerr << "If you want zero defects, set defect strength to zero." << endl;
    abort();
  }
  if ( (L1!=L2) && ( (abs(defect_quad) > pow(10,-6)) && (abs(defect_octo) > pow(10,-6)) ) ){
    std::cerr << "Only use L1 != L2 if defect strength is zero." << endl;
    abort();
  }
  if ( ((L1<6) || (L2<6)) && ( (abs(defect_quad) > pow(10,-6)) && (abs(defect_octo) > pow(10,-6)) )) {
    std::cerr << "Only use L1,L2<6 if defect strength is zero." << endl;
    abort();
  }

  SublatticeOffset.resize(NumSublattices);

  Cluster.resize(NumSites);
  ClusterInfo.resize(NumSites);
  CreateRhombicCluster();

  Defects.resize(NumDefects);
  CreateDefectPositions();
  AddDefectHamiltonia();

  CreateStripySignMatrices();
}

void Triangular::CreateDefectPositions(){
  // only to be trusted for l1 = l2= 6*n!!
  // only to be trusted for NumDefects = 1,3,9 only!

  double pos = (double)L1/6.0;
  if (NumDefects >= 1){
    Defects[0] = std::make_tuple(3*(uint)pos,3*(uint)pos,0);
  }
  if (NumDefects >= 3){
    Defects[1] = std::make_tuple((uint)pos,    (uint)pos,0);
    Defects[2] = std::make_tuple(5*(uint)pos,5*(uint)pos,0);
  }
  if (NumDefects == 9){
    Defects[3] = std::make_tuple(3*(uint)pos,  (uint)pos,0);
    Defects[4] = std::make_tuple(5*(uint)pos,  (uint)pos,0);
    Defects[5] = std::make_tuple((uint)pos,  3*(uint)pos,0);
    Defects[6] = std::make_tuple(5*(uint)pos,3*(uint)pos,0);
    Defects[7] = std::make_tuple((uint)pos,  5*(uint)pos,0);
    Defects[8] = std::make_tuple(3*(uint)pos,5*(uint)pos,0);
  }
}

void Triangular::CreateRhombicCluster()
// T1 = a1-a2, T2 = a1
{
  Translation1 = a1-a2;
  Translation2 = a1;

  SublatticeOffset[0] = (Translation1+Translation2)/3.0;

  Vector3LD some_vec;
  int higher_x, higher_y, lower_x, lower_y;
  uint flat_index;
  Vector2LD position;
  for (int y=0; y<L2; ++y){
    for (int x=0; x<L1; ++x){
      for (uint sub =0; sub<NumSublattices; ++sub){
        // cout << "unit cell: " << x << " " << y << endl;
        higher_x = (x+1)%L1;
        higher_y = (y+1)%L2;
        //if index becomes -1, wrap it around to l-1
        lower_x = (x-1);
        lower_y = (y-1);
        if (lower_x >= 0){lower_x = lower_x%L1;}
        else if (lower_x < 0){lower_x = L1-1;}
        if (lower_y >= 0){lower_y = lower_y%L2;}
        else if (lower_y < 0){lower_y = L2-1;}
        // cout << lower_x << " " << lower_y << endl;
        position = x*Translation1+y*Translation2+SublatticeOffset[sub];
        BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);
        ClusterInfo[flat_index] = {
                                    {std::make_tuple(higher_x,        y, sub, HamInfo.Hz, x+1, y  ),
                                     std::make_tuple(       x, higher_y, sub, HamInfo.Hy, x  , y+1),
                                     std::make_tuple( lower_x, higher_y, sub, HamInfo.Hx, x-1, y+1),
                                     std::make_tuple( lower_x,        y, sub, HamInfo.Hz, x-1, y  ),
                                     std::make_tuple(       x,  lower_y, sub, HamInfo.Hy, x  , y-1),
                                     std::make_tuple(higher_x,  lower_y, sub, HamInfo.Hx, x+1, y-1)},
                                     position
        };
        Cluster[flat_index] = some_vec;
      }
    }
  }
}

void Triangular::AddDefectHamiltonia()
{
  Vector2LD rdefect, rsite, rNN, rbond, rseparation;
  auto function = Lorentzian;
  long double strength;
  //displacements from the defect position
  // first is the real defect
  // then follows the 6 nearest image charges from PBC
  // further image charges are assumed to be non-zero. this should be changed if defect
  // is long ranged
  std::vector<std::tuple<int, int>> defect_index = {
                                                  {  0,  0},
                                                  {  1,  0},
                                                  { -1,  0},
                                                  {  0,  1},
                                                  {  0, -1},
                                                  { -1,  1},
                                                  {  1, -1},
                                                  {  1,  1},
                                                  { -1,  2},
                                                  { -2,  1},
                                                  {  2, -1},
                                                  {  1, -2},
                                                  { -1, -1}
                                                };
  for (auto& indices:Defects){
    rdefect = get<0>(indices)*Translation1+get<1>(indices)*Translation2+SublatticeOffset[get<2>(indices)]; //+
    for (uint i=0;i<Cluster.size();i++){
      rsite = ClusterInfo[i].Position;
      //check length would go here
      for (auto& nninfo:ClusterInfo[i].NearestNeighbours){
        //I'm using the non-wrapped position indices here.
        //this is to deal with the (periodic) bonds that go beyond the finite size cluster
        rNN = get<4>(nninfo)*Translation1+get<5>(nninfo)*Translation2+SublatticeOffset[get<2>(nninfo)];
        rbond = (rsite+rNN)/2.0;
        rseparation = rbond - rdefect;
        for (auto &j:defect_index){
            strength = function(1,
            (rseparation - get<0>(j)*L1*Translation1 - get<1>(j)*L2*Translation2).norm(),
            DefectLengthScale);
            get<3>(nninfo)(1,1) += +1*DefectOcto*strength;
            get<3>(nninfo)(0,0) += +1*DefectQuad*strength;
            get<3>(nninfo)(2,2) += +1*DefectQuad*strength;
        }
      }
    }
  }
}

void Triangular::CreateStripySignMatrices()
{
  ArrayXXLD signs_Y(L1,L2);
  ArrayXXLD signs_Z(L1,L2);
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      x%2 == 0 ? signs_Y(x,y) = 1 : signs_Y(x,y) = -1;
      y%2 == 0 ? signs_Z(x,y) = 1 : signs_Z(x,y) = -1;
    }
  }
  ArrayXXLD signs_X = signs_Y*signs_Z;

  Eigen::Map<ArrayXLD> ssx(signs_X.data(), signs_X.size()); //reshapes matrix (column-major)
  Eigen::Map<ArrayXLD> ssy(signs_Y.data(), signs_Y.size());
  Eigen::Map<ArrayXLD> ssz(signs_Z.data(), signs_Z.size());

  StripySignsX = ssx;
  StripySignsY = ssy;
  StripySignsZ = ssz;
}

void Triangular::CalculateClusterOP()
{
  uint flat_index;
  Vector3LD spin, fmop;
  Eigen::Matrix<long double, 3, 2> stripyopmatrix;
  fmop << 0,0,0;
  stripyopmatrix << 0, 0,
                    0, 0,
                    0, 0;

  for (uint flat_index=0; flat_index<NumSites; ++flat_index){
    // local_energy=0;
    //calculate energy
    spin = Cluster[flat_index];
    fmop += spin;
    //calculate three stripy OPs
    // cout << fmop << endl;
    stripyopmatrix(0,0) += StripySignsX(flat_index)*spin(0); stripyopmatrix(0,1) += StripySignsX(flat_index)*spin(2);
    stripyopmatrix(1,0) += StripySignsY(flat_index)*spin(0); stripyopmatrix(1,1) += StripySignsY(flat_index)*spin(2);
    stripyopmatrix(2,0) += StripySignsZ(flat_index)*spin(0); stripyopmatrix(2,1) += StripySignsZ(flat_index)*spin(2);
    // cout << fmop << endl;
  }

  ClusterFMOP = fmop;
  ClusterStripyOPMatrix = stripyopmatrix;

  // SelectStripyOP();
  AverageStripyOP();
}

void Triangular::SelectStripyOP()
{
  Eigen::MatrixXd::Index indices[1];
  long double stripymax = ClusterStripyOPMatrix.rowwise().norm().maxCoeff(&indices[0]);
  ClusterStripyOP = ClusterStripyOPMatrix.row(indices[0]);

  Vector3LD combinedop;
  combinedop << ClusterStripyOP(0), ClusterFMOP(1),ClusterStripyOP(1);
  ClusterCombinedOP = combinedop;
}

void Triangular::AverageStripyOP()
{
  ClusterStripyOP = ClusterStripyOPMatrix.colwise().sum().transpose();

  Vector3LD combinedop;
  combinedop << ClusterStripyOP(0), ClusterFMOP(1),ClusterStripyOP(1);
  ClusterCombinedOP = combinedop;
}

void Triangular::PrintLatticeParameters(std::ostream &out){
  const uint type = 2;
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Which cluster type?\n";
  out << type << "\n";
  out << "Number of Sublattices/Unit Cells (s, l1, l2)\n";
  out << NumSublattices << " " << L1 << " " << L2 << " " << "\n";
}

void Triangular::PrintHamiltonianParameters(std::ostream &out){
  out << HamInfo.ParameterOutput;
  out << "Defect quad/defect octo/defect length scale/number of defects\n";
  out << DefectQuad << "/" << DefectOcto << "/" << DefectLengthScale << "/" << NumDefects << "\n";
}

void Triangular::PrintOP(std::ostream &out)
{
  // SelectStripyOP();
  AverageStripyOP();
  out << "-------------------------------Final configuration observables--------------------------------\n";
  out << "Order parameters. R1: (FM_x, FM_y, FM_z), R2: (stripy_x, stripy_z),  R3: (stripy_x, FM_y, stripy_z)\n";
  out << std::setprecision(14) << ClusterFMOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterStripyOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterCombinedOP.transpose()/(long double)NumSites << "\n";
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// void TriangularLattice::PrintThermalObservables(std::ostream &out){
//   long double ns = NumSites;
//   long double ns2 = pow(NumSites,2);
//   long double ns3 = pow(NumSites,3);
//   long double ns4 = pow(NumSites,4);
//
//   out << "-------------------------------Thermal-averaged observables-----------------------------\n";
//   out << "Energy cumulants (C: E, E2, E3, E4) \n";
//   out << std::setprecision(14) << EBar/ns << " " << E2Bar/ns2 << " " << E3Bar/ns3 << " " << E4Bar/ns4 << "\n";
//   out << "Order parameter cumulants (R: FM, Perp, Par, Combined), (C: |m|, |m|2, |m|4)\n";
//   out << std::setprecision(14) << FMNorm/ns << " " << FMNorm2/ns2 << " " << FMNorm4/ns4 << "\n";
//   out << std::setprecision(14) << PerpNorm/ns << " " << PerpNorm2/ns2 << " " <<  PerpNorm4/ns4 << "\n";
//   out << std::setprecision(14) << ParNorm/ns << " " << ParNorm2/ns2 << " " << ParNorm4/ns4 << "\n";
//   out << std::setprecision(14) << CombinedNorm/ns << " " << CombinedNorm2/ns2 << " " << CombinedNorm4/ns4 << "\n";
// }
//
//
// void TriangularLattice::SampleConfiguration(double &temp, const uint& max_sweeps,
//                                             const uint& sampling_time){
//   // cout << temp << " " << temp << endl;
//   long double m_e = 0;
//   long double m_e2 = 0;
//   long double m_e3 = 0;
//   long double m_e4 = 0;
//
//   long double m_fm_norm =0;
//   long double m_perp_norm =0;
//   long double m_par_norm =0;
//   long double m_combined_norm =0;
//
//   long double m_fm_norm2 =0;
//   long double m_perp_norm2 =0;
//   long double m_par_norm2 =0;
//   long double m_combined_norm2 =0;
//
//   long double m_fm_norm4 =0;
//   long double m_perp_norm4 =0;
//   long double m_par_norm4 =0;
//   long double m_combined_norm4 =0;
//
//   long double energy, fm_norm, perp_norm, par_norm, combined_norm;
//   uint sweep = 0;
//   uint samples = 0;
//   uint accept;
//
//   while (sweep < max_sweeps){
//     DoTheSweeps(temp, accept);
//
//     if (sweep%sampling_time == 0){
//       // cout << samples << " " << ClusterEnergy/NumSites << endl;
//
//       // SelectStripyOP();
//       CalculateClusterEnergyandOP();
//       energy = ClusterEnergy;
//       AverageStripyOP();
//
//       m_e  += energy;
//       m_e2 += pow(energy,2);
//       m_e3 += pow(energy,3);
//       m_e4 += pow(energy,4);
//
//       fm_norm   = ClusterFMOP.norm();
//       perp_norm = ClusterStripyOP.norm();
//       par_norm  = abs(ClusterFMOP(1));
//       combined_norm = ClusterCombinedOP.norm();
//
//       m_fm_norm    += fm_norm;
//       m_perp_norm  += perp_norm;
//       m_par_norm   += par_norm;
//       m_combined_norm += combined_norm;
//
//       m_fm_norm2   += pow(fm_norm,2);
//       m_perp_norm2 += pow(perp_norm,2);
//       m_par_norm2  += pow(par_norm,2);
//       m_combined_norm2  += pow(combined_norm,2);
//
//       m_fm_norm4   += pow(fm_norm,4);
//       m_perp_norm4 += pow(perp_norm,4);
//       m_par_norm4  += pow(par_norm,4);
//       m_combined_norm4  += pow(combined_norm,4);
//       ++samples;
//
//     }
//     ++sweep;
//   }
//
//   EBar  = m_e/(long double)samples;
//   // cout << "....." << endl;
//   // cout << EBar << endl;
//   // cout << "....." << endl;
//
//   E2Bar = m_e2/(long double)samples;
//   E3Bar = m_e3/(long double)samples;
//   E4Bar = m_e4/(long double)samples;
//   // cout << std::setprecision(14) << 0  << " " << ebar << " " << e2bar << endl;
//
//   FMNorm       = m_fm_norm/(long double)samples;
//   PerpNorm     = m_perp_norm/(long double)samples;
//   ParNorm      = m_par_norm/(long double)samples;
//   CombinedNorm = m_combined_norm/(long double)samples;
//
//   FMNorm2       = m_fm_norm2/(long double)samples;
//   PerpNorm2     = m_perp_norm2/(long double)samples;
//   ParNorm2      = m_par_norm2/(long double)samples;
//   CombinedNorm2 = m_combined_norm2/(long double)samples;
//
//   FMNorm4       = m_fm_norm4/(long double)samples;
//   PerpNorm4     = m_perp_norm4/(long double)samples;
//   ParNorm4      = m_par_norm4/(long double)samples;
//   CombinedNorm4 = m_combined_norm4/(long double)samples;
// }
