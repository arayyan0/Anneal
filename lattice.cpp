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

  // DebugHamiltonians();

  SSFPoints = {gamma,
               gammapx, gammapy, gammapz,
               kx, ky, kz,
               kx/2.0, ky/2.0, kz/2.0,
               mx, my, mz,
               2.0/3.0*mx, 2.0/3.0*my, 2.0/3.0*mz};

  SSFPointsLabel = {"Gamma   ",
                    "GammaPx ",
                    "GammaPy ",
                    "GammaPz ",
                    "Kx      ",
                    "Ky      ",
                    "Kz      ",
                    "Kx/2    ",
                    "Ky/2    ",
                    "Kz/2    ",
                    "Mx      ",
                    "My      ",
                    "Mz      ",
                    "2Mx/3   ",
                    "2My/3   ",
                    "2Mz/3   "
  };

}

void Triangular::DebugHamiltonians(){
  uint flat1, flat2, s2, x2, y2;
  uint wrongham=0;
  double wrongamout=0.0;
  Matrix3LD hamdiff;
  for (uint y1=0; y1<L2; ++y1){
    for (uint x1=0; x1<L1; ++x1){
      for (uint s1=0; s1<NumSublattices;++s1){
        BravaisIndicesToFlat(s1, x1, NumSublattices, y1, L1*NumSublattices, flat1);
        std::cerr << "site: " << x1 << " " << y1 << " " << s1 << endl;
        std::cerr << "---------entering NN loop --------" << endl;
        for (auto &j:ClusterInfo[flat1].NearestNeighbours){
          x2 = get<0>(j);
          y2 = get<1>(j);
          s2 = get<2>(j);
          BravaisIndicesToFlat(s2, x2, NumSublattices, y2, L1*NumSublattices, flat2);
          std::cerr << "NN: " << x2 << " " << y2 << " " << s2 << endl;
          /////////////////////output bond Hamiltonian
          if ((x2 == 7) && (y2 == 6)){
            std::cerr << get<3>(j)<< endl;
          }
          // ///////////////////find other pair's Hamiltonian and output diff
          //
          for (auto &jj:ClusterInfo[flat2].NearestNeighbours){
            if ((x1 == get<0>(jj)) && (y1 == get<1>(jj)) && (s1 == get<2>(jj))){
              hamdiff = get<3>(j) - get<3>(jj);
              if (hamdiff.norm()> 0){
                // std::cerr << "oof, wrong hamiltonian for: (" << x1 <<","<<y1<<") and ("<< x2 <<","<<y2<<") w/ norm "<< hamdiff.norm()<< endl;
                wrongham++;
                wrongamout+=hamdiff.norm();
              }
            }
          }
        //
        }
        std::cerr << "---------finished NN loop --------" << endl;
      }
    }
  }
  std::cerr << (double)wrongham/6.0/(double)NumSites*100 << "% wrong bonds, total wrong is " << wrongamout << endl;
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
                                                  // {  1,  0},
                                                  // { -1,  0},
                                                  // {  0,  1},
                                                  // {  0, -1},
                                                  // { -1,  1},
                                                  // {  1, -1},
                                                  // {  1,  1},
                                                  // { -1,  2},
                                                  // { -2,  1},
                                                  // {  2, -1},
                                                  // {  1, -2},
                                                  // { -1, -1}
                                                };
  for (auto& indices:Defects){
    rdefect = get<0>(indices)*Translation1+get<1>(indices)*Translation2+SublatticeOffset[get<2>(indices)]; //+
    for (uint i=0;i<Cluster.size();i++){
    // for (uint i=0;i<1;i++){
      rsite = ClusterInfo[i].Position;
      // cout << rsite.transpose() << endl;
      //check length would go here
      // std::cerr << "site: " << i << endl;
      // std::cerr << "_______________________________entering NN loop" << endl;
      for (auto& nninfo:ClusterInfo[i].NearestNeighbours){
        //I'm using the non-wrapped position indices here.
        //this is to deal with the (periodic) bonds that go beyond the finite size cluster
        rNN = get<4>(nninfo)*Translation1+get<5>(nninfo)*Translation2+SublatticeOffset[get<2>(nninfo)];
        rbond = (rsite+rNN)/2.0;
        rseparation = rbond - rdefect;
        // std::cerr << rseparation.transpose() << " " << rseparation.norm() << endl;
        for (auto &j:defect_index){
            strength = function(
            (rseparation - get<0>(j)*L1*Translation1 - get<1>(j)*L2*Translation2).norm(),
            DefectLengthScale);
            get<3>(nninfo)(1,1) += +1*DefectOcto*strength;
            get<3>(nninfo)(0,0) += +1*DefectQuad*strength;
            get<3>(nninfo)(2,2) += +1*DefectQuad*strength;
        }
      }
      // std::cerr << "_______________________________leaving NN loop" << endl;
    }
  }
}


void Triangular::CalculateClusterOP()
{
  Vector3LD SS_ij;
  Vector2LD r_ij;

  Vector3LDc SS_q;

  vector<Vector3LD> clusterssf(SSFPoints.size(),Vector3LD::Zero());
  for (uint q=0; q<SSFPoints.size(); q++){
    SS_q = Vector3LDc::Zero();
    for (uint i =0; i<NumSites;i++){
      for (uint j =0; j<NumSites;j++){
        SS_ij = Cluster[i].cwiseProduct(Cluster[j]);
        r_ij = ClusterInfo[i].Position - ClusterInfo[j].Position;
        SS_q += std::polar((long double)1.0,-SSFPoints[q].dot(r_ij))*SS_ij;
      }
    }
    clusterssf[q]=SS_q.cwiseAbs();
  }
  ClusterSSf = clusterssf;
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
  out << "-------------------------------Final configuration observables--------------------------------\n";
  Vector3LD whoa;
  out << "             Quad             Octo\n";
  for(uint j=0; j <SSFPoints.size(); j++){
    whoa = ClusterSSf[j]/(long double)NumSites/(long double)NumSites;
    out << SSFPointsLabel[j] << sqrt(whoa(0) +whoa(2)) << " " << sqrt(whoa(1))  << endl;
  }
}

void Triangular::PrintThermalObservables(std::ostream &out){
  long double ns = NumSites;
  long double ns2 = pow(NumSites,2);
  Vector3LD whoa;
  out << std::setprecision(14);
  out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  out << "Energy cumulants: <E>, <E2> \n";
  out << EBar/ns << " " << E2Bar/ns2 << "\n";
  out << "Magnetic order parameter: M(Q) = 1/N sqrt( sum_ij exp[-i Q.(Ri - Rj)] <Si^a Sj^a> )\n";
  out << "             Quad             Octo\n";
  for(uint j=0; j <SSFPoints.size(); j++){
    whoa = SSfBar[j]/ns2;
    out << SSFPointsLabel[j] << sqrt(whoa(0) +whoa(2)) << " " << sqrt(whoa(1))  << endl;
  }
}
