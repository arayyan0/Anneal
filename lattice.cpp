///  @file     lattice.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the types of Lattices
#include "lattice.hpp"

// Honeycomb::Honeycomb(const uint& hc_or_kek, const uint& type,
//                      const uint& num_sublattices, const uint& l1, const uint& l2,
//                      Hamiltonia& haminfo
//                    ):
//  HcOrKekule(hc_or_kek), ClusterType(type), NumSublattices(num_sublattices), L1(l1), L2(l2),
//  NumUnitCells(l1*l2), NumSites(l1*l2*num_sublattices),
//  HamInfo(haminfo), hField(haminfo.hField)
// {
//   SublatticeOffset.resize(NumSublattices);
//
//   Cluster = Eigen::Matrix<long double, 3, Eigen::Dynamic>::Zero(3, NumSites);
//   Position = Eigen::Matrix<long double, 2, Eigen::Dynamic>::Zero(2, NumSites);
//
//   NumNeighbours = 3;
//   FlatIndex.resize(NumNeighbours*NumSites);
//   Hamiltonians.resize(NumNeighbours*NumSites);
//
//   ClusterInfo.resize(NumSites);
//
//   if (HcOrKekule == 0){
//     if      (NumSublattices == 2) CreateRhombicCluster();
//     else if (NumSublattices == 4) CreateRectangularCluster();
//     else if (NumSublattices == 6) CreateC3Cluster();
//   }
//   else if (HcOrKekule == 1){
//     if (NumSublattices == 6) CreateKekuleCluster();
//   }
//
// }
//
// void Honeycomb::CreateRhombicCluster()
// {
//   Matrix3LD h1, h2, h3;
//   if (ClusterType == 1){ //rhombus oriented around z-bond
//     Translation1 = a1;
//     Translation2 = a2;
//     h1 = HamInfo.Hx;
//     h2 = HamInfo.Hy;
//     h3 = HamInfo.Hz;
//   } else if (ClusterType == 2){ //rhombus oriented around x-bond
//     Translation1 = a1-a2;
//     Translation2 = a1;
//     h1 = HamInfo.Hy;
//     h2 = HamInfo.Hz;
//     h3 = HamInfo.Hx;
//   } else if (ClusterType == 3){ //rhombus oriented around y-bond
//     Translation1 = a2;
//     Translation2 = a2-a1;
//     h1 = HamInfo.Hz;
//     h2 = HamInfo.Hx;
//     h3 = HamInfo.Hy;
//   }
//
//   SublatticeOffset[0] = (1.0*Translation1+1.0*Translation2)/3.0;
//   SublatticeOffset[1] = (2.0*Translation1+2.0*Translation2)/3.0;
//
//   Vector3LD some_vec;
//   int higher_x, higher_y, lower_x, lower_y;
//   uint flat_index;
//   uint nn1, nn2, nn3, nn4, nn5, nn6;
//   Vector2LD position;
//   for (int y=0; y<L2; ++y){
//     for (int x=0; x<L1; ++x){ //both x,y should NOT be ints as they may become negative
//       for (uint sub =0; sub<NumSublattices; ++sub){
//         // cout << "unit cell: " << x << " " << y << endl;
//         higher_x = (x+1)%L1;
//         higher_y = (y+1)%L2;
//         //if index becomes -1, wrap it around to l-1
//         lower_x = (x-1);
//         lower_y = (y-1);
//         if (lower_x >= 0){lower_x = lower_x%L1;}
//         else if (lower_x < 0){lower_x = L1-1;}
//         if (lower_y >= 0){lower_y = lower_y%L2;}
//         else if (lower_y < 0){lower_y = L2-1;}
//         // cout << lower_x << " " << lower_y << endl;
//         position = x*Translation1+y*Translation2+SublatticeOffset[sub];
//         BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);
//         Position.col(flat_index) = position;
//
//         BravaisIndicesToFlat(1, lower_x, NumSublattices,       y, L1*NumSublattices, nn1);
//         BravaisIndicesToFlat(1,       x, NumSublattices, lower_y, L1*NumSublattices, nn2);
//         BravaisIndicesToFlat(1,       x, NumSublattices,       y, L1*NumSublattices, nn3);
//         BravaisIndicesToFlat(0, higher_x, NumSublattices,        y, L1*NumSublattices, nn4);
//         BravaisIndicesToFlat(0,        x, NumSublattices, higher_y, L1*NumSublattices, nn5);
//         BravaisIndicesToFlat(0,        x, NumSublattices,        y, L1*NumSublattices, nn6);
//
//         if (flat_index%NumSublattices==0){
//           ClusterInfo[flat_index+0] = {
//                                       {std::make_tuple(lower_x,       y, 1, h1, x-1,   y, nn1),
//                                        std::make_tuple(      x, lower_y, 1, h2,   x, y-1, nn2),
//                                        std::make_tuple(      x,       y, 1, h3,   x,   y, nn3)}
//                                      };
//           ClusterInfo[flat_index+1] = {
//                                       {std::make_tuple(higher_x,        y, 0, h1, x+1,   y, nn4),
//                                        std::make_tuple(       x, higher_y, 0, h2,   x, y+1, nn5),
//                                        std::make_tuple(       x,        y, 0, h3,   x,   y, nn6)}
//                                      };
//         }
//       }
//     }
//   }
// }
//
// void Honeycomb::CreateRectangularCluster()
// {
//   Matrix3LD h1, h2, h3;
//   if (ClusterType == 1){        //rectangle oriented around y-bond
//     Translation1 = a1;
//     Translation2 = 2*a2-a1;
//     h1 = HamInfo.Hx;
//     h2 = HamInfo.Hy;
//     h3 = HamInfo.Hz;
//   } else if (ClusterType == 2){ //rectangle oriented around z-bond
//     Translation1 = a1-a2;
//     Translation2 = a1+a2;
//     h1 = HamInfo.Hy;
//     h2 = HamInfo.Hz;
//     h3 = HamInfo.Hx;
//   } else if (ClusterType == 3){ //rectangle oriented around x-bond
//     Translation1 = a2;
//     Translation2 = a2-2*a1;
//     h1 = HamInfo.Hz;
//     h2 = HamInfo.Hx;
//     h3 = HamInfo.Hy;
//   }
//
//   SublatticeOffset[0] = (3.0*Translation1+1.0*Translation2)/6.0;
//   SublatticeOffset[1] = (6.0*Translation1+2.0*Translation2)/6.0;
//   SublatticeOffset[2] = (6.0*Translation1+4.0*Translation2)/6.0;
//   SublatticeOffset[3] = (3.0*Translation1+5.0*Translation2)/6.0;
//
//   Vector3LD some_vec;
//   int higher_x, higher_y, lower_x, lower_y;
//   uint flat_index;
//   Vector2LD position;
//   uint nn1, nn2, nn3, nn4, nn5, nn6, nn7, nn8, nn9, nn10, nn11, nn12;
//   for (int y=0; y<L2; ++y){
//     for (int x=0; x<L1; ++x){
//       for (uint sub =0; sub<NumSublattices; ++sub){
//         // cout << "unit cell: " << x << " " << y << endl;
//         higher_x = (x+1)%L1;
//         higher_y = (y+1)%L2;
//         //if index becomes -1, wrap it around to l-1
//         lower_x = (x-1);
//         lower_y = (y-1);
//         if (lower_x >= 0){lower_x = lower_x%L1;}
//         else if (lower_x < 0){lower_x = L1-1;}
//         if (lower_y >= 0){lower_y = lower_y%L2;}
//         else if (lower_y < 0){lower_y = L2-1;}
//         // cout << lower_x << " " << lower_y << endl;
//         position = x*Translation1+y*Translation2+SublatticeOffset[sub];
//         BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);
//         Position.col(flat_index) = position;
//
//         BravaisIndicesToFlat(1, lower_x, NumSublattices,       y, L1*NumSublattices, nn1);
//         BravaisIndicesToFlat(3,       x, NumSublattices, lower_y, L1*NumSublattices, nn2);
//         BravaisIndicesToFlat(1,       x, NumSublattices,       y, L1*NumSublattices, nn3);
//         BravaisIndicesToFlat(0, higher_x, NumSublattices, y, L1*NumSublattices, nn4);
//         BravaisIndicesToFlat(2,        x, NumSublattices, y, L1*NumSublattices, nn5);
//         BravaisIndicesToFlat(0,        x, NumSublattices, y, L1*NumSublattices, nn6);
//         BravaisIndicesToFlat(3,        x, NumSublattices, y, L1*NumSublattices, nn7);
//         BravaisIndicesToFlat(1,        x, NumSublattices, y, L1*NumSublattices, nn8);
//         BravaisIndicesToFlat(3, higher_x, NumSublattices, y, L1*NumSublattices, nn9);
//         BravaisIndicesToFlat(2,       x, NumSublattices,        y, L1*NumSublattices, nn10);
//         BravaisIndicesToFlat(0,       x, NumSublattices, higher_y, L1*NumSublattices, nn11);
//         BravaisIndicesToFlat(2, lower_x, NumSublattices,        y, L1*NumSublattices, nn12);
//
//         if (flat_index%NumSublattices==0){
//           ClusterInfo[flat_index+0] = {
//                                       {std::make_tuple(lower_x,       y, 1, h1, x-1,   y,nn1),
//                                        std::make_tuple(      x, lower_y, 3, h2,   x, y-1,nn2),
//                                        std::make_tuple(      x,       y, 1, h3,   x,   y,nn3)}
//                                      };
//           ClusterInfo[flat_index+1] = {
//                                       {std::make_tuple(higher_x,        y, 0, h1, x+1,   y,nn4),
//                                        std::make_tuple(       x,        y, 2, h2,   x,   y,nn5),
//                                        std::make_tuple(       x,        y, 0, h3,   x,   y,nn6)}
//                                      };
//           ClusterInfo[flat_index+2] = {
//                                       {std::make_tuple(       x,       y, 3, h1,   x,   y,nn7),
//                                        std::make_tuple(       x,       y, 1, h2,   x,   y,nn8),
//                                        std::make_tuple(higher_x,       y, 3, h3, x+1,   y,nn9)}
//                                      };
//           ClusterInfo[flat_index+3] = {
//                                       {std::make_tuple(       x,        y, 2, h1,   x,   y,nn10),
//                                        std::make_tuple(       x, higher_y, 0, h2,   x, y+1,nn11),
//                                        std::make_tuple( lower_x,        y, 2, h3, x-1,   y,nn12)}
//                                      };
//         }
//       }
//     }
//   }
// }
//
// void Honeycomb::CreateC3Cluster()
// {
//   //C3 cluster with z-bonds arranged vertically
//   Translation1 = 2*a1-a2;
//   Translation2 = 2*a2-a1;
//
//   SublatticeOffset[0] = (1.0*Translation1+1.0*Translation2)/6.0;
//   SublatticeOffset[1] = (3.0*Translation1+1.0*Translation2)/6.0;
//   SublatticeOffset[2] = (5.0*Translation1+3.0*Translation2)/6.0;
//   SublatticeOffset[3] = (5.0*Translation1+5.0*Translation2)/6.0;
//   SublatticeOffset[4] = (3.0*Translation1+5.0*Translation2)/6.0;
//   SublatticeOffset[5] = (1.0*Translation1+3.0*Translation2)/6.0;
//
//   Vector3LD some_vec;
//   int higher_x, higher_y, lower_x, lower_y;
//   uint flat_index;
//   Vector2LD position;
//   uint nn1, nn2, nn3, nn4, nn5, nn6, nn7, nn8, nn9, nn10, nn11, nn12, nn13, nn14, nn15, nn16, nn17, nn18;
//   for (int y=0; y<L2; ++y){
//     for (int x=0; x<L1; ++x){
//       for (uint sub =0; sub<NumSublattices; ++sub){
//         // cout << "unit cell: " << x << " " << y << endl;
//         higher_x = (x+1)%L1;
//         higher_y = (y+1)%L2;
//         //if index becomes -1, wrap it around to l-1
//         lower_x = (x-1);
//         lower_y = (y-1);
//         if (lower_x >= 0){lower_x = lower_x%L1;}
//         else if (lower_x < 0){lower_x = L1-1;}
//         if (lower_y >= 0){lower_y = lower_y%L2;}
//         else if (lower_y < 0){lower_y = L2-1;}
//         // cout << lower_x << " " << lower_y << endl;
//         position = x*Translation1+y*Translation2+SublatticeOffset[sub];
//         BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);
//         Position.col(flat_index) = position;
//
//         BravaisIndicesToFlat(1,       x, NumSublattices,       y, L1*NumSublattices, nn1);
//         BravaisIndicesToFlat(5,       x, NumSublattices,       y, L1*NumSublattices, nn2);
//         BravaisIndicesToFlat(3, lower_x, NumSublattices, lower_y, L1*NumSublattices, nn3);
//         BravaisIndicesToFlat(0, x, NumSublattices,       y, L1*NumSublattices, nn4);
//         BravaisIndicesToFlat(4, x, NumSublattices, lower_y, L1*NumSublattices, nn5);
//         BravaisIndicesToFlat(2, x, NumSublattices,       y, L1*NumSublattices, nn6);
//         BravaisIndicesToFlat(5, higher_x, NumSublattices, y, L1*NumSublattices, nn7);
//         BravaisIndicesToFlat(3,        x, NumSublattices, y, L1*NumSublattices, nn8);
//         BravaisIndicesToFlat(1,        x, NumSublattices, y, L1*NumSublattices, nn9);
//         BravaisIndicesToFlat(4,        x, NumSublattices,        y, L1*NumSublattices, nn10);
//         BravaisIndicesToFlat(2,        x, NumSublattices,        y, L1*NumSublattices, nn11);
//         BravaisIndicesToFlat(0, higher_x, NumSublattices, higher_y, L1*NumSublattices, nn12);
//         BravaisIndicesToFlat(3, x, NumSublattices,        y, L1*NumSublattices, nn13);
//         BravaisIndicesToFlat(1, x, NumSublattices, higher_y, L1*NumSublattices, nn14);
//         BravaisIndicesToFlat(5, x, NumSublattices,        y, L1*NumSublattices, nn15);
//         BravaisIndicesToFlat(2, lower_x, NumSublattices, y, L1*NumSublattices, nn16);
//         BravaisIndicesToFlat(0,       x, NumSublattices, y, L1*NumSublattices, nn17);
//         BravaisIndicesToFlat(4,       x, NumSublattices, y, L1*NumSublattices, nn18);
//
//         if (flat_index%NumSublattices==0){
//           ClusterInfo[flat_index+0] = {
//                                       {std::make_tuple(      x,       y, 1, HamInfo.Hx,   x,   y,nn1),
//                                        std::make_tuple(      x,       y, 5, HamInfo.Hy,   x,   y,nn2),
//                                        std::make_tuple(lower_x, lower_y, 3, HamInfo.Hz, x-1, y-1,nn3)}
//                                      };
//           ClusterInfo[flat_index+1] = {
//                                       {std::make_tuple(       x,        y, 0, HamInfo.Hx,   x,   y,nn4),
//                                        std::make_tuple(       x,  lower_y, 4, HamInfo.Hy,   x, y-1,nn5),
//                                        std::make_tuple(       x,        y, 2, HamInfo.Hz,   x,   y,nn6)}
//                                      };
//           ClusterInfo[flat_index+2] = {
//                                       {std::make_tuple(higher_x,       y, 5, HamInfo.Hx, x+1,   y,nn7),
//                                        std::make_tuple(       x,       y, 3, HamInfo.Hy,   x,   y,nn8),
//                                        std::make_tuple(       x,       y, 1, HamInfo.Hz,   x,   y,nn9)}
//                                      };
//           ClusterInfo[flat_index+3] = {
//                                       {std::make_tuple(       x,        y, 4, HamInfo.Hx,   x,   y,nn10),
//                                        std::make_tuple(       x,        y, 2, HamInfo.Hy,   x,   y,nn11),
//                                        std::make_tuple(higher_x, higher_y, 0, HamInfo.Hz, x+1, y+1,nn12)}
//                                      };
//           ClusterInfo[flat_index+4] = {
//                                      {std::make_tuple(       x,       y, 3, HamInfo.Hx,   x,   y,nn13),
//                                       std::make_tuple(       x,higher_y, 1, HamInfo.Hy,   x, y+1,nn14),
//                                       std::make_tuple(       x,       y, 5, HamInfo.Hz,   x,   y,nn15)}
//                                     };
//           ClusterInfo[flat_index+5] = {
//                                      {std::make_tuple( lower_x,        y, 2, HamInfo.Hx, x-1,   y,nn16),
//                                       std::make_tuple(       x,        y, 0, HamInfo.Hy,   x,   y,nn17),
//                                       std::make_tuple(       x,        y, 4, HamInfo.Hz,   x,   y,nn18)}
//                                     };
//         }
//       }
//     }
//   }
// }
//
// void Honeycomb::CreateKekuleCluster(){}                      //to be implemented later
// void Honeycomb::CalculateClusterOP(){}                       //to be implemented later
// void Honeycomb::PrintOP(std::ostream &out){}                 //to be implemented later
// void Honeycomb::PrintThermalObservables(std::ostream &out){} //to be implemented later
//
// void Honeycomb::PrintLatticeParameters(std::ostream &out){
//   out << "-------------------------Simulation Parameters-------------------------\n";
//   out << "Honeycomb/Kekule? Which cluster type?\n";
//   out << HcOrKekule << " " << ClusterType << "\n";
//   out << "Number of Sublattices/Unit Cells (s, l1, l2)\n";
//   out << NumSublattices << " " << L1 << " " << L2 << " " << "\n";
// }
//
// void Honeycomb::PrintHamiltonianParameters(std::ostream &out){
//   out << HamInfo.ParameterOutput;
// }


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
  ClusterType = 2;

  Cluster = Eigen::Matrix<long double, 3, Eigen::Dynamic>::Zero(3, NumSites);
  Position = Eigen::Matrix<long double, 2, Eigen::Dynamic>::Zero(2, NumSites);

  NumNeighbours = 6;
  FlatIndex.resize(NumNeighbours*NumSites);
  Hamiltonians.resize(NumNeighbours*NumSites);

  ClusterInfo.resize(NumSites);
  CreateRhombicCluster();

  Defects.resize(NumDefects);
  CreateDefectPositions();
  AddDefectHamiltonia();

  // DebugHamiltonians();

  SSFPoints = {gamma,
               gammapx,    gammapy,    gammapz,
               kx,         ky,         kz,
               kx/2.0,     ky/2.0,     kz/2.0,
               mx,         my,         mz,
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

  Phases.resize(SSFPoints.size());
  VectorXLD v;
  for (uint i=0; i<SSFPoints.size(); i++){
    //create v
    Phases[i] = v;
  }


}

void Triangular::DebugHamiltonians(){
  uint flat1, fnn, fnn2;
  uint wrongham=0;
  double wrongamout=0.0;
  Matrix3LD hamdiff;
  for (uint y1=0; y1<L2; ++y1){
    for (uint x1=0; x1<L1; ++x1){
      for (uint s1=0; s1<NumSublattices;++s1){
        BravaisIndicesToFlat(s1, x1, NumSublattices, y1, L1*NumSublattices, flat1);
        std::cerr << "site1: " << x1 << " " << y1 << " " << s1 << endl;
        std::cerr << "flat site1:" << flat1 << endl;
        std::cerr << "---------entering NN loop --------" << endl;
        for(uint i=0; i<NumNeighbours; i++){
          fnn = FlatIndex[NumNeighbours*flat1+i];
          std::cerr << "flat siteNN: " << fnn << endl;
          /////////////////////output bond Hamiltonian
          // if (FlatIndex[NumNeighbours*flat1+i] == 10){
          //   std::cerr << Hamiltonians[NumNeighbours*flat1+i] << endl;
          // }
          // ///////////////////find other pair's Hamiltonian and output diff
          //
          for (uint z=0; z<NumNeighbours; z++){
            fnn2 = FlatIndex[NumNeighbours*fnn+z];
            if (flat1 == fnn2){
              hamdiff = Hamiltonians[NumNeighbours*flat1+i] - Hamiltonians[NumNeighbours*fnn+z];
              if (hamdiff.norm()> 0){
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

void Triangular::CreateRhombicCluster()
{
  Matrix3LD h1, h2, h3;
  if (ClusterType == 1){
    Translation1 = a1;
    Translation2 = a2;
    h1 = HamInfo.Hx;
    h2 = HamInfo.Hy;
    h3 = HamInfo.Hz;
  } else if (ClusterType == 2){
    Translation1 = a1-a2;
    Translation2 = a1;
    h1 = HamInfo.Hy;
    h2 = HamInfo.Hz;
    h3 = HamInfo.Hx;
  } else if (ClusterType == 3){
    Translation1 = a2;
    Translation2 = a2-a1;
    h1 = HamInfo.Hz;
    h2 = HamInfo.Hx;
    h3 = HamInfo.Hy;
  }

  SublatticeOffset[0] = (Translation1+Translation2)/3.0;

  int higher_x, higher_y, lower_x, lower_y;
  uint flat_index, nn;

  for (int y=0; y<L2; ++y){
    PBCIndices(y, L2, lower_y, higher_y);
    for (int x=0; x<L1; ++x){
      PBCIndices(x, L1, lower_x, higher_x);
      for (uint sub =0; sub<NumSublattices; ++sub){
        BravaisIndicesToFlat(sub, x, NumSublattices, y, L1*NumSublattices, flat_index);
        Position.col(flat_index) = x*Translation1+y*Translation2+SublatticeOffset[sub];

        vector< int>     real_x = {     x+1,        x,      x-1,     x-1,       x,      x+1};
        vector< int>      PBC_x = {higher_x,        x,  lower_x, lower_x,       x, higher_x};
        vector< int>     real_y = {       y,      y+1,      y+1,       y,     y-1,      y-1};
        vector< int>      PBC_y = {       y, higher_y, higher_y,       y, lower_y,  lower_y};
        vector<uint>       subs = {     sub,      sub,      sub,     sub,     sub,      sub};
        vector<Matrix3LD>  hams = {      h2,       h1,       h3,      h2,      h1,       h3};

        for (uint i=0; i<NumNeighbours; i++){
          BravaisIndicesToFlat(subs[i], PBC_x[i], NumSublattices, PBC_y[i], L1*NumSublattices, nn);

          (ClusterInfo[flat_index].NearestNeighbours).push_back(std::make_tuple(real_x[i], real_y[i], subs[i]));
          FlatIndex[NumNeighbours*flat_index+i] = nn;
          Hamiltonians[NumNeighbours*flat_index+i] = hams[i];
        }
      }
    }
  }
}

void Triangular::CreateDefectPositions(){
  // only to be trusted for l1 = l2= 6*n!!
  // only to be trusted for NumDefects = 1,3,9 only!
  uint flat1,flat2,flat3,flat4,flat5,flat6,flat7,flat8,flat9;
  double pos = (double)L1/6.0;
  if (NumDefects >= 1){
    BravaisIndicesToFlat(0, 3*(uint)pos, NumSublattices, 3*(uint)pos, L1*NumSublattices, flat1);
    Defects[0] = flat1;
  }
  if (NumDefects >= 3){
    BravaisIndicesToFlat(0,   (uint)pos, NumSublattices,   (uint)pos, L1*NumSublattices, flat2);
    BravaisIndicesToFlat(0, 5*(uint)pos, NumSublattices, 5*(uint)pos, L1*NumSublattices, flat3);
    Defects[1] = flat2;
    Defects[2] = flat3;
  }
  if (NumDefects == 9){
    BravaisIndicesToFlat(0, 3*(uint)pos, NumSublattices,   (uint)pos, L1*NumSublattices, flat4);
    BravaisIndicesToFlat(0, 5*(uint)pos, NumSublattices,   (uint)pos, L1*NumSublattices, flat5);
    BravaisIndicesToFlat(0,   (uint)pos, NumSublattices, 3*(uint)pos, L1*NumSublattices, flat6);
    BravaisIndicesToFlat(0, 5*(uint)pos, NumSublattices, 3*(uint)pos, L1*NumSublattices, flat7);
    BravaisIndicesToFlat(0,   (uint)pos, NumSublattices, 5*(uint)pos, L1*NumSublattices, flat8);
    BravaisIndicesToFlat(0, 3*(uint)pos, NumSublattices, 5*(uint)pos, L1*NumSublattices, flat9);
    Defects[3] = flat4;
    Defects[4] = flat5;
    Defects[5] = flat6;
    Defects[6] = flat7;
    Defects[7] = flat8;
    Defects[8] = flat9;
  }
}


void Triangular::AddDefectHamiltonia()
{
  Vector2LD rdefect, rsite, rNN, rbond, rseparation;
  auto function = Lorentzian;
  long double strength;
  for (uint d=0; d<Defects.size(); d++){
    rdefect = Position.col(Defects[d]); //+
    for (uint flat_index=0;flat_index<NumSites;flat_index++){
      rsite = Position.col(flat_index);
      for (uint i=0; i<NumNeighbours; i++){
        //I'm using the non-wrapped position indices here.
        //this is to deal with the (periodic) bonds that go beyond the finite size cluster
        auto nninfo = ClusterInfo[flat_index].NearestNeighbours[i];
        rNN = get<0>(nninfo)*Translation1+get<1>(nninfo)*Translation2+SublatticeOffset[get<2>(nninfo)];
        rbond = (rsite+rNN)/2.0;
        rseparation = rbond - rdefect;
        strength = function(rseparation.norm(),DefectLengthScale);
        Hamiltonians[NumNeighbours*flat_index+i](1,1) += +1*DefectOcto*strength;
        Hamiltonians[NumNeighbours*flat_index+i](0,0) += +1*DefectQuad*strength;
        Hamiltonians[NumNeighbours*flat_index+i](2,2) += +1*DefectQuad*strength;
      }
    }
  }
}


void Triangular::CalculateClusterOP()
{
  // Vector3LD SS_ij;
  // Vector2LD r_ij;
  //
  // Vector3LDc SS_q;
  //
  // vector<Vector3LD> clusterssf(SSFPoints.size(),Vector3LD::Zero());
  // for (uint q=0; q<SSFPoints.size(); q++){
  //   SS_q = Vector3LDc::Zero();
  //   for (uint i =0; i<NumSites;i++){
  //     for (uint j =0; j<NumSites;j++){
  //       SS_ij = Cluster[i].cwiseProduct(Cluster[j]);
  //       r_ij = ClusterInfo[i].Position - ClusterInfo[j].Position;
  //       SS_q += std::polar((long double)1.0,-SSFPoints[q].dot(r_ij))*SS_ij;
  //     }
  //   }
  //   clusterssf[q]=SS_q.cwiseAbs();
  // }
  // ClusterSSf = clusterssf;
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
  out << "Magnetic order parameter:\n";
  out << "             Quad             Octo\n";
  for(uint j=0; j <SSFPoints.size(); j++){
    whoa = SSfBar[j]/ns2;
    out << SSFPointsLabel[j] << sqrt(whoa(0) +whoa(2)) << " " << sqrt(whoa(1))  << endl;
  }
}
