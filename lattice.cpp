///  @file     lattice.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the types of Lattices
#include "lattice.hpp"

Lattice::Lattice(const uint& which_lattice,
                       const uint& shape,
                       const uint& l1, const uint& l2, const uint& l3,
                       Hamiltonia& haminfo
                      ):
WhichLattice(which_lattice),
L1(l1), L2(l2), L3(l3),
NumUnitCells(l1*l2*l3),
HamInfo(haminfo),
hField(haminfo.hField)
{
  SpecifyLattice(shape);
  // CheckPositions();
  if ((L1 == 0) || (L2 == 0) || (L3 == 0)){
    std::cerr << "Lattice dimensions should be greater than 1." << endl;
    abort();
  }
  if ((LatticeDimensions != 2) && (LatticeDimensions != 3)){
    std::cerr << "Can only handle 2- or 3-dimensional lattices for now." << endl;
    abort();
  }
  if ((LatticeDimensions == 2) && (L3 != 1)){
    std::cerr << "Two-dimensional lattices require L3 = 1." << endl;
    abort();
  }
  CreateLatticePositions();
  Reciprocal recip(PrimitiveCell, WhichLattice);
  QPoints = recip.QPoints;
  QPointsLabel = recip.QPointsLabel;
  CalculatePhases();
  //
  FindNearestNeighbours();

  // CheckNNBonds();
  Cluster.resize(3, NumSites);
  // CheckPositions();
}

void Lattice::SpecifyLattice(const uint& shape)
{
  int volumefactor;
  Matrix3XLD s_matrix, nn_matrix;
  Matrix3LD t_matrix, xyztoabc, tri_spatial_basis, permute_matrix_1,permute_matrix_2;
  tri_spatial_basis <<  + 1.0, +          0.5, 0.0,
                        + 0.0, +sqrt(3.0)*0.5, 0.0,
                          0.0,            0.0, 1.0;
  permute_matrix_1 <<  0.0, 1.0, 0.0,
                      -1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0,
  permute_matrix_2 <<  0.0, 1.0, 0.0,
                      -1.0, 0.0, 0.0,
                       0.0, 0.0, 1.0,
  xyztoabc << A_Dir, B_Dir, C_Dir;
  switch (WhichLattice){ //0 for triangular, 1 for honeycomb, 2 for FCC
    case 0: //triangular
      LatticeDimensions  = 2;
      CoordinationNumber = 6;
      NumSublattices     = 1;
      nn_matrix.resize(3,CoordinationNumber);
      //            z  y   x   z   y   x
      nn_matrix <<  1, 0, -1, -1,  0,  1,
                    0, 1,  1,  0, -1, -1,
                    0, 0,  0,  0,  0,  0;
      NNBondType.resize(CoordinationNumber);
      NNBondType = {2, 1, 0, 2, 1, 0};
      LocalOctahedra = permute_matrix_2*xyztoabc.transpose();
      PrimitiveCell = tri_spatial_basis;
      switch (shape){ //0 for rhom (primitive), 1 for rect, 2 for sqrt
        case 0:
          t_matrix = Matrix3LD::Identity();
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/3.0,
                       1/3.0,
                         0.0;
          break;
        case 1:
          t_matrix << 1, -1, 0,
                      0,  2, 0,
                      0,  0, 1;
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/2.0, 1/1.0,
                       1/6.0, 4/6.0,
                         0.0,   0.0;
          break;
        case 2:
          t_matrix << 1, -2, 0,
                      1,  1, 0,
                      0,  0, 1;
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/6.0, 5/6.0, 3/6.0,
                       1/6.0, 3/6.0, 5/6.0,
                         0.0,   0.0,   0.0;
          break;
      }
      ConventionalCell = PrimitiveCell * t_matrix;
      SublatticeVectors = ConventionalCell * s_matrix;
      NNVectors = PrimitiveCell * nn_matrix;
      break;
    case 1: //honeycomb
      LatticeDimensions  = 2;
      CoordinationNumber = 3;
      NumSublattices     = 2;
      nn_matrix.resize(3,2*CoordinationNumber); //due to honeycomb sublattice
      //            x        z         y         x         z          y
      nn_matrix <<  1.0/3.0, 1.0/3.0, -2.0/3.0, -1.0/3.0, -1.0/3.0,   2.0/3.0,
                    1.0/3.0,-2.0/3.0,  1.0/3.0, -1.0/3.0,  2.0/3.0,  -1.0/3.0,
                    0.0/1.0, 0.0/1.0,  0.0/1.0,  0.0/1.0,  0.0/1.0,   0.0/1.0;
      NNBondType.resize(2*CoordinationNumber);
      NNBondType = {0, 2, 1, 0, 2, 1};
      LocalOctahedra = xyztoabc.transpose();
      //lattice constant is equal to 1/sqrt(3)
      PrimitiveCell = tri_spatial_basis;
      switch (shape){ //0 for rhom (primitive), 1 for rect, 2 for sqrt
        case 0:
          t_matrix = Matrix3LD::Identity();
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/3.0, 2/3.0,
                       1/3.0, 2/3.0,
                         0.0,   0.0;
          break;
        case 1:
          t_matrix << 1, -1, 0,
                      0,  2, 0,
                      0,  0, 1;
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/2.0, 1/1.0, 1/1.0, 1/2.0,
                       1/6.0, 1/3.0, 2/3.0, 5/6.0,
                         0.0,   0.0,   0.0,   0.0;
          break;
        case 2:
          t_matrix << 1, -2, 0,
                      1,  1, 0,
                      0,  0, 1;
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<  1/6.0, 3/6.0, 5/6.0, 5/6.0, 3/6.0, 1/6.0,
                       1/6.0, 1/6.0, 3/6.0, 5/6.0, 5/6.0, 3/6.0,
                         0.0,   0.0,   0.0,   0.0,   0.0,   0.0;
          break;
      }
      ConventionalCell = PrimitiveCell * t_matrix;
      SublatticeVectors = ConventionalCell * s_matrix;
      NNVectors = PrimitiveCell * nn_matrix;
      break;
    case 2: //FCC
      LatticeDimensions  =  3;
      CoordinationNumber = 12;
      NumSublattices     =  1;
      nn_matrix.resize(3,CoordinationNumber);
                  //x        y        z         x         y         z         x         y         z         x         y         z
      nn_matrix <<  1.0/1.0, 0.0/1.0, 0.0/1.0, -1.0/1.0,  0.0/1.0,  0.0/1.0,  0.0/1.0, +1.0/1.0, -1.0/1.0,  0.0/1.0, -1.0/1.0, +1.0/1.0,
                    0.0/1.0, 1.0/1.0, 0.0/1.0, -0.0/1.0, -1.0/1.0,  0.0/1.0, -1.0/1.0,  0.0/1.0, +1.0/1.0, +1.0/1.0,  0.0/1.0, -1.0/1.0,
                    0.0/1.0, 0.0/1.0, 1.0/1.0, -0.0/1.0,  0.0/1.0, -1.0/1.0, +1.0/1.0, -1.0/1.0,  0.0/1.0, -1.0/1.0, +1.0/1.0,  0.0/1.0;
      NNBondType.resize(CoordinationNumber);
      NNBondType = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
      LocalOctahedra = permute_matrix_1 * xyztoabc.transpose();
      PrimitiveCell = LocalOctahedra * (Matrix3LD::Ones() - Matrix3LD::Identity())/sqrt(2.0);
      switch (shape){
        case 0: //0 for FCC primitive (xyz), 1 for FCC conventional (ABC)
          t_matrix << Matrix3LD::Identity();
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<   0.0,
                        0.0,
                        0.0;
          break;
        case 1:
          t_matrix << -1,  1,  1,
                       1, -1,  1,
                       1,  1, -1;
          volumefactor = lrint(t_matrix.determinant());
          NumTotalSublattices = NumSublattices*volumefactor;
          s_matrix.resize(3,NumTotalSublattices);
          s_matrix <<   0.0, 0.0, 0.5, 0.5,
                        0.0, 0.5, 0.0, 0.5,
                        0.0, 0.5, 0.5, 0.0;
          break;
      }
      ConventionalCell = PrimitiveCell * t_matrix;
      SublatticeVectors = ConventionalCell * s_matrix;
      NNVectors = PrimitiveCell * nn_matrix;
      break;
  }
  NumSites = NumUnitCells*NumTotalSublattices;
}

void Lattice::CreateLatticePositions()
{
  Positions.resize(3,NumSites);
  uint flat_index;
  for (int z=0; z<L3; ++z){
    for (int y=0; y<L2; ++y){
      for (int x=0; x<L1; ++x){
        for (uint sub =0; sub<NumTotalSublattices; ++sub){
          ThreeDBravaisIndicesToFlat(sub, x, NumTotalSublattices,
                                          y, NumTotalSublattices*L1,
                                          z, NumTotalSublattices*L1*L2,
                                          flat_index);
          Positions.col(flat_index) = x*ConventionalCell.col(0)
                                     +y*ConventionalCell.col(1)
                                     +z*ConventionalCell.col(2)
                                     +SublatticeVectors.col(sub);
                                     // cout<<"wow"<<endl;
            // cout << x << " " << y << " " << z << " " << sub << " " << flat_index
            //      << " " << Positions.col(flat_index).transpose() << endl;
        }
      }
    }
  }
}

void Lattice::FindNearestNeighbours()
//function that uses the positions of each site to
//1) find the nearest neighbours and
//2) associate a bond Hamiltonian to each pair of nearest neighbours
//it does so by checking if ri - rj = d, where
//d is one of the columns of NNVectors. however, we need to
//include the case where ri and rj cross the edge of the periodic cluster.
//so we must check if ri - rj = d up to translations of the whole cluster, ie.
// (ri - rj) - d = (T1, T2, T3) . (n1/L1, n2/L2, n3/L3) AND
// n1/L1, n2/L2, n3/L3 are all INTEGERS
//it then checks if the nearest neighbours were found by checking that, for each site,
//the number of nearest neighbours = the coordination number of the crystal lattice
{
  NNSiteInfo.resize(NumSites);

  Vector3LD ri, rj, R_sep;

  Matrix3LD inv_translation = ConventionalCell.inverse();
  Vector3LD solpluss, solp, ls;
  ls << 1/(double)L1, 1/(double)L2, 1/(double)L3;
  bool ispbond;
  long double eps = 1e-6;

  SiteInfo site_info;
  vector<Matrix3LD> hams = {HamInfo.Hx, HamInfo.Hy, HamInfo.Hz};
  // cout << hams[0]<< endl;
  for (uint j=0; j<NumSites; ++j){                      ///for eacb site j...
    rj = Positions.col(j);
    site_info = {};
    for (uint i=0; i<NumSites; ++i){                    ///...check if site i...
      // cout << "(i, j) = " << i << "," <<j << endl;
      ri = Positions.col(i);
      R_sep = ri - rj; //points from site j to site i
      for (uint k=0;k<NNVectors.cols();k++){          ///...is a NN of type BondType[k]...
        solpluss = inv_translation * (R_sep - NNVectors.col(k));
        solp = (solpluss.array()*ls.array()).matrix();
        ispbond = IsThisFloatAnIntegerWithinEpsilon(solp(0),eps) &&
                  IsThisFloatAnIntegerWithinEpsilon(solp(1),eps) &&
                  IsThisFloatAnIntegerWithinEpsilon(solp(2),eps);
        if (ispbond){
          // cout << "direction? " << k << " type? " << NNBondType[k] << endl;
          site_info.NNInfo.push_back(std::make_tuple(hams[NNBondType[k]], NNBondType[k], i, k));
        }
      }
    }
    //first sign something is wrong: each site has a different number of nearest neighbours
    //so, let us assert here that the number of nearest neighbours found
    //is equal to the coordination number.
    if (site_info.NNInfo.size() != CoordinationNumber){
      std::cerr << "Number of nearest neighbours found is not equal to the coordination number." << endl;
      std::cerr << "Please check your lattice creation." << endl;
      abort();
    } else{
      NNSiteInfo[j] = site_info;
    }
  }
}

void Lattice::CheckPositions()
{
  std::cerr << ConventionalCell << endl;
  std::cerr << SublatticeVectors << endl;
  std::cerr << NNVectors << endl;
}

void Lattice::CheckNNBonds()
{
  for (uint i=0; i<NumSites; i++){
    std::cerr << "site: " << i << endl;
    SiteInfo site_info =  NNSiteInfo[i];
    std::cerr << "NN list: " << endl;
    for (auto& j : site_info.NNInfo){
      std::cerr << "NN " << get<2>(j) << " with type "<< get<1>(j) << " and H = "<< endl;
      std::cerr << get<0>(j)<< endl;
    }
    std::cerr << "______________________" << endl;
  }
}

void Lattice::CheckHamiltonians()
{
  uint wrongham=0;
  long double wrongamount=0.0;
  Matrix3LD hamdiff;
  SiteInfo siteinfo1, siteinfo2;
  for (uint i=0; i<NumSites; i++){
    std::cerr << "site1: " << i << endl;
    siteinfo1 = NNSiteInfo[i];
    std::cerr << "------entering NN loop------" << endl;
    for (uint s1=0; s1 < CoordinationNumber; s1++){
      auto nninfo1 = siteinfo1.NNInfo[s1];
      ////////////////////////////////////////////////////////
      // output the neighbours and bond types
      std::cerr << "site1 NN: " << get<2>(nninfo1) << " of type " << get<1>(nninfo1) << endl;
      std::cerr << get<0>(nninfo1) << endl;
      ////////////////////////////////////////////////////////
      // find other pair's Hamiltonian and output diff
      // siteinfo2 = NNSiteInfo[get<2>(nninfo1)];
      // for (uint s2=0; s2 < CoordinationNumber; s2++){
      //   auto nninfo2 = siteinfo2.NNInfo[s2];
      //   if (get<2>(nninfo2) == i){
      //     hamdiff = get<0>(nninfo1) - get<0>(nninfo2);
      //     if (hamdiff.norm()> 0){
      //       wrongham++;
      //       wrongamount+=hamdiff.norm();
      //     }
      //   }
      // }
    }
    std::cerr << "------finished NN loop------" << endl;
  }
  std::cerr << wrongham << " wrong bonds, total wrong is " << wrongamount << endl;
}

void Lattice::PrintLatticeParameters(std::ostream &out)
{
  out << "-------------------------Simulation Parameters-------------------------\n";
  out << "Which lattice?\n";
  out << WhichLattice << endl;
  out << "Linear dimensions (L1, L2, L3)\n";
  out << L1 << " " << L2 << " " << L3 << endl;
  out << "Translation vectors (col corresponds to T_1, T_2, T_3)\n";
  out << ConventionalCell << endl;
  out << "Sublattice vectors (col corresponds to s_1, s_2,...s_Ns)\n";
  out << SublatticeVectors << endl;
}
void Lattice::PrintHamiltonianParameters(std::ostream &out)
{
  out << HamInfo.ParameterOutput;
  out << "Spin basis (matrix from spin basis to global basis)\n";
  out << LocalOctahedra * HamInfo.SpinBasis << endl;
  if (Defects.size() > 0){
    out << "Defect quad/defect octo/defect length scale/number of defects\n";
    out << DefectQuad << " " << DefectOcto << " " << DefectLengthScale << " " << Defects.size() << "\n";
  }
}

void Lattice::CalculatePhases()
{
  Phases.resize(QPoints.cols(),NumSites);
  std::complex<long double> im(0.0,1.0);
  Phases << exp((-im * QPoints.transpose() * Positions).array()).matrix();
}

void Lattice::CalculateOP()
{
  ArrayXXLD sxsx = (Cluster.row(0).transpose() * Cluster.row(0)).array();
  ArrayXXLD sysy = (Cluster.row(1).transpose() * Cluster.row(1)).array();
  ArrayXXLD szsz = (Cluster.row(2).transpose() * Cluster.row(2)).array();

  SSFComponents.resize(3,Phases.rows());

  for (uint i=0; i<Phases.rows(); i++){
    auto sxx = ((Phases.row(i).adjoint() * Phases.row(i)).array() * sxsx).sum();
    auto syy = ((Phases.row(i).adjoint() * Phases.row(i)).array() * sysy).sum();
    auto szz = ((Phases.row(i).adjoint() * Phases.row(i)).array() * szsz).sum();
    SSFComponents.col(i) << sxx.real(), syy.real(), szz.real();
  }
}
void Lattice::PrintOP(std::ostream &out)
{
  long double ns2 = pow(NumSites,2);
  out << "-------------------------------Final configuration observables--------------------------------\n";
  out << "beware: only trust values for q-points accessible by the cluster.\n";
  out << "             sxx             syy             szz\n";
  for(uint j=0; j<QPoints.cols(); j++){
    out << QPointsLabel[j] << SSFComponents.col(j)(0)/ns2
                    << " " << SSFComponents.col(j)(1)/ns2
                    << " " << SSFComponents.col(j)(2)/ns2  << endl;
  }
}
void Lattice::PrintThermalObservables(std::ostream &out)
{
  long double ns = NumSites;
  long double ns2 = pow(NumSites,2);
  out << std::setprecision(14);
  out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  out << "Energy cumulants: <E>, <E2> \n";
  out << EBar/ns << " " << E2Bar/ns2 << "\n";
  // out << "Magnetic order parameter:\n";
  // out << "             Quad             Octo\n";
  // for(uint j=0; j <SSFPoints.size(); j++){
  //   whoa = SSfBar[j]/ns2;
  //   out << SSFPointsLabel[j] << sqrt(whoa(0) +whoa(2)) << " " << sqrt(whoa(1))  << endl;
  // }
}

void Lattice::CreateDefectPositions()
{
  Defects.resize(1);

  uint flat_index;
  uint x = (int)(L1/2.0);
  uint y = (int)(L2/2.0);
  uint z = (int)(L3/2.0);
  uint sub = 0;

  // cout << x << " "  << y << " " <<  z << endl;

  ThreeDBravaisIndicesToFlat(sub, x, NumTotalSublattices,
                                  y, NumTotalSublattices*L1,
                                  z, NumTotalSublattices*L1*L2,
                                  flat_index);

  Defects[0] = flat_index;
}

void Lattice::AddDefectHamiltonia(const long double& defect_quad       ,
                                     const long double& defect_octo       ,
                                     const long double& defect_lengthscale)
{
  DefectQuad        = defect_quad       ;
  DefectOcto        = defect_octo       ;
  DefectLengthScale = defect_lengthscale;
  Vector3LD rdefect, rsite, rNN, rbond, rseparation;
  auto function = Lorentzian;
  long double strength;
  for (uint d=0; d<Defects.size(); d++){
    rdefect = Positions.col(Defects[d]); //+
    for (uint site=0;site<NumSites;site++){
      rsite = Positions.col(site);
      for (uint nn=0; nn<CoordinationNumber; nn++){
        // rNN = Positions.col(get<2>(NNSiteInfo[site].NNInfo[nn]));
        rNN = rsite+NNVectors.col(get<3>(NNSiteInfo[site].NNInfo[nn]));
        rbond = (rsite+rNN)/2.0;
        rseparation = rbond - rdefect;
        strength = function(rseparation.norm(),DefectLengthScale);
        get<0>(NNSiteInfo[site].NNInfo[nn])(0,0) += +1*DefectQuad*strength;
        get<0>(NNSiteInfo[site].NNInfo[nn])(1,1) += +1*DefectOcto*strength;
        get<0>(NNSiteInfo[site].NNInfo[nn])(2,2) += +1*DefectQuad*strength;
      }
    }
  }
}

Reciprocal::Reciprocal(Matrix3LD& primitive_cell, const uint& which_lattice)
{
  CreateReciprocalVectors(primitive_cell, CrystalReciprocalVectors);
  CreateSymmetryPoints(which_lattice);
}

void Reciprocal::CreateReciprocalVectors(const Matrix3LD& translation_vectors,
                                               Matrix3LD& recip_vectors)
{
  Matrix3LD inv_translation = (translation_vectors.transpose()).inverse();
  recip_vectors = 2*pi*inv_translation*Matrix3LD::Identity();
}

void Reciprocal::CreateSymmetryPoints(const uint& which_lattice)
{
  Vector3LD b1, b2, b3;
  b1 = CrystalReciprocalVectors.col(0);
  b2 = CrystalReciprocalVectors.col(1);
  b3 = CrystalReciprocalVectors.col(2);
  if ( (which_lattice==0) || (which_lattice==1) ){
    Vector3LD bB1 =  b1+b2;
    Vector3LD bB2 = -b1;

    Vector3LD   gamma = (0.0*bB1+0.0*bB2)/1.0;
    Vector3LD gammapx = (1.0*bB1+0.0*bB2)/1.0;
    Vector3LD gammapy = (0.0*bB1+1.0*bB2)/1.0;
    Vector3LD gammapz = (1.0*bB1+1.0*bB2)/1.0;
    Vector3LD      kx = (1.0*bB1+2.0*bB2)/3.0;
    Vector3LD      ky = (2.0*bB1+1.0*bB2)/3.0;
    Vector3LD      kz = (1.0*bB1-1.0*bB2)/3.0;
    Vector3LD      mx = (1.0*bB1+0.0*bB2)/2.0;
    Vector3LD      my = (0.0*bB1+1.0*bB2)/2.0;
    Vector3LD      mz = (1.0*bB1+1.0*bB2)/2.0;

    int num_q = 16;
    QPoints.resize(3,num_q);
    QPoints << gamma,
               gammapx,  gammapy,  gammapz,
               kx,       ky,       kz,
               kx/2.0,   ky/2.0,   kz/2.0,
               mx,       my,       mz,
               2*mx/3.0, 2*my/3.0, 2*mz/3.0;
    QPointsLabel.resize(num_q);
    QPointsLabel = {"Gamma   ",
                    "GammaPx ", "GammaPy ", "GammaPz ",
                    "Kx      ", "Ky      ", "Kz      ",
                    "Kx/2    ", "Ky/2    ", "Kz/2    ",
                    "Mx      ", "My      ", "Mz      ",
                    "2Mx/3   ", "2My/3   ", "2Mz/3   "
                  };
  }else if (which_lattice==2){
    Vector3LD   gamma = (0.0*b1+0.0*b2+0.0*b3)/1.0;
    Vector3LD       l = (1.0*b1+1.0*b2+1.0*b3)/2.0;
    Vector3LD       k = (3.0*b1+3.0*b2+6.0*b3)/8.0;
    Vector3LD       u = (5.0*b1+2.0*b2+5.0*b3)/8.0;
    Vector3LD       w = (2.0*b1+1.0*b2+3.0*b3)/4.0;
    Vector3LD       x = (1.0*b1+0.0*b2+1.0*b3)/2.0;

    int num_q = 6;
    QPoints.resize(3,num_q);
    QPoints << gamma,
               l,
               k,
               u,
               w,
               x;
    QPointsLabel.resize(num_q);
    QPointsLabel = {"Gamma   ",
                    "L       ",
                    "K       ",
                    "U       ",
                    "W       ",
                    "X       "
                  };
  }
}
