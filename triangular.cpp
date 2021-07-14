///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#include "triangular.hpp"

TriangularLattice::TriangularLattice(const uint& l1, const uint& l2, const uint& num_defects,
                                     const long double& jtau, const long double& lambda,
                                     const long double& ising_y,
                                     const long double& defect,
                                     const long double& h,
                                     Vector3LD& hdir):
L1(l1), L2(l2), NumSites(l1*l2), NumDefects(num_defects),
JTau(jtau), Lambda(lambda), IsingY(ising_y), Defect(defect),
HField(h), HDirection(hdir)
{
  // defining the bond Hamiltonians
  FixMPHamiltonians();

  //changing size of Cluster to L1*L2
  Cluster.resize(L1);
  for(uint i = 0; i < L1; ++i){
      Cluster[i].resize(L2);
  }
  CreateClusterPBC();

  Defects.resize(NumDefects);
  for(uint i = 0; i < NumDefects; ++i){
      Defects[i].resize(2);
  }
  CreateDefectPositions();
  Hdefect1 << 0,      0, 0,
             0, Defect, 0,
             0,      0, 0;
  Hdefect1=JTau*Hdefect1;
  Hdefect2 = 0.0*Hdefect1;
  AddDefectHamiltonia();

  // cout << Hdefect1 << endl;
  // cout << Hdefect2 << endl;

  // InitializeRandomSpins();
  InitializeFMSpins(pi/2.0,pi/2.0);

  CreateStripySignMatrices();
  CalculateClusterEnergyandOP();

  std::uniform_int_distribution<uint> l1d(0, L1-1);
  std::uniform_int_distribution<uint> l2d(0, L2-1);
  L1Dist = l1d;
  L2Dist = l2d;

  overrelaxMCratio = 10;
}

void TriangularLattice::CreateClusterPBC()
// T1 = a1-a2, T2 = a1
{
  Spin some_spin;
  int higher_x, higher_y, lower_x, lower_y;
  Matrix3LD hamiltonian = Matrix3LD::Zero();
  for (int y=0; y<L2; ++y)
  {
    for (int x=0; x<L1; ++x)
    {
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
      Cluster[x][y] = {
                        {std::make_tuple(higher_x,        y, Hz),
                         std::make_tuple(       x, higher_y, Hy),
                         std::make_tuple( lower_x, higher_y, Hx),
                         std::make_tuple( lower_x,        y, Hz),
                         std::make_tuple(       x,  lower_y, Hy),
                         std::make_tuple(higher_x,  lower_y, Hx)},
                         some_spin
                      };
    }
  }
}

void TriangularLattice::CreateDefectPositions(){
  // only to be trusted for l1 = l2= 6*n!!
  // only to be trusted for NumDefects = 1,3,9 only!

  double pos = (double)L1/6.0;
  if (NumDefects >= 1){
    Defects[0] = {3*(uint)pos,3*(uint)pos};
  }
  if (NumDefects >= 3){
    Defects[1] = {(uint)pos,(uint)pos};
    Defects[2] = {5*(uint)pos,5*(uint)pos};
  }
  if (NumDefects == 9){
    Defects[3] = {3*(uint)pos,(uint)pos};
    Defects[4] = {5*(uint)pos,(uint)pos};
    Defects[5] = {(uint)pos,3*(uint)pos};
    Defects[6] = {5*(uint)pos,3*(uint)pos};
    Defects[7] = {(uint)pos,5*(uint)pos};
    Defects[8] = {3*(uint)pos,5*(uint)pos};
  }
  // uint i = 0;
  // for (uint y=0; y<L2; ++y){
  //   for (uint x=0; x<L1; ++x){
  //     i = y+x*L1;
  //     // Defects[i] = {x};
  //   }
  // }

  // cout << "...Defects list..." << endl;
  // for (auto i:Defects){
  //   cout << i[0] << " " << i[1] << endl;
  // }
}

bool TriangularLattice::CheckIfPoisoned(uint lx, uint ly){
  // only to be trusted for l1 = l2= 6*n!!
  // only to be trusted for NumDefects = 1,3,9 only!
  uint mult_six = (uint)((double)L1/6.0);
  //the 6 above comes from the location of defects in CreateDefectPositions

  if ((lx%mult_six!=0) or (ly%mult_six!=0)) {
    return false;
  }
  if ( (lx%(2*mult_six)==0) or (ly%(2*mult_six)==0) ) {
    return false;
  }
  // cout << lx << " " << ly<< endl;
  for (auto defectposition:Defects){
    if ((defectposition[0]==lx) and (defectposition[1]==ly)) {
      // cout << lx << " " << ly << endl;
      return true;
    }
  }
  return false;
}

void TriangularLattice::AddDefectHamiltonia()
{
  bool poisoned_site, poisoned_neighbour;
  Site current_site;

  for (int y=0; y<L2; ++y){
    for (int x=0; x<L1; ++x){
      poisoned_site = CheckIfPoisoned(x, y);
      for (uint i=0; i<6; ++i){
        auto [nn_x, nn_y, old_hamiltonian] = Cluster[x][y].NearestNeighbours[i];

        if (poisoned_site == true){
          std::get<2>(Cluster[x][y].NearestNeighbours[i]) += Hdefect1;
        } else {
          poisoned_neighbour = CheckIfPoisoned(nn_x,nn_y);
          if (poisoned_neighbour == true){
            std::get<2>(Cluster[x][y].NearestNeighbours[i]) += Hdefect1;
            uint ibefore = (i-1)%6;
            uint iafter = (i+1)%6;
            std::get<2>(Cluster[x][y].NearestNeighbours[ibefore]) += 0.5*Hdefect2;
            std::get<2>(Cluster[x][y].NearestNeighbours[iafter]) += 0.5*Hdefect2;
          }
        }

      }

    }
  }

}

Matrix3LD TriangularLattice::ReturnMPHamiltonian(const long double& angle)
//Eqn 17 of arXiv:2105.09334v1
{
  long double turn_jtau_off = 1;
  Matrix3LD matrix1, matrix2, matrix3;
  //bond-independent: lambda term + Ising in y direction
  matrix1  <<  -Lambda,      0,       0,
                     0, Lambda,       0,
                     0,      0, -Lambda;

  matrix2  <<             0,      0, 0,
                          0, IsingY, 0,
                          0,      0, 0;

  //bond-dependent: tau term
  long double c = cos(angle);
  long double s = sin(angle);
  matrix3  <<  1-c, 0,  -s,
                 0, 0,   0,
                -s, 0, 1+c;

  return matrix1/2.0 + matrix2 + turn_jtau_off*matrix3/2.0;
}

void TriangularLattice::FixMPHamiltonians()
{
  long double pz = 0;
  long double px = 2*pi/3.0;
  long double py = 4*pi/3.0;

  Hz = ReturnMPHamiltonian(pz)*JTau;
  Hx = ReturnMPHamiltonian(px)*JTau;
  Hy = ReturnMPHamiltonian(py)*JTau;
}

void TriangularLattice::InitializeFMSpins(const long double& theta, const long double& phi)
{
  Spin spin(theta, phi);
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      Cluster[x][y].OnsiteSpin = spin;
    }
  }
}

void TriangularLattice::InitializeRandomSpins()
{
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      SpherePointPicker(Cluster[x][y].OnsiteSpin);
    }
  }
}

void TriangularLattice::MolecularField(const Site& site, Vector3LD& molec)
{
  Vector3LD v = Vector3LD::Zero();
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto j : site.NearestNeighbours){
    auto [nn_x, nn_y, ham] = j;
    spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
    v += -spin_j.VectorXYZ.transpose()*ham;
  }
  molec = v.transpose()+ JTau*HField*HDirection.transpose();
}

void TriangularLattice::CalculateLocalEnergy(const Site& site, long double& energy)
{
  Vector3LD molec;
  MolecularField(site, molec);
  energy = -site.OnsiteSpin.VectorXYZ.dot(molec+ JTau*HField*HDirection);
  // cout << energy<< endl;
}

void TriangularLattice::CalculateClusterEnergy(){
  long double e=0;
  long double local_energy;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      //calculate energy
      CalculateLocalEnergy(Cluster[n1][n2], local_energy);
      e += local_energy;
    }
  }
  ClusterEnergy = e/2.0;
}

void TriangularLattice::CalculateClusterEnergyandOP()
{
  long double e=0;
  long double local_energy;

  Vector3LD spin, fmop;
  Eigen::Matrix<long double, 3, 2> stripyopmatrix;
  fmop << 0,0,0;
  stripyopmatrix << 0, 0,
                    0, 0,
                    0, 0;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      // local_energy=0;
      //calculate energy
      CalculateLocalEnergy(Cluster[n1][n2], local_energy);
      e += local_energy;
      //calculate three stripy OPs

      spin = Cluster[n1][n2].OnsiteSpin.VectorXYZ;
      fmop += spin;
      //
      stripyopmatrix(0,0) += StripySignsX(n1,n2)*spin(0); stripyopmatrix(0,1) += StripySignsX(n1,n2)*spin(2);
      stripyopmatrix(1,0) += StripySignsY(n1,n2)*spin(0); stripyopmatrix(1,1) += StripySignsY(n1,n2)*spin(2);
      stripyopmatrix(2,0) += StripySignsZ(n1,n2)*spin(0); stripyopmatrix(2,1) += StripySignsZ(n1,n2)*spin(2);
      }
    }

    ClusterEnergy = e/2.0;

    ClusterFMOP = fmop;
    ClusterStripyOPMatrix = stripyopmatrix;

    SelectStripyOP();
}

void TriangularLattice::OverrelaxationSweep(){
  Site *chosen_site_ptr;
  Vector3LD old_spin_vec,new_spin_vec;
  Vector3LD molec_field, normalized_field,spindiff;

  for (uint n1 =0; n1<L1; n1++){
    for (uint n2 =0; n2<L2; n2++){
      chosen_site_ptr = &Cluster[n1][n2];
      old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

      MolecularField(*chosen_site_ptr, molec_field);
      normalized_field = molec_field.normalized();

      new_spin_vec = -old_spin_vec + 2.0*normalized_field.dot(old_spin_vec)*normalized_field;
      Spin new_spin(new_spin_vec);
      chosen_site_ptr->OnsiteSpin = new_spin;

      //no need to update energy dynamically since these sweeps perserve energy
      //updating OP dynamically
      spindiff = new_spin_vec - old_spin_vec;
      ClusterFMOP+=spindiff;
      ClusterStripyOPMatrix(0,0) += StripySignsX(n1,n2)*spindiff(0); ClusterStripyOPMatrix(0,1) += StripySignsX(n1,n2)*spindiff(2);
      ClusterStripyOPMatrix(1,0) += StripySignsY(n1,n2)*spindiff(0); ClusterStripyOPMatrix(1,1) += StripySignsY(n1,n2)*spindiff(2);
      ClusterStripyOPMatrix(2,0) += StripySignsZ(n1,n2)*spindiff(0); ClusterStripyOPMatrix(2,1) += StripySignsZ(n1,n2)*spindiff(2);
    }
  }
}

void TriangularLattice::MetropolisSweep(const double& temperature, uint& accept)
{
  uint n1, n2;
  Spin old_spin_at_chosen_site;
  long double old_local_energy, energydiff, new_local_energy;
  double r, ppp;
  Site *chosen_site_ptr;
  Vector3LD spindiff;
  uint flip = 0;

  accept=0;

  while (flip < NumSites){
    n1 = L1Dist(MyRandom::RNG);
    n2 = L2Dist(MyRandom::RNG);

    chosen_site_ptr = &Cluster[n1][n2];

    CalculateLocalEnergy(*chosen_site_ptr, old_local_energy);
    old_spin_at_chosen_site = chosen_site_ptr->OnsiteSpin;
    //selection of angle, currently uniform update
    SpherePointPicker(chosen_site_ptr->OnsiteSpin);

    CalculateLocalEnergy(*chosen_site_ptr, new_local_energy);

    energydiff = new_local_energy-old_local_energy;
    r = MyRandom::unit_interval(MyRandom::RNG);
    ppp = std::min(exp(-energydiff/temperature),1.0);

    if (r < ppp) {
      accept++;
      //update energy dynamically
      ClusterEnergy+=energydiff;
      //update OP dynamically
      spindiff = (chosen_site_ptr->OnsiteSpin).VectorXYZ - old_spin_at_chosen_site.VectorXYZ;
      ClusterFMOP+=spindiff;
      ClusterStripyOPMatrix(0,0) += StripySignsX(n1,n2)*spindiff(0); ClusterStripyOPMatrix(0,1) += StripySignsX(n1,n2)*spindiff(2);
      ClusterStripyOPMatrix(1,0) += StripySignsY(n1,n2)*spindiff(0); ClusterStripyOPMatrix(1,1) += StripySignsY(n1,n2)*spindiff(2);
      ClusterStripyOPMatrix(2,0) += StripySignsZ(n1,n2)*spindiff(0); ClusterStripyOPMatrix(2,1) += StripySignsZ(n1,n2)*spindiff(2);
    } else {
      chosen_site_ptr->OnsiteSpin = old_spin_at_chosen_site;
    }
    ++flip;
  }

}

void TriangularLattice::DoTheSweeps(double& temp, uint& accept){
  for (uint i=0;i<overrelaxMCratio;i++){
    OverrelaxationSweep();
  }
  MetropolisSweep(temp,accept);
  cout << ClusterEnergy/(double)NumSites << endl;
}


void TriangularLattice::ThermalizeConfiguration(double& temp, const uint& max_sweeps)
{
    // cout << temp << " " << temp << endl;

    uint sweep = 0;

    uint accept;
    uint Metroaccept = 0;
    // uint overrelaxMCratio = 1;
    while (sweep < max_sweeps){
      DoTheSweeps(temp, accept);
      Metroaccept += accept;
      ++sweep;
    }
    // std::cerr << temp << " " << (double)Metroaccept/(double)max_sweeps/(double)NumSites << endl;
}

void TriangularLattice::SimulatedAnnealing(const uint& max_sweeps,
                                          double& initial_T, double& final_T, double& rate)
{
  double scale = rate;
  double temp_T = initial_T;
  while(temp_T >= final_T){
    // std::cout << "temp " << temp_T << endl;
    ThermalizeConfiguration(temp_T,max_sweeps);

    // CalculateClusterEnergyandOP();
    // cout << temp_T << " " << ClusterEnergy/NumSites << endl;

    temp_T = scale*temp_T;
  }
  FinalT = temp_T;
}

void TriangularLattice::DeterministicSweeps(const uint& max_sweeps)
{
  uint sweep = 0;
  uint align;
  uint uc_x, uc_y;
  Site *chosen_site_ptr;
  Spin new_spin;

  Vector3LD old_spin_vec;
  Vector3LD molec_field;

  // long double x, norm;
  // long double eps = pow(10,-20);

  while (sweep < max_sweeps){
    // x = 0;
    align = 0;
    while (align < NumSites){
        uc_x = L1Dist(MyRandom::RNG);
        uc_y = L2Dist(MyRandom::RNG);

        chosen_site_ptr = &Cluster[uc_x][uc_y];
        old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

        MolecularField(*chosen_site_ptr, molec_field);
        Vector3LD mf_t = molec_field;
        Spin new_spin(mf_t);
        chosen_site_ptr->OnsiteSpin = new_spin;
        align++;
    }
    sweep++;
  }
  ActualDetSweeps = sweep;
  CalculateClusterEnergyandOP(); //updates energy and OP statically
  // cout << "final: " << std::setprecision(14) << ClusterEnergy/NumSites << endl;
}

void TriangularLattice::PrintConfiguration(std::ostream &out)
{
  out << "--------------------------------Final Configuration--------------------------------\n";
  out << "Energy per site\n";
  out << std::setprecision(14) << ClusterEnergy/NumSites << "\n";
  out << "Spin configuration\n";
  for (uint y=0; y<L2; ++y){
    for (uint x=0; x<L1; ++x){
      out << std::setprecision(14) << x << " " << y << " " << 0 << " " << Cluster[x][y].OnsiteSpin.VectorXYZ.transpose() << "\n";
    }
  }
  SelectStripyOP();
  out << "-------------------------------Final configuration observables--------------------------------\n";
  out << "Order parameters. R1: (FM_x, FM_y, FM_z), R2: (stripy_x, stripy_z),  R3: (stripy_x, FM_y, stripy_z)\n";
  out << std::setprecision(14) << ClusterFMOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterStripyOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterCombinedOP.transpose()/(long double)NumSites << "\n";
}

void TriangularLattice::PrintThermalObservables(std::ostream &out){
  long double ns = NumSites;
  long double ns2 = pow(NumSites,2);
  long double ns3 = pow(NumSites,3);
  long double ns4 = pow(NumSites,4);

  out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  out << "Energy cumulants (C: E, E2, E3, E4) \n";
  out << std::setprecision(14) << EBar/ns << " " << E2Bar/ns2 << " " << E3Bar/ns3 << " " << E4Bar/ns4 << "\n";
  out << "Order parameter cumulants (R: FM, Perp, Par, Combined), (C: |m|, |m|2, |m|4)\n";
  out << std::setprecision(14) << FMNorm/ns << " " << FMNorm2/ns2 << " " << FMNorm4/ns4 << "\n";
  out << std::setprecision(14) << PerpNorm/ns << " " << PerpNorm2/ns2 << " " <<  PerpNorm4/ns4 << "\n";
  out << std::setprecision(14) << ParNorm/ns << " " << ParNorm2/ns2 << " " << ParNorm4/ns4 << "\n";
  out << std::setprecision(14) << CombinedNorm/ns << " " << CombinedNorm2/ns2 << " " << CombinedNorm4/ns4 << "\n";
}


void TriangularLattice::SampleConfiguration(double &temp, const uint& max_sweeps,
                                            const uint& sampling_time){
      // cout << temp << " " << temp << endl;
      long double m_e = 0;
      long double m_e2 = 0;
      long double m_e3 = 0;
      long double m_e4 = 0;

      long double m_fm_norm =0;
      long double m_perp_norm =0;
      long double m_par_norm =0;
      long double m_combined_norm =0;

      long double m_fm_norm2 =0;
      long double m_perp_norm2 =0;
      long double m_par_norm2 =0;
      long double m_combined_norm2 =0;

      long double m_fm_norm4 =0;
      long double m_perp_norm4 =0;
      long double m_par_norm4 =0;
      long double m_combined_norm4 =0;

      long double energy, fm_norm, perp_norm, par_norm, combined_norm;
      uint sweep = 0;
      uint samples = 0;
      uint accept;
      // cout << NumSites << " " << 0<<" " << 0<<" " << 0<<" " << 0<< endl;
      // uint overrelaxMCratio = 2;
      while (sweep < max_sweeps){
        DoTheSweeps(temp, accept);

        if (sweep%sampling_time == 0){
          energy = ClusterEnergy;
          // cout << samples << " " << ClusterEnergy/NumSites << endl;

          SelectStripyOP();

          m_e  += energy;
          m_e2 += pow(energy,2);
          m_e3 += pow(energy,3);
          m_e4 += pow(energy,4);

          fm_norm   = ClusterFMOP.norm();
          perp_norm = ClusterStripyOP.norm();
          par_norm  = abs(ClusterFMOP(1));
          combined_norm = ClusterCombinedOP.norm();

          m_fm_norm    += fm_norm;
          m_perp_norm  += perp_norm;
          m_par_norm   += par_norm;
          m_combined_norm += combined_norm;

          m_fm_norm2   += pow(fm_norm,2);
          m_perp_norm2 += pow(perp_norm,2);
          m_par_norm2  += pow(par_norm,2);
          m_combined_norm2  += pow(combined_norm,2);

          m_fm_norm4   += pow(fm_norm,4);
          m_perp_norm4 += pow(perp_norm,4);
          m_par_norm4  += pow(par_norm,4);
          m_combined_norm4  += pow(combined_norm,4);
          ++samples;

        }
        ++sweep;
      }

      EBar  = m_e/(long double)samples;
      // cout << "....." << endl;
      // cout << EBar << endl;
      // cout << "....." << endl;

      E2Bar = m_e2/(long double)samples;
      E3Bar = m_e3/(long double)samples;
      E4Bar = m_e4/(long double)samples;
      // cout << std::setprecision(14) << 0  << " " << ebar << " " << e2bar << endl;

      FMNorm       = m_fm_norm/(long double)samples;
      PerpNorm     = m_perp_norm/(long double)samples;
      ParNorm      = m_par_norm/(long double)samples;
      CombinedNorm = m_combined_norm/(long double)samples;

      FMNorm2       = m_fm_norm2/(long double)samples;
      PerpNorm2     = m_perp_norm2/(long double)samples;
      ParNorm2      = m_par_norm2/(long double)samples;
      CombinedNorm2 = m_combined_norm2/(long double)samples;

      FMNorm4       = m_fm_norm4/(long double)samples;
      PerpNorm4     = m_perp_norm4/(long double)samples;
      ParNorm4      = m_par_norm4/(long double)samples;
      CombinedNorm4 = m_combined_norm4/(long double)samples;
}

void TriangularLattice::CreateStripySignMatrices()
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

  StripySignsX = signs_X;
  StripySignsY = signs_Y;
  StripySignsZ = signs_Z;
}

void TriangularLattice::SelectStripyOP()
{
  Eigen::MatrixXd::Index indices[1];
  long double stripymax = ClusterStripyOPMatrix.rowwise().norm().maxCoeff(&indices[0]);
  ClusterStripyOP = ClusterStripyOPMatrix.row(indices[0]);

  Vector3LD combinedop;
  combinedop << ClusterStripyOP(0), ClusterFMOP(1),ClusterStripyOP(1);
  ClusterCombinedOP = combinedop;
}
