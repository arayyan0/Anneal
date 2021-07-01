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

  Hdefect2 = 0.5*Hdefect1;
  AddDefectHamiltonia();

  // cout << Hdefect1 << endl;
  // cout << Hdefect2 << endl;

  InitializeRandomSpins();

  std::uniform_int_distribution<uint> l1d(0, L1-1);
  std::uniform_int_distribution<uint> l2d(0, L2-1);

  L1Dist = l1d;
  L2Dist = l2d;

  CreateStripySignMatrices();
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
          std::get<2>(Cluster[x][y].NearestNeighbours[i]) = old_hamiltonian+Hdefect1;
        } else {
          poisoned_neighbour = CheckIfPoisoned(nn_x,nn_y);
          if (poisoned_neighbour == true){
            std::get<2>(Cluster[x][y].NearestNeighbours[i]) = old_hamiltonian+Hdefect1;

            uint ibefore = (i-1)%6;
            uint iafter = (i+1)%6;
            // oldbefore = Cluster[x][y].NearestNeighbours[ibefore];
            // oldafter = Cluster[x][y].NearestNeighbours[iafter ];
            std::get<2>(Cluster[x][y].NearestNeighbours[ibefore]) += Hdefect2;
            std::get<2>(Cluster[x][y].NearestNeighbours[iafter]) += Hdefect2;
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

void TriangularLattice::CalculateLocalEnergy(const Site& site, long double& energy)
{
  long double e = 0;
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto nn_info : site.NearestNeighbours){
    auto [nn_x, nn_y, ham] = nn_info;
    spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
    e += spin_j.VectorXYZ.transpose()*(ham*spin_i.VectorXYZ);
        - JTau*HField*HDirection.dot(spin_i.VectorXYZ + spin_j.VectorXYZ)/6.0;
  }
  energy = e;
}

void TriangularLattice::MolecularField(const Site& site, Vector3LDTrans& molec)
{
  Vector3LDTrans v = Vector3LD::Zero();
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto j : site.NearestNeighbours){
    auto [nn_x, nn_y, ham] = j;
    spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
    v += -spin_j.VectorXYZ.transpose()*ham + JTau*HField*HDirection.transpose()/6.0;
  }
  molec = v;
}

void TriangularLattice::CalculateClusterOP(){
  Vector3LD spin, FMOP, CombinedOP;
  Eigen::Matrix<long double, 3, 2> StripyOP;

  FMOP << 0,0,0;
  StripyOP << 0, 0,
              0, 0,
              0, 0;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      //calculate three stripy OPs
      spin = Cluster[n1][n2].OnsiteSpin.VectorXYZ;
      FMOP += spin;
      StripyOP(0,0) += StripySignsX(n1,n2)*spin(0); StripyOP(0,1) +=  StripySignsX(n1,n2)*spin(2);
      StripyOP(1,0) += StripySignsY(n1,n2)*spin(0); StripyOP(1,1) +=  StripySignsY(n1,n2)*spin(2);
      StripyOP(2,0) += StripySignsZ(n1,n2)*spin(0); StripyOP(2,1) +=  StripySignsZ(n1,n2)*spin(2);
    }
  }
  ClusterFMOP = FMOP;

  Eigen::MatrixXd::Index indices[1];
  long double stripymax = StripyOP.rowwise().norm().maxCoeff(&indices[0]);

  ClusterStripyOP = StripyOP.row(indices[0]);
  CombinedOP << ClusterStripyOP(0), ClusterFMOP(1),ClusterStripyOP(1);
  ClusterCombinedOP = CombinedOP;
}

void TriangularLattice::CalculateClusterEnergy()
{
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

  Vector3LD spin, FMOP, CombinedOP;
  Eigen::Matrix<long double, 3, 2> StripyOP;
  FMOP << 0,0,0;
  StripyOP << 0, 0,
              0, 0,
              0, 0;
  CombinedOP << 0,0,0;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      //calculate energy
      CalculateLocalEnergy(Cluster[n1][n2], local_energy);
      e += local_energy;
      //calculate three stripy OPs
      spin = Cluster[n1][n2].OnsiteSpin.VectorXYZ;

      FMOP += spin;

      StripyOP(0,0) += StripySignsX(n1,n2)*spin(0); StripyOP(0,1) += StripySignsX(n1,n2)*spin(2);
      StripyOP(1,0) += StripySignsY(n1,n2)*spin(0); StripyOP(1,1) += StripySignsY(n1,n2)*spin(2);
      StripyOP(2,0) += StripySignsZ(n1,n2)*spin(0); StripyOP(2,1) += StripySignsZ(n1,n2)*spin(2);
      }
    }

    ClusterEnergy = e/2.0;

    ClusterFMOP = FMOP;

    Eigen::MatrixXd::Index indices[1];
    long double stripymax = StripyOP.rowwise().norm().maxCoeff(&indices[0]);

    ClusterStripyOP = StripyOP.row(indices[0]);
    CombinedOP << ClusterStripyOP(0), ClusterFMOP(1),ClusterStripyOP(1);
    ClusterCombinedOP = CombinedOP;
}

void TriangularLattice::OverrelaxationFlip(){
  uint uc_x, uc_y;
  Site *chosen_site_ptr;
  Vector3LD old_spin_vec,new_spin_vec;

  Vector3LDTrans molec_field, normalized_field;

  uc_x = L1Dist(MyRandom::RNG);
  uc_y = L2Dist(MyRandom::RNG);

  chosen_site_ptr = &Cluster[uc_x][uc_y];
  old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

  cout << uc_x << " " << uc_y << endl;
  cout << "old spin " << old_spin_vec.transpose() << " with norm " << old_spin_vec.norm() << endl;


  MolecularField(*chosen_site_ptr, molec_field);
  normalized_field = molec_field.normalized();


  cout << "norm molec field " << normalized_field << " with norm " << normalized_field.norm() << endl;

  long double coeff = normalized_field*old_spin_vec;
  new_spin_vec = -old_spin_vec + 2.0*coeff*normalized_field.transpose();
  Spin new_spin(new_spin_vec);
  chosen_site_ptr->OnsiteSpin = new_spin;

  cout << "new spin " << (chosen_site_ptr->OnsiteSpin).VectorXYZ.transpose() << " with norm " << (chosen_site_ptr->OnsiteSpin).VectorXYZ.norm() << endl;

}

void TriangularLattice::OverrelaxationSweep(){
  uint flip = 0;
  while (flip < NumSites){
    CalculateClusterEnergy();
    long double e_i = ClusterEnergy;

    OverrelaxationFlip();

    CalculateClusterEnergy();
    long double e_f = ClusterEnergy;

    cout << flip << " " << e_i - e_f << endl;
    cout <<  "......................" << endl;
    ++flip;
  }
}

void TriangularLattice::MetropolisFlip(
  uint& uc_x, uint& uc_y,
  Spin&  old_spin_at_chosen_site,
  long double &old_local_energy, long double &new_local_energy, long double &energy_diff,
  double &r, double &pd, const double& temperature
)
{
  Site *chosen_site_ptr;
  uc_x = L1Dist(MyRandom::RNG);
  uc_y = L2Dist(MyRandom::RNG);

  chosen_site_ptr = &Cluster[uc_x][uc_y];

  CalculateLocalEnergy(*chosen_site_ptr, old_local_energy);
  old_spin_at_chosen_site = chosen_site_ptr->OnsiteSpin;

  SpherePointPicker(chosen_site_ptr->OnsiteSpin);

  CalculateLocalEnergy(*chosen_site_ptr, new_local_energy);
  energy_diff = new_local_energy - old_local_energy;
  if (energy_diff < 0){}
  else{
    r = MyRandom::unit_interval(MyRandom::RNG);
    pd = exp(-energy_diff/temperature);
    if (r < pd){}
    else chosen_site_ptr->OnsiteSpin = old_spin_at_chosen_site;
  }
}

void TriangularLattice::MetropolisSweep(const double& temperature)
{
  uint uc_x, uc_y;
  Spin old_spin_at_chosen_site;
  long double old_local_energy, new_local_energy, energy_diff;
  double r, pd;
  Site *chosen_site_ptr;

  uint flip = 0;

  while (flip < NumSites){
    MetropolisFlip(
      uc_x, uc_y,
      old_spin_at_chosen_site,
      old_local_energy, new_local_energy, energy_diff, r,
      pd,
      temperature
    );
    ++flip;
  }
}

void TriangularLattice::SimulatedAnnealing(const uint& max_sweeps,
                                          double& initial_T, double& final_T, double& rate)
{
  double scale = rate;
  double temp_T = initial_T;
  while(temp_T >= final_T){
    // std::cout << "temp " << temp_T << endl;
    uint sweep = 0;
    while (sweep < max_sweeps){
      // std::cout << "sweep " << sweep << endl;

      // for (uint i=0;i<10;i++){
      //   OverrelaxationSweep();
      // }

      MetropolisSweep(temp_T);
      ++sweep;
    }
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
  Vector3LDTrans molec_field;

  // long double x, norm;
  long double eps = pow(10,-20);

  while (sweep < max_sweeps){
    // x = 0;
    align = 0;
    while (align < NumSites){
        uc_x = L1Dist(MyRandom::RNG);
        uc_y = L2Dist(MyRandom::RNG);

        chosen_site_ptr = &Cluster[uc_x][uc_y];
        old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

        MolecularField(*chosen_site_ptr, molec_field);
        Vector3LD mf_t = molec_field.transpose();
        Spin new_spin(mf_t);
        chosen_site_ptr->OnsiteSpin = new_spin;

        // norm = (new_spin.VectorXYZ - old_spin_vec).norm();
        // if (norm > x){x = norm;}
        // else{}
        align++;
    }
    // if (x > eps){}
    // else{
      // ActualDetFlips = sweep;
      // break;
    // }
    sweep++;
  }
  ActualDetSweeps = sweep;
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
  out << "-------------------------------Final configuration observables--------------------------------\n";
  out << "Order parameters. R1: (FM_x, FM_y, FM_z), R2: (stripy_x, stripy_z),  R3: (stripy_x, FM_y, stripy_z)\n";
  out << std::setprecision(14) << ClusterFMOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterStripyOP.transpose()/(long double)NumSites << "\n";
  out << std::setprecision(14) << ClusterCombinedOP.transpose()/(long double)NumSites << "\n";
  //which=0 is ground state, which=1 is thermal
}

void TriangularLattice::PrintThermalObservables(std::ostream &out){
  long double ns = NumSites;
  long double ns2 = pow(NumSites,2);
  long double ns3 = pow(NumSites,3);
  long double ns4 = pow(NumSites,4);

  // long double ns =1;
  // long double ns2 = 1;
  // long double ns3 = 1;
  // long double ns4 = 1;

  out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  out << "Energy cumulants (C: E, E2, E3, E4) \n";
  out << std::setprecision(14) << EBar/ns << " " << E2Bar/ns2 << " " << E3Bar/ns3 << " " << E4Bar/ns4 << "\n";
  out << "Order parameter cumulants (R: FM, Perp, Par, Combined), (C: |m|, |m|2, |m|4)\n";
  out << std::setprecision(14) << FMNorm/ns << " " << FMNorm2/ns2 << " " << FMNorm4/ns4 << "\n";
  out << std::setprecision(14) << PerpNorm/ns << " " << PerpNorm2/ns2 << " " <<  PerpNorm4/ns4 << "\n";
  out << std::setprecision(14) << ParNorm/ns << " " << ParNorm2/ns2 << " " << ParNorm4/ns4 << "\n";
  out << std::setprecision(14) << CombinedNorm/ns << " " << CombinedNorm2/ns2 << " " << CombinedNorm4/ns4 << "\n";
  // out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  // out << "Energy cumulants (C: E, E2, E3, E4) \n";
  // out << std::setprecision(30) << EBar << " " << E2Bar << " " << E3Bar << " " << E4Bar << "\n";
  // out << "Order parameter cumulants (R: FM, Perp, Par, Combined), (C: |m|, |m|2, |m|4)\n";
  // out << std::setprecision(30) << FMNorm << " " << FMNorm2 << " " << FMNorm4 << "\n";
  // out << std::setprecision(30) << PerpNorm << " " << PerpNorm2 << " " <<  PerpNorm4 << "\n";
  // out << std::setprecision(30) << ParNorm << " " << ParNorm2 << " " << ParNorm4 << "\n";
  // out << std::setprecision(30) << CombinedNorm << " " << CombinedNorm2 << " " << CombinedNorm4 << "\n";
}

void TriangularLattice::ThermalizeConfiguration(double& temp, const uint& max_sweeps)
{
    // cout << temp << " " << temp << endl;

    uint sweep = 0;
    while (sweep < max_sweeps){
      MetropolisSweep(temp);

      // CalculateClusterEnergyandMagnetization();
      // cout << sweep << " " << ClusterEnergy/(double)NumSites  << endl;

      ++sweep;
    }
}

void TriangularLattice::SampleConfiguration(double &temp, const uint& max_sweeps, const uint& sampling_time){
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
      // cout << NumSites << " " << 0<<" " << 0<<" " << 0<<" " << 0<< endl;
      while (sweep < max_sweeps){
        MetropolisSweep(temp);
        if (sweep%sampling_time == 0){
          CalculateClusterEnergyandOP();
          energy = ClusterEnergy;

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
