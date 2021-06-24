///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#include "triangular.hpp"

TriangularLattice::TriangularLattice(const uint& l1, const uint& l2, const uint& num_defects,
                                     const double& jtau, const double& lambda,
                                     const double& ising_y,
                                     const double& defect,
                                     const double& h,
                                     Eigen::Vector3d& hdir):
L1(l1), L2(l2), NumSites(l1*l2), NumDefects(num_defects),
JTau(jtau), Lambda(lambda), IsingY(ising_y), Defect(defect),
HField(h), HDirection(hdir)
{

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

  InitializeRandomSpins();

  // defining the bond Hamiltonians
  FixMPHamiltonians();
  Hdefect << 0,      0, 0,
             0, Defect, 0,
             0,      0, 0;
  Hdefect=JTau*Hdefect;

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
                        {std::make_tuple(higher_x,        y, 2),
                         std::make_tuple(       x, higher_y, 1),
                         std::make_tuple( lower_x, higher_y, 0),
                         std::make_tuple( lower_x,        y, 2),
                         std::make_tuple(       x,  lower_y, 1),
                         std::make_tuple(higher_x,  lower_y, 0)},
                         some_spin
                      };
    }
  }
}

void TriangularLattice::CreateDefectPositions(){
  double pos;
  //only to be trusted when L1=L2=L

  for (uint d=0; d<NumDefects; ++d){
    pos = d*L1/(double)NumDefects + L1/2.0/(double)NumDefects;
    Defects[d] = {(uint)pos,(uint)pos};
  }

  // for (auto i:Defects){
  //   cout << i[0] << endl;
  // }
}

bool TriangularLattice::CheckIfPoisoned(uint lx, uint ly){
  bool truth = false;
  if (lx!=ly){
    return truth;
  }
  for (auto defectposition:Defects){
    if (defectposition[0]==lx) {
      truth = true;
      return truth;
    }
  }
  return truth;
}

Eigen::Matrix3d TriangularLattice::ReturnMPHamiltonian(const double& angle)
//Eqn 17 of arXiv:2105.09334v1
{
  double turn_jtau_off = 1;
  Eigen::Matrix3d matrix1, matrix2, matrix3;
  //bond-independent: lambda term + Ising in y direction
  matrix1  <<  -Lambda,      0,       0,
                     0, Lambda,       0,
                     0,      0, -Lambda;

  matrix2  <<             0,      0, 0,
                          0, IsingY, 0,
                          0,      0, 0;

  //bond-dependent: tau term
  double c = cos(angle);
  double s = sin(angle);
  matrix3  <<  1-c, 0,  -s,
                 0, 0,   0,
                -s, 0, 1+c;

  return matrix1/2.0 + matrix2 + turn_jtau_off*matrix3/2.0;
}

void TriangularLattice::FixMPHamiltonians()
{
  double pz = 0;
  double px = 2*pi/3.0;
  double py = 4*pi/3.0;

  Hz = ReturnMPHamiltonian(pz)*JTau;
  Hx = ReturnMPHamiltonian(px)*JTau;
  Hy = ReturnMPHamiltonian(py)*JTau;
}

void TriangularLattice::InitializeFMSpins(const double& theta, const double& phi)
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

  //first checks if I selected poisoned site or not
  if ( PoisonedSite_Flag == false ){
    //selected site is not poisoned
    //check each neighbour to see if that one IS poisoned
    double value;
    bool poisoned_neighbour;
    for (auto nn_info : site.NearestNeighbours){
      auto [nn_x, nn_y, bond_type] = nn_info;
      //check if neighbour is poisoned
      poisoned_neighbour = CheckIfPoisoned(nn_x, nn_y);
      //if neighbour is poisoned,
      double value = poisoned_neighbour ? 1.0 : 0.0;
      spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
      if (bond_type == 0){
        Ham = Hx + value*Hdefect;
      } else if (bond_type == 1){
        Ham = Hy + value*Hdefect;
      } else if (bond_type == 2){
        Ham = Hz + value*Hdefect;
      }
      e += spin_i.VectorXYZ.transpose().dot(Ham*spin_j.VectorXYZ)
          - JTau*HField*HDirection.transpose().dot(spin_i.VectorXYZ + spin_j.VectorXYZ)/6.0;
    }
  } else{
    //selected site IS poisoned
    //add defect to all neighbouring sites
    for (auto nn_info : site.NearestNeighbours){
      auto [nn_x, nn_y, bond_type] = nn_info;
      spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
      if (bond_type == 0){
        Ham = Hx + Hdefect;
      } else if (bond_type == 1){
        Ham = Hy + Hdefect;
      } else if (bond_type == 2){
        Ham = Hz + Hdefect;
      }
      e += spin_i.VectorXYZ.transpose().dot(Ham*spin_j.VectorXYZ)
          - JTau*HField*HDirection.transpose().dot(spin_i.VectorXYZ + spin_j.VectorXYZ)/6.0;
    }
  }
  energy = e;
}


void TriangularLattice::MolecularField(const Site& site, Eigen::Vector3d& molec)
{
  Eigen::Vector3d v = Eigen::Vector3d::Zero();
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  if (PoisonedSite_Flag == false) {
    double value;
    bool poisoned_neighbour;
    for (auto i : site.NearestNeighbours){
      auto [nn_x, nn_y, bond_type] = i;
      //check if neighbour is poisoned
      poisoned_neighbour = CheckIfPoisoned(nn_x, nn_y);
      //if neighbour is poisoned, return one. otherwise, return 0.
      double value = poisoned_neighbour ? 1.0 : 0.0;
      spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
      if (bond_type == 0){
        Ham = Hx+value*Hdefect;
      } else if (bond_type == 1){
        Ham = Hy+value*Hdefect;
      } else if (bond_type == 2){
        Ham = Hz+value*Hdefect;
      }
      v += -Ham*spin_j.VectorXYZ + JTau*HField*HDirection/6.0;
    }
  }else{
    for (auto i : site.NearestNeighbours){
      auto [nn_x, nn_y, bond_type] = i;
      spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
      if (bond_type == 0){
        Ham = Hx+Hdefect;
      } else if (bond_type == 1){
        Ham = Hy+Hdefect;
      } else if (bond_type == 2){
        Ham = Hz+Hdefect;
      }
      v += -Ham*spin_j.VectorXYZ+ JTau*HField*HDirection/6.0;
    }
  }
  molec = v;
}

void TriangularLattice::CalculateClusterOP(){
  Eigen::Vector3d spin, FMOP;
  Eigen::Matrix<double, 3, 2> StripyOP;
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
  Eigen::MatrixXd::Index index[1];
  double max = StripyOP.rowwise().norm().maxCoeff(&index[0]);

  ClusterFMOP = FMOP;
  ClusterStripyOP = StripyOP.row(index[0]);
}

void TriangularLattice::CalculateClusterEnergy()
{
  long double e=0;
  long double local_energy;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      //calculate energy
      PoisonedSite_Flag = CheckIfPoisoned(n1, n2);
      CalculateLocalEnergy(Cluster[n1][n2], local_energy);
      e += local_energy;
      PoisonedSite_Flag = false;
    }
  }
  ClusterEnergy = e/2.0;
}

void TriangularLattice::CalculateClusterEnergyandOP()
{
  long double e=0;
  long double local_energy;

  Eigen::Vector3d spin, FMOP;
  Eigen::Matrix<double, 3, 2> StripyOP;
  FMOP << 0,0,0;
  StripyOP << 0, 0,
              0, 0,
              0, 0;

  for (uint n1=0; n1<L1; ++n1){
    for (uint n2=0; n2<L2; ++n2){
      //calculate energy
      PoisonedSite_Flag = CheckIfPoisoned(n1, n2);
      CalculateLocalEnergy(Cluster[n1][n2], local_energy);
      e += local_energy;
      PoisonedSite_Flag = false;
      //calculate three stripy OPs
      spin = Cluster[n1][n2].OnsiteSpin.VectorXYZ;

      FMOP += spin;
      StripyOP(0,0) += StripySignsX(n1,n2)*spin(0); StripyOP(0,1) +=  StripySignsX(n1,n2)*spin(2);
      StripyOP(1,0) += StripySignsY(n1,n2)*spin(0); StripyOP(1,1) +=  StripySignsY(n1,n2)*spin(2);
      StripyOP(2,0) += StripySignsZ(n1,n2)*spin(0); StripyOP(2,1) +=  StripySignsZ(n1,n2)*spin(2);
      }
    }

  ClusterEnergy = e/2.0;
  Eigen::MatrixXd::Index index[1];
  long double max = StripyOP.rowwise().norm().maxCoeff(&index[0]);

  ClusterFMOP = FMOP;
  ClusterStripyOP = StripyOP.row(index[0]);
}


void TriangularLattice::MetropolisFlip(
  uint& uc_x, uint& uc_y,
  Spin&  old_spin_at_chosen_site,
  long double &old_local_energy, long double &new_local_energy, long double &energy_diff, double &r,
  double &pd,
  const double& temperature
)
{
  Site *chosen_site_ptr;
  uc_x = L1Dist(MyRandom::RNG);
  uc_y = L2Dist(MyRandom::RNG);

  PoisonedSite_Flag = CheckIfPoisoned(uc_x, uc_y);

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
  PoisonedSite_Flag = false;
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
  double x, norm;

  Eigen::Vector3d old_spin_vec, molec_field;
  double eps = pow(10,-20);

  while (sweep < max_sweeps){
    x = 0;
    align = 0;
    while (align < NumSites){
        uc_x = L1Dist(MyRandom::RNG);
        uc_y = L2Dist(MyRandom::RNG);

        PoisonedSite_Flag = CheckIfPoisoned(uc_x, uc_y);

        chosen_site_ptr = &Cluster[uc_x][uc_y];
        old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

        MolecularField(*chosen_site_ptr, molec_field);
        Spin new_spin(molec_field);
        chosen_site_ptr->OnsiteSpin = new_spin;

        norm = (new_spin.VectorXYZ - old_spin_vec).norm();
        if (norm > x){x = norm;}
        else{}
        PoisonedSite_Flag = false;
        align++;
    }
    // if (x > eps){}
    // else{
      // ActualDetFlips = sweep;
      // break;
    // }
    sweep++;
  }
  ActualDetFlips = sweep;
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
  out << "Order parameters (R1: (FM_x, FM_y, FM_z), R2: (stripy_x, stripy_z)\n";
  out << std::setprecision(14) << ClusterFMOP.transpose()/(double)NumSites << "\n";
  out << std::setprecision(14) << ClusterStripyOP.transpose()/(double)NumSites << "\n";
  //which=0 is ground state, which=1 is thermal
}

void TriangularLattice::PrintThermalObservables(std::ostream &out){
  out << "-------------------------------Thermal-averaged observables-----------------------------\n";
  out << "Energy cumulants (C: E, E2, E3, E4) \n";
  out << std::setprecision(14) << EBar << " " << E2Bar << " " << E3Bar << " " << E4Bar << "\n";
  out << "Order parameter cumulants (R: FM, Perp norm, Par norm), (C: |m|, |m|2, |m|4)\n";
  out << std::setprecision(14) << FMNorm << " " << FMNorm2 << " " << FMNorm4 << "\n";
  out << std::setprecision(14) << PerpNorm << " " << PerpNorm2 << " " <<  PerpNorm4 << "\n";
  out << std::setprecision(14) << ParNorm << " " << ParNorm2 << " " << ParNorm4 << "\n";
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
      long double e = 0;
      long double e2 = 0;
      long double e3 = 0;
      long double e4 = 0;

      long double m_fm_norm =0;
      long double m_perp_norm =0;
      long double m_par_norm =0;

      long double m_fm_norm2 =0;
      long double m_perp_norm2 =0;
      long double m_par_norm2 =0;

      long double m_fm_norm4 =0;
      long double m_perp_norm4 =0;
      long double m_par_norm4 =0;

      long double energydensity, fm_norm, perp_norm, par_norm;
      uint sweep = 0;
      uint samples = 0;
      while (sweep < max_sweeps){
        MetropolisSweep(temp);
        if (sweep%sampling_time == 0){
          CalculateClusterEnergyandOP();
          energydensity = ClusterEnergy/(double)NumSites;

          e  += energydensity;
          e2 += pow(energydensity,2);
          e3 += pow(energydensity,3);
          e4 += pow(energydensity,4);
          // cout << std::setprecision(14) << " " << samples << " " << pow(energydensity,3)  << " " << pow(energydensity,4)<< endl;

          fm_norm   = ClusterFMOP.norm()/(double)NumSites;
          perp_norm = ClusterStripyOP.norm()/(double)NumSites;
          par_norm  = abs(ClusterFMOP(1))/(double)NumSites;

          m_fm_norm    += fm_norm;
          m_perp_norm  += perp_norm;
          m_par_norm   += par_norm;
          m_fm_norm2   += pow(fm_norm,2);
          m_perp_norm2 += pow(perp_norm,2);
          m_par_norm2  += pow(par_norm,2);
          m_fm_norm4   += pow(fm_norm,4);
          m_perp_norm4 += pow(perp_norm,4);
          m_par_norm4  += pow(par_norm,4);
          ++samples;
        }
        ++sweep;
      }

      EBar  = e/((double)samples);
      E2Bar = e2/((double)samples);
      E3Bar = e3/((double)samples);
      E4Bar = e4/((double)samples);
      // cout << std::setprecision(14) << 0  << " " << ebar << " " << e2bar << endl;

      FMNorm = m_fm_norm/((double)samples);
      PerpNorm = m_perp_norm/((double)samples);
      ParNorm = m_par_norm/((double)samples);

      FMNorm2 = m_fm_norm2/((double)samples);
      PerpNorm2 = m_perp_norm2/((double)samples);
      ParNorm2 = m_par_norm2/((double)samples);

      FMNorm4 = m_fm_norm4/((double)samples);
      PerpNorm4 = m_perp_norm4/((double)samples);
      ParNorm4 = m_par_norm4/((double)samples);
}

void TriangularLattice::CreateStripySignMatrices()
{
  Eigen::ArrayXXd signs_Y(L1,L2);
  Eigen::ArrayXXd signs_Z(L1,L2);
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      x%2 == 0 ? signs_Y(x,y) = 1 : signs_Y(x,y) = -1;
      y%2 == 0 ? signs_Z(x,y) = 1 : signs_Z(x,y) = -1;
    }
  }
  Eigen::ArrayXXd signs_X = signs_Y*signs_Z;

  StripySignsX = signs_X;
  StripySignsY = signs_Y;
  StripySignsZ = signs_Z;
}
