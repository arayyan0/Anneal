///  @file     triangular.hpp
///  @author   Ahmed Rayyan
///  @date     June 7, 2021
///  @brief    defining the triangular Lattice class
#include "triangular.hpp"

TriangularLattice::TriangularLattice(const uint& l1, const uint& l2,
                                     const double& jtau, const double& lambda):
L1(l1), L2(l2), NumSites(l1*l2), JTau(jtau), Lambda(lambda)
{
  //changing size of Cluster to L1*L2
  Cluster.resize(L1);
  for(uint i = 0; i < L1; ++i){
      Cluster[i].resize(L2);
  }

  // defining the bond Hamiltonians
  FixMPHamiltonians(jtau, lambda);

  std::uniform_int_distribution<uint> l1d(0, L1-1);
  std::uniform_int_distribution<uint> l2d(0, L2-1);

  L1Dist = l1d;
  L2Dist = l2d;
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

Eigen::Matrix3d TriangularLattice::ReturnMPHamiltonian(const double& angle,
                                                       const double& lambda)
//Eqn 17 of arXiv:2105.09334v1
{
  double c = cos(angle);
  double s = sin(angle);
  Eigen::Matrix3d matrix;
  matrix  <<  1-c-lambda,      0,          -s,
                       0, lambda,           0,
                      -s,      0,  1+c-lambda;
  return matrix;
}

void TriangularLattice::FixMPHamiltonians(const double& jtau, const double& lambda)
{
  double pz = 0;
  double px = 2*pi/3.0;
  double py = 4*pi/3.0;

  Hz = ReturnMPHamiltonian(pz, lambda)*jtau/2.0;
  Hx = ReturnMPHamiltonian(px, lambda)*jtau/2.0;
  Hy = ReturnMPHamiltonian(py, lambda)*jtau/2.0;
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

void TriangularLattice::CalculateLocalEnergy(const Site& site, double& energy)
{
  double e = 0;
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto nn_info : site.NearestNeighbours){
    auto [nn_x, nn_y, bond_type] = nn_info;
    spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
    // Bond bond(spin_i, spin_j, bond_type, jtau);
    // e += bond.BondEnergy;
    if (bond_type == 0){
      e += spin_i.VectorXYZ.transpose().dot(Hx*spin_j.VectorXYZ);
    } else if (bond_type == 1){
      e += spin_i.VectorXYZ.transpose().dot(Hy*spin_j.VectorXYZ);
    } else if (bond_type == 2){
      e += spin_i.VectorXYZ.transpose().dot(Hz*spin_j.VectorXYZ);
    }
  }
  energy = e;
}


void TriangularLattice::MolecularField(const Site& site, Eigen::Vector3d& molec)
{
  Eigen::Vector3d v = Eigen::Vector3d::Zero();
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto i : site.NearestNeighbours){
    auto [nn_x, nn_y, bond_type] = i;
    spin_j = Cluster[nn_x][nn_y].OnsiteSpin;
    // Bond bond(spin_i, spin_j, bond_type, p);
    // v += bond.MolecFieldContribution;
    if (bond_type == 0){
      v += -Hx*spin_j.VectorXYZ;
    } else if (bond_type == 1){
      v += -Hy*spin_j.VectorXYZ;
    } else if (bond_type == 2){
      v += -Hz*spin_j.VectorXYZ;
    }
  }
  molec = v;
}

void TriangularLattice::CalculateClusterEnergy()
{
  double e=0;
  double local_energy;
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      CalculateLocalEnergy(Cluster[x][y], local_energy);
      e += local_energy;
    }
  }
  ClusterEnergy = e/2.0;
}

void TriangularLattice::MetropolisSweep(const double& temperature)
{
  uint uc_x, uc_y;
  Spin old_spin_at_chosen_site;
  double old_local_energy, new_local_energy, energy_diff, r, pd;
  Site *chosen_site_ptr;

  uint flip = 0;
  double counter = 0;
  while (flip < NumSites){
    uc_x = L1Dist(MyRandom::RNG);
    uc_y = L2Dist(MyRandom::RNG);

    chosen_site_ptr = &Cluster[uc_x][uc_y];

    CalculateLocalEnergy(*chosen_site_ptr, old_local_energy);
    old_spin_at_chosen_site = chosen_site_ptr->OnsiteSpin;

    SpherePointPicker(chosen_site_ptr->OnsiteSpin);

    CalculateLocalEnergy(*chosen_site_ptr, new_local_energy);
    energy_diff = new_local_energy - old_local_energy;
    if (energy_diff < 0){++counter;}
    else{
      r = MyRandom::unit_interval(MyRandom::RNG);
      pd = exp(-energy_diff/temperature);
      if (r < pd){++counter;}
      else chosen_site_ptr->OnsiteSpin = old_spin_at_chosen_site;
    }
    ++flip;
  }
  // cout << temperature << " " << counter/(double)NumSites << endl;
  // cout << "..." << endl;
}

void TriangularLattice::SimulatedAnnealing(const uint& max_sweeps, double& initial_T,
                                          double& final_T)
{
  double scale = 0.95;
  double temp_T = initial_T;
  while(temp_T >= final_T){
    CalculateClusterEnergy();
    //std::cerr << temp_T << " " << std::setprecision(15) << ClusterEnergy/NumSites << endl;
    uint sweep = 0;
    while (sweep < max_sweeps){
      MetropolisSweep(temp_T);
      ++sweep;
    }
    temp_T = scale*temp_T;
  }
  final_T = temp_T/scale;
  FinalT = temp_T/scale;
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

        chosen_site_ptr = &Cluster[uc_x][uc_y];
        old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

        MolecularField(*chosen_site_ptr, molec_field);
        Spin new_spin(molec_field);
        chosen_site_ptr->OnsiteSpin = new_spin;

        norm = (new_spin.VectorXYZ - old_spin_vec).norm();
        if (norm > x){x = norm;}
        else{}
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
  out << "--------------------------------Results--------------------------------\n";
  out << "Energy per site\n";
  out << std::setprecision(14) << ClusterEnergy/NumSites << "\n";
  out << "Spin configuration\n";
  for (uint y=0; y<L2; ++y){
    for (uint x=0; x<L1; ++x){
      out << std::setprecision(14) << x << " " << y << " " << " " << Cluster[x][y].OnsiteSpin.VectorXYZ.transpose() << "\n";
    }
  }
}
