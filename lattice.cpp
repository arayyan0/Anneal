///  @file     lattice.cpp
///  @author   Ahmed Rayyan
///  @date     December 2, 2019
///  @brief    constructing the lattice + simulated annealing algorithm
#include "lattice.hpp"

Lattice::Lattice(const uint& type, const uint& l1, const uint& l2, const uint& number_of_sublattices):
Type(type), L1(l1), L2(l2), NumSublattices(number_of_sublattices), NumUnitCells(l1*l2),
NumSites(l1*l2*number_of_sublattices)
{
  //changing size of Cluster to L1*L2*NumSublattices
  Cluster.resize(L1);
  for(uint i = 0; i < L1; ++i){
      Cluster[i].resize(L2);
      for (uint j = 0; j < L2; ++j) Cluster[i][j].resize(NumSublattices);
  }
  //selecting which cluster to create
  if (Type == 0){
    if (NumSublattices == 2) CreateRhombicCluster1();
    else if (NumSublattices == 4) CreateRectangularCluster2();
  }
  else if (Type == 1){
    if (NumSublattices == 6) CreateKekuleCluster();
  }

  std::uniform_int_distribution<uint> l1d(0, l1-1);
  std::uniform_int_distribution<uint> l2d(0, l2-1);
  std::uniform_int_distribution<uint> sd(0, number_of_sublattices-1);

  L1Dist = l1d;
  L2Dist = l2d;
  SubDist = sd;
}

void Lattice::CreateRhombicCluster1()
// T1 = a1, T2 = a2
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
      //see pg 2 of "Monte Carlo Cluster rules"
      Cluster[x][y][0] = {
                          {std::make_tuple(lower_x,       y, 1, 0),
                           std::make_tuple(      x, lower_y, 1, 1),
                           std::make_tuple(      x,       y, 1, 2)},
                          some_spin
                         };
      Cluster[x][y][1] = {
                          {std::make_tuple(higher_x,        y, 0, 0),
                           std::make_tuple(       x, higher_y, 0, 1),
                           std::make_tuple(       x,        y, 0, 2)},
                          some_spin
                         };
    }
  }
}

void Lattice::CreateRhombicCluster2()
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
      //see pg 2 of "Monte Carlo Cluster rules"
      Cluster[x][y][0] = {
                          {std::make_tuple(lower_x,       y, 1, 1),
                           std::make_tuple(      x, lower_y, 1, 2),
                           std::make_tuple(      x,       y, 1, 0)},
                          some_spin
                         };
      Cluster[x][y][1] = {
                          {std::make_tuple(higher_x,        y, 0, 1),
                           std::make_tuple(       x, higher_y, 0, 2),
                           std::make_tuple(       x,        y, 0, 0)},
                          some_spin
                         };
    }
  }
}

void Lattice::CreateRectangularCluster1()
// T1 = a1, T2 = 2a2-a1
{
  Spin some_spin;
  int higher_x, higher_y, lower_x, lower_y;
  for (int y=0; y<L2; ++y)
  {
    for (int x=0; x<L1; ++x)
    {
      higher_x = (x+1)%L1;
      higher_y = (y+1)%L2;
      //if index becomes -1, wrap it around to l-1
      lower_x = (x-1);
      lower_y = (y-1);
      if (lower_x >= 0){lower_x = lower_x%L1;}
      else if (lower_x < 0){lower_x = L1-1;}
      if (lower_y >= 0){lower_y = lower_y%L2;}
      else if (lower_y < 0){lower_y = L2-1;}
      // cout << lower_x << " " << lower_y << endl;      //see pg 2 of "Monte Carlo Cluster rules"
      Cluster[x][y][0] = {
                          {std::make_tuple(lower_x,       y, 1, 0),
                           std::make_tuple(      x, lower_y, 3, 1),
                           std::make_tuple(      x,       y, 1, 2)},
                          some_spin
                         };
      Cluster[x][y][1] = {
                          {std::make_tuple(higher_x,        y, 0, 0),
                           std::make_tuple(       x,        y, 2, 1),
                           std::make_tuple(       x,        y, 0, 2)},
                          some_spin
                         };
      Cluster[x][y][2] = {
                          {std::make_tuple(       x,       y, 3, 0),
                           std::make_tuple(       x,       y, 1, 1),
                           std::make_tuple(higher_x,       y, 3, 2)},
                          some_spin
                         };
      Cluster[x][y][3] = {
                          {std::make_tuple(      x,        y, 2, 0),
                           std::make_tuple(      x, higher_y, 0, 1),
                           std::make_tuple(lower_x,        y, 2, 2)},
                          some_spin
                         };
    }
  }
}

void Lattice::CreateRectangularCluster2()
// T1 = a1-a2, T2 = a1+a2
{
  Spin some_spin;
  int higher_x, higher_y, lower_x, lower_y;
  for (int y=0; y<L2; ++y)
  {
    for (int x=0; x<L1; ++x)
    {
      higher_x = (x+1)%L1;
      higher_y = (y+1)%L2;
      //if index becomes -1, wrap it around to l-1
      lower_x = (x-1);
      lower_y = (y-1);
      if (lower_x >= 0){lower_x = lower_x%L1;}
      else if (lower_x < 0){lower_x = L1-1;}
      if (lower_y >= 0){lower_y = lower_y%L2;}
      else if (lower_y < 0){lower_y = L2-1;}
      // cout << lower_x << " " << lower_y << endl;      //see pg 2 of "Monte Carlo Cluster rules"
      Cluster[x][y][0] = {
                          {std::make_tuple(lower_x,       y, 1, 1),
                           std::make_tuple(      x, lower_y, 3, 2),
                           std::make_tuple(      x,       y, 1, 0)},
                          some_spin
                         };
      Cluster[x][y][1] = {
                          {std::make_tuple(higher_x,        y, 0, 1),
                           std::make_tuple(       x,        y, 2, 2),
                           std::make_tuple(       x,        y, 0, 0)},
                          some_spin
                         };
      Cluster[x][y][2] = {
                          {std::make_tuple(       x,       y, 3, 1),
                           std::make_tuple(       x,       y, 1, 2),
                           std::make_tuple(higher_x,       y, 3, 0)},
                          some_spin
                         };
      Cluster[x][y][3] = {
                          {std::make_tuple(      x,        y, 2, 1),
                           std::make_tuple(      x, higher_y, 0, 2),
                           std::make_tuple(lower_x,        y, 2, 0)},
                          some_spin
                         };
    }
  }
}

void Lattice::CreateKekuleCluster()
{
  Spin some_spin;
  int higher_x, higher_y, lower_x, lower_y;
  for (int y=0; y<L2; ++y)
  {
    for (int x=0; x<L1; ++x)
    {
      higher_x = (x+1)%L1;
      higher_y = (y+1)%L2;
      //if index becomes -1, wrap it around to l-1
      lower_x = (x-1);
      lower_y = (y-1);
      if (lower_x >= 0){lower_x = lower_x%L1;}
      else if (lower_x < 0){lower_x = L1-1;}
      if (lower_y >= 0){lower_y = lower_y%L2;}
      else if (lower_y < 0){lower_y = L2-1;}
      //see pg 2 of "Monte Carlo Kekule Cluster rules"
      Cluster[x][y][0] = {
                          {std::make_tuple(x,        y, 1, 0),
                           std::make_tuple(higher_x, y, 3, 1),
                           std::make_tuple(x,        y, 5, 2)},
                          some_spin
                         };
      Cluster[x][y][1] = {
                          {std::make_tuple(x,       y, 0, 0),
                           std::make_tuple(x, lower_y, 4, 1),
                           std::make_tuple(x,       y, 2, 2)},
                          some_spin
                         };
      Cluster[x][y][2] = {
                          {std::make_tuple(      x,       y, 3, 0),
                           std::make_tuple(lower_x, lower_y, 5, 1),
                           std::make_tuple(      x,       y, 1, 2)},
                          some_spin
                         };
      Cluster[x][y][3] = {
                          {std::make_tuple(      x, y, 2, 0),
                           std::make_tuple(lower_x, y, 0, 1),
                           std::make_tuple(      x, y, 4, 2)},
                          some_spin
                         };
      Cluster[x][y][4] = {
                          {std::make_tuple(x,        y, 5, 0),
                           std::make_tuple(x, higher_y, 1, 1),
                           std::make_tuple(x,        y, 3, 2)},
                          some_spin
                          };
      Cluster[x][y][5] = {
                          {std::make_tuple(       x,        y, 4, 0),
                           std::make_tuple(higher_x, higher_y, 2, 1),
                           std::make_tuple(       x,        y, 0, 2)},
                          some_spin
                          };
    }
  }
}

void Lattice::InitializeRandomSpins()
{
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      for (uint sub=0; sub<NumSublattices; ++sub){
        SpherePointPicker(Cluster[x][y][sub].OnsiteSpin);
      }
    }
  }
}

void Lattice::InitializeFMSpins(const double& theta, const double& phi)
{
  Spin spin(theta, phi);
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      for (uint sub=0; sub<NumSublattices; ++sub){
        Cluster[x][y][sub].OnsiteSpin = spin;
      }
    }
  }
}

void Lattice::InitializeFromFile(const std::vector<std::string>& file_lines)
{
  uint nx, ny, s;
  double sx, sy, sz;
  for (uint i = 32; i < 32+NumSites; i++){
    std::stringstream line_ss(file_lines[i]);
    std::vector<double> line_double;
    double number;
    while (line_ss >> number)
      line_double.push_back(number);
    nx = std::round(line_double[0]);
    ny = std::round(line_double[1]);
    s = std::round(line_double[2]);
    sx = line_double[3];
    sy = line_double[4];
    sz = line_double[5];

    Eigen::Vector3d vector(sx, sy, sz);
    Spin spin(vector);
    Cluster[nx][ny][s].OnsiteSpin = spin;
  }
}

void Lattice::CalculateLocalEnergy(const Site& site, const Parameters& p, double& energy)
{
  double e = 0;
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto i : site.NearestNeighbours){
    auto [nn_x, nn_y, nn_sub, bond_type] = i;
    spin_j = Cluster[nn_x][nn_y][nn_sub].OnsiteSpin;
    Bond bond(spin_i, spin_j, bond_type, p);
    e += bond.BondEnergy;
  }
  energy = e;
}

void Lattice::MolecularField(const Site& site, const Parameters& p, Eigen::Vector3d& molec)
{
  Eigen::Vector3d v = Eigen::Vector3d::Zero();
  Spin spin_i = site.OnsiteSpin; Spin spin_j;
  for (auto i : site.NearestNeighbours){
    auto [nn_x, nn_y, nn_sub, bond_type] = i;
    spin_j = Cluster[nn_x][nn_y][nn_sub].OnsiteSpin;
    Bond bond(spin_i, spin_j, bond_type, p);
    v += bond.MolecFieldContribution;
  }
  molec = v;
}

void Lattice::CalculateClusterEnergy(const Parameters& p)
{
  double e=0;
  double local_energy;
  for (uint x=0; x<L1; ++x){
    for (uint y=0; y<L2; ++y){
      for (uint sub=0; sub<NumSublattices; ++sub){
        CalculateLocalEnergy(Cluster[x][y][sub], p, local_energy);
        e += local_energy;
      }
    }
  }
  ClusterEnergy = e/2.0;
}

void Lattice::MetropolisSweep(const Parameters& p, const double& temperature)
{
  uint uc_x, uc_y, sublattice;
  Spin old_spin_at_chosen_site;
  double old_local_energy, new_local_energy, energy_diff, r, pd;
  Site *chosen_site_ptr;

  uint flip = 0;
  double counter = 0;
  while (flip < NumSites){
    uc_x = L1Dist(MyRandom::RNG);
    uc_y = L2Dist(MyRandom::RNG);
    sublattice = SubDist(MyRandom::RNG);

    chosen_site_ptr = &Cluster[uc_x][uc_y][sublattice];

    CalculateLocalEnergy(*chosen_site_ptr, p, old_local_energy);
    old_spin_at_chosen_site = chosen_site_ptr->OnsiteSpin;

    SpherePointPicker(chosen_site_ptr->OnsiteSpin);

    CalculateLocalEnergy(*chosen_site_ptr, p, new_local_energy);
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

void Lattice::SimulatedAnnealing(const Parameters& p, const uint& max_sweeps,
                                 double& initial_T, double& final_T)
{
  double scale = 0.9;
  double temp_T = initial_T;
  while(temp_T >= final_T){
    CalculateClusterEnergy(p);
    cout << temp_T << " " << std::setprecision(14) << ClusterEnergy/NumSites << endl;
    uint sweep = 0;
    while (sweep < max_sweeps){
      MetropolisSweep(p, temp_T);
      ++sweep;
    }
    temp_T = scale*temp_T;
  }
  final_T = temp_T/scale;
  FinalT = temp_T/scale;
}

void Lattice::DeterministicSweeps(const Parameters& p, const uint& max_sweeps)
{
  uint sweep = 0;
  uint align;
  uint uc_x, uc_y, sublattice;
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
        sublattice = SubDist(MyRandom::RNG);

        chosen_site_ptr = &Cluster[uc_x][uc_y][sublattice];
        old_spin_vec = (chosen_site_ptr->OnsiteSpin).VectorXYZ;

        MolecularField(*chosen_site_ptr, p, molec_field);
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

void Lattice::PrintConfiguration(std::ostream &out)
{
  out << "--------------------------------Results--------------------------------\n";
  out << "Energy per site\n";
  out << std::setprecision(14) << ClusterEnergy/NumSites << "\n";
  out << "Spin configuration\n";
  for (uint y=0; y<L2; ++y){
    for (uint x=0; x<L1; ++x){
      for (uint sub=0; sub<NumSublattices; ++sub){
        out << std::setprecision(14) << x << " " << y << " " << sub << " " << Cluster[x][y][sub].OnsiteSpin.VectorXYZ.transpose() << "\n";
      }
    }
  }
}
