#include "fix.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "update.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace LAMMPS_NS;

class FixUpsideBridge : public Fix {
 public:
  FixUpsideBridge(LAMMPS *lmp, int narg, char **arg);
  ~FixUpsideBridge() override = default;
  int setmask() override;
  void post_integrate() override;

 private:
  int protein_groupbit, lipid_groupbit;
  std::string upside_cmd;
  void export_lipid_curvature();
  void import_protein_coordinates();
  void run_upside();
};

FixUpsideBridge::FixUpsideBridge(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg) {
  if (narg < 6) error->all(FLERR, "Illegal fix upside/bridge command");

  // Group names
  protein_groupbit = group->find(arg[3]);
  lipid_groupbit = group->find(arg[4]);
  if (protein_groupbit == -1 || lipid_groupbit == -1)
    error->all(FLERR, "Group ID not found");

  // Join command string (everything after arg[5])
  upside_cmd = "upside ";
  for (int i = 5; i < narg; ++i) {
    upside_cmd += std::string(arg[i]) + " ";
  }
  upside_cmd += "input.up";  // always append the input file
}

int FixUpsideBridge::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixUpsideBridge::post_integrate() {
  export_lipid_curvature();
  run_upside();
  import_protein_coordinates();
}

void FixUpsideBridge::export_lipid_curvature() {
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  std::vector<double> z_coords;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & lipid_groupbit) {
      z_coords.push_back(x[i][2]);
    }
  }

  double avg_z = 0.0;
  for (double z : z_coords) avg_z += z;
  avg_z /= z_coords.size();

  double curvature = 0.0;
  for (double z : z_coords) curvature += (z - avg_z) * (z - avg_z);
  curvature = sqrt(curvature / z_coords.size());

  std::ofstream fout("curvature_to_upside.dat");
  fout << curvature << std::endl;
  fout.close();
}

void FixUpsideBridge::run_upside() {
  int ret = system(upside_cmd.c_str());
  if (ret != 0) {
    char msg[512];
    snprintf(msg, sizeof(msg), "Upside execution failed with code %d", ret);
    error->warning(FLERR, msg);
  }
}

void FixUpsideBridge::import_protein_coordinates() {
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  std::ifstream fin("protein_coords_from_upside.dat");
  if (!fin.is_open()) return;

  std::string line;
  for (int i = 0; i < nlocal && std::getline(fin, line); i++) {
    if (!(mask[i] & protein_groupbit)) continue;

    std::istringstream ss(line);
    double newx, newy, newz;
    ss >> newx >> newy >> newz;
    x[i][0] = newx;
    x[i][1] = newy;
    x[i][2] = newz;
  }

  fin.close();
}