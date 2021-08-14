// Dien N., Aug-11, 2021
#ifndef __DCFID_SIDIS_H__
#define __DCFID_SIDIS_H__

// I want to write a DC-fid class for all particle
// Testing first then clean up later
// Trying only XY first
// Now fid depends on the part's ID
//

#include "TF1.h"
#include "TVector3.h"
#include <unordered_map>
#include <utility>

class DCfid_SIDIS {
public:
  DCfid_SIDIS();
  virtual ~DCfid_SIDIS();

  TVector3 rotate(double x, double y, int sector) const;
  bool DC_fid_xy_sidis(int pid, double x, double y, int sector, int layer,
                       bool bending) const;

  bool DC_fid_th_ph_sidis(int pid, double x, double y, double z, int sector,
                          int layer, bool bending) const;

  std::pair<double, double> cal_th_ph(double x, double y, double z,
                                      int sector) const;

  TF1 *get_fmin_xy_in(int pid, int sector, int layer) const;
  TF1 *get_fmax_xy_in(int pid, int sector, int layer) const;
  TF1 *get_fmin_xy_out(int pid, int sector, int layer) const;
  TF1 *get_fmax_xy_out(int pid, int sector, int layer) const;

  TF1 *get_fmin_ph_in(int pid, int sector, int layer) const;
  TF1 *get_fmax_ph_in(int pid, int sector, int layer) const;

protected:
  TF1 *fmin_S_L_xy_in[6][6][3];
  TF1 *fmax_S_L_xy_in[6][6][3];
  TF1 *fmin_S_L_xy_out[6][6][3];
  TF1 *fmax_S_L_xy_out[6][6][3];

  TF1 *fmin_S_L_ph_in[6][6][3];
  TF1 *fmax_S_L_ph_in[6][6][3];

  static const double maxparams_xy_in[6][6][3][2];
  static const double minparams_xy_in[6][6][3][2];
  static const double maxparams_xy_out[6][6][3][2];
  static const double minparams_xy_out[6][6][3][2];

  static const double maxparams_th_ph_in[6][6][3][4];
  static const double minparams_th_ph_in[6][6][3][4];

  static const std::unordered_map<int, int> pid_map;

  // const double Pival = 3.14;
};

#endif
