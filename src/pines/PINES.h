/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2021, Andrea Arsiccio

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_PINES_vec_PINES_h
#define __PLUMED_PINES_vec_PINES_h
#include <queue>

using namespace std;

namespace PLMD {
namespace PINES {
// Ideally core/Colvar.h should be moved to this directory and Colvar should stay in namespace PLMD::Sasa
// With this trick, PLMD::Colvar is visible as PLMD::Sasa::Colvar
using PLMD::Colvar;

class PINES      : public Colvar
{
private:
  bool pbc, serial, timer;
  ForwardDecl<Stopwatch> stopwatch_fwd;
  Stopwatch& stopwatch=*stopwatch_fwd;
  // Added NL_const_size to fix solute-solvent elements as constant size
  unsigned NL_const_size;
  int updatePINES, N_Blocks, steps_since_update, nstride;
  size_t Nprec;
  unsigned Natm,Nlist,NLsize;
  double Fvol,Vol0;
  double m_PINESdistance;
  unsigned solv_blocks;
  std::string ref_file;
  NeighborList *nlall;
  NeighborList *nlreduced;
  std::vector<SwitchingFunction> sfs;
  std::vector<std:: vector<double> > rPINES;
  std::vector<double> scaling,r00;
  std::vector<double> nl_skin;
  std::vector<double> fmass;
  std::vector<bool> dosort;
  std::vector<Vector> compos;
  std::vector<string> sw;
  std::vector<NeighborList *> nl;
  std::vector<NeighborList *> nl_small;
  std::vector<NeighborList *> nlcom;
  std::vector<Vector> m_deriv;
  std::vector<std:: vector<double> > PIV;

  // ds_array is the 1D array (dv(r)/dr) of the switching function --NH
  std::vector<double> ds_array;
  // ann_deriv is the 3D array (dv(r)/dxyz) passed to the plumed core --NH
  std::vector<std:: vector<Vector> > ann_deriv;

  // The PINES_Pair vectors record the atom IDs for the PINES elements that are passed to the VAE --NH
  std::vector<int> PINES_Pair0;
  std::vector<int> PINES_Pair1;

  std::vector<std::vector<AtomNumber>> Plist;
  std::vector<AtomNumber> listall;
  std::vector<unsigned> AtomToResID_Dict;
  std::vector<unsigned> NList_OW_blocks;
  std::vector<unsigned> NList_HW_blocks;
  std::vector<AtomNumber> listreduced;
  std::vector<AtomNumber> listnonwater;
  std::vector<string> atype;
  std::vector<double> nl_cut;
  std::vector<int> nl_st;

  // adding a flag (cart2PINES) for post-processing a trajectory in cartesian coordinates to a PINES representation
  // -- SD flag for writing a single file containing PINES values when using plumed driver [writePINEStraj,writestride];
  bool Svol,cross,direct,doneigh,test,CompDer,writePINEStraj,writestride;
  // -- SD variables to control output PINES and ANN PINES derivative file during simulation.
  int writePINESstride, writeannstride; 
  bool cart2PINES, com; // NOTE com not tested/used -- SD

  // -- SD variables in prepare() function.
  bool invalidateList,firsttime,enableLog;
  
  Tensor m_virial; // [NOT SUPPORTED/USED -- SD]

  std::vector<string> block_params;
  std::vector<std::vector<std::vector<AtomNumber>>> block_groups_atom_list;
  std::vector<int> block_lengths;
  std::vector<int> G1_limits;
  std::vector<int> Buffer_Pairs;
  std::vector<int> tot_num_pairs;
  std::vector<std::vector<std::pair<AtomNumber,AtomNumber>>> Exclude_Pairs;
  std::vector<std::vector<std::vector<AtomNumber>>> ID_list;
  std::vector<std::vector<std::vector<int>>> ResID_list;
  std::vector<std::vector<std::vector<string>>> Name_list;
  std::vector<std::vector<std::vector<bool>>> filters;
  std::vector<std::vector<std::pair<AtomNumber,AtomNumber>>>pairlist;
  
  //std::vector<std::priority_queue<std::pair<double,std::pair<AtomNumber,AtomNumber>>,std::vector<std::pair<double,std::pair<AtomNumber,AtomNumber>>>,CompareDist>>maxHeapVec;


  // dr_dxyz_array is the 3D array (dr/dxyz) used to build ann_deriv and ANN_sum_array --NH [NOT USED -- SD]
  // std::vector<std:: vector<Vector> > dr_dxyz_array;
  // ANN_sum_array is the 1D array (sum dv_d/dv_n) written to an output file for use by the ANN code --NH
  //std::vector<double> ANN_sum_array;
  // ANN PINES derivatives array written to output file for use by ANN code --SD
  // std::vector<std::vector<double>> ANN_PINES_deriv; [NOT USED -- SD]
public:
  static void registerKeywords( Keywords& keys );                                                                       
  explicit PINES(const ActionOptions&); 
  //~PINES();                                                                                                               
  // active methods:
  struct CompareDist {
    bool operator()(const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p1, const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p2) {
      return p1.first < p2.first; // Max heap
    }
  };                                                                                             
  virtual void calculate();
  void checkFieldsAllowed() {}                                                                                           
  // -- SD prepare to requestAtoms during simulation 
  void prepare() override;
  // -- SD ANN SUM DERIVATIVE                                                                                              
  std::vector<vector<double>> get_ann_sum_derivative();
  void Update_NL();
};

}
}

#endif
