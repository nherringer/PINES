/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 of Nicholas Herringer and Siva Dasetty.

The PINES module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The PINES module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

// -- SD header file for PINES
#include "PINES.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <queue>

using namespace std;

namespace PLMD
{
namespace PINES
{

PLUMED_REGISTER_ACTION(PINES,"PINES")

void PINES::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered","SWITCH","The switching functions parameter."
           "You should specify a Switching function for all PINES blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("numbered","BLOCK","Each block of the PIV");
  keys.add("compulsory","PRECISION","the precision for approximating reals with integers in sorting.");
  keys.add("compulsory","REF_FILE","PDB file name that contains the information about system connectivity and labels.");
  keys.add("compulsory","N_BLOCKS","Number of blocks in PIV");
  keys.add("compulsory","SIZE","Length of each PIV Block");
  keys.add("optional","ATOMID","AtomIDs");
  keys.add("optional","RESID","ResIDs");
  keys.add("optional","NAME","Atom Names");
  keys.add("optional","EXLCUDE_PAIRS","Excluded pairs");
  keys.add("optional","LIMIT_PER_G1","NN for each group 1 atom");
  keys.add("optional","BUFFER","Number of additional pairwise distances to include as a buffer for each PIV block");
  keys.add("optional","UPDATEPINES","Frequency (in steps) at which the PINES is updated.");
  keys.addFlag("DERIVATIVES",false,"Activate the calculation of the PINES for every class (needed for numerical derivatives).");
  keys.add("optional","NL_CUTOFF","Neighbor lists cutoff.");
  keys.add("optional","NL_STRIDE","Update neighbor lists every NL_STRIDE steps. When post-processing a trajectory a stride of 1 should always be used.");
  keys.add("optional","NL_SKIN","The maximum atom displacement tolerated for the neighbor lists update.");
  // -- SD Flag for writing PINES values in a single file when using plumed driver.
  keys.addFlag("WRITEPINESTRAJ",false,"Flag to enable or disable writing PINES_representation when using plumed driver.");
  // -- SD Variables to control frequency of writing PINES values and ANN PINES derivatives during simulation.
  keys.add("optional","WRITEPINESSTRIDE","STRIDE to write PINES_representation.");
  componentsAreNotOptional(keys);
  // Changing "COMPONENTS" to "default" and slightly modifying the name. Added components for ANN_SUM_DERIV
  keys.addOutputComponent("ELEMENT", "default", "Elements of the PINES block. The position in the N choose 2 interactions (i) and the neighbor in the neighbor list (j) is given as PINES-i-j.");
  //keys.addOutputComponent("ANNSUMDERIV", "default", "2D array of PINES element partial derivatives (used with ANN module).");
  keys.addFlag("ENABLE_LOG",false,"Flag to enable logging.");
  keys.reset_style("SWITCH","compulsory");
}

PINES::PINES(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  updatePINES(1),
  Nprec(1000),
  Natm(1),
  Nlist(1),
  NLsize(1),
  Fvol(1.),
  m_PINESdistance(0.),
  solv_blocks(1),
  r00(std:: vector<double>(Nlist)),
  nl_skin(std:: vector<double>(Nlist)),
  fmass(std:: vector<double>(Nlist)),
  sw(std:: vector<string>(Nlist)),
  nl(std:: vector<NeighborList *>(Nlist)),
  nl_small(std:: vector<NeighborList *>(Nlist)),
  m_deriv(std:: vector<Vector>(1)),
  ds_array(std:: vector<double>(1)),
  ann_deriv(std:: vector<std:: vector<Vector> >(1)),
  PIV(std:: vector<std:: vector<double> >(1)),
  PINES_Pair0(std:: vector<int>(1)),
  PINES_Pair1(std:: vector<int>(1)),
  Plist(std:: vector<std:: vector<AtomNumber> >(1)),
  listall(std:: vector<AtomNumber>(1)),
  listreduced(std:: vector<AtomNumber>(1)),
  nl_cut(std:: vector<double>(Nlist)),
  nl_st(std:: vector<int>(Nlist)),
  CompDer(false),
  // SD -- local variables corresponding to user defined flags.
  writePINEStraj(false),
  writestride(false),
  writePINESstride(-1),
  cart2PINES(true),
  // SD -- used in prepare function.
  invalidateList(true),
  firsttime(true),
  enableLog(false),
  block_params(std:: vector<string>(1)),
  block_groups_atom_list(std:: vector<std:: vector<std:: vector<AtomNumber> > >(1)),
  block_lengths(std:: vector<int>(1)),
  G1_limits(std:: vector<int>(1)),
  Buffer_Pairs(std:: vector<int>(1)),
  Exclude_Pairs(std:: vector<std:: vector<std:: pair<AtomNumber> > >(1)),
  ID_list(std:: vector<std:: vector<std:: vector<AtomNumber> > >(1)),
  ResID_list(std:: vector<std:: vector<std:: vector<int> > >(1)),
  Name_list(std:: vector<std:: vector<std:: vector<string> > >(1)),
  filters(std:: vector<std:: vector<std:: vector<bool> > >(1)),
  pairlist(std:: vector<std:: pair<AtomNumber, AtomNumber> > (1)),
  N_Blocks(1)
{
  if(keywords.exists("ENABLE_LOG")) {
    parseFlag("ENABLE_LOG",enableLog);
  }

  if (enableLog) {
  log << "Starting PINES Constructor\n";
  }

  // Precision on the real-to-integer transformation for the sorting
  parse("PRECISION",Nprec);
  if(Nprec<2) error("Precision must be => 2");

  // PBC
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if (enableLog) {
    if(pbc) {
      log << "Using Periodic Boundary Conditions\n";
    } else  {
      log << "Isolated System (NO PBC)\n";
    }
  }

  // SERIAL/PARALLEL
  parseFlag("SERIAL",serial);
  if (enableLog) {
    if(serial) {
      log << "Serial PINES construction\n";
    } else     {
      log << "Parallel PINES construction\n";
    }
  }

  // Derivatives
  parseFlag("DERIVATIVES",CompDer);
  if (enableLog) {
  if(CompDer) log << "Computing Derivatives\n";
  }

  // Timing
  parseFlag("TIMER",timer);
  if(timer) {
    if (enableLog) {
      log << "Timing analysis\n";
    }
    stopwatch.start();
    stopwatch.pause();
  }

  // Test
  parseFlag("TEST",test);

  // Constant Neighbor List Size
  if(keywords.exists("NL_CONSTANT_SIZE")) {
    parse("NL_CONSTANT_SIZE",NL_const_size);
  }

  // UPDATEPINES
  if(keywords.exists("UPDATEPINES")) {
    parse("UPDATEPINES",updatePINES);
  }

  // Volume Scaling
  parse("VOLUME",Vol0);
  if (Vol0>0) {
    Svol=true;
  }

  

  // PINES scaled option
  scaling.resize(Nlist);
  for(unsigned j=0; j<Nlist; j++) {
    scaling[j]=1.;
  }

  if(keywords.exists("SFACTOR")) {
    parseVector("SFACTOR",scaling);
  }

  // Added STRIDE to write PINES representation and ANN sum derivatives -- SD
  if(keywords.exists("WRITEPINESTRAJ")){
    parseFlag("WRITEPINESTRAJ",writePINEStraj);
  }
  if(keywords.exists("WRITEPINESSTRIDE")) { 
    parse("WRITEPINESSTRIDE",writePINESstride);
  }
  if (writePINESstride != -1) {
    writestride=true;
  }

  // Reference PDB file from which atom names, types, ids, and initial positions are determined
  parse("REF_FILE",ref_file);
  PDB mypdb;
  FILE* fp=fopen(ref_file.c_str(),"r");
  if (fp!=NULL) {
    if (enableLog) {
      log<<"Opening PDB file with reference frame: "<<ref_file.c_str()<<"\n";
    }
    mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    fclose (fp);
  } else {
    error("Error in reference PDB file");
  }

  // Create variable to get number of blocks
  parse("N_BLOCKS",N_Blocks);

  // These all need to be initialized to take the relevant parameters
  // from all blocks
  block_params.resize(N_Blocks);
  block_groups_atom_list.resize(N_Blocks);
  block_lengths.resize(N_Blocks);
  G1_limits.resize(N_Blocks);
  Exclude_Pairs.resize(N_Blocks);
  ID_list.resize(N_Blocks);
  ResID_list.resize(N_Blocks);
  Name_list.resize(N_Blocks);

  for (unsigned n=0; n<N_Blocks; n++) {
    G1_limits[n]=0;
    ID_list[n].resize(2);
    ResID_list[n].resize(2);
    Name_list[n].resize(2);
  }

  // Check that the correct number of Blocks are specified
  for (unsigned n=0; n<N_Blocks; n++) {
    if( !parseNumbered( "BLOCK", n+1, block_params[n] ) ) break;
  }

  FILE *parse_test = NULL;
                                                                              
  string parse_test_fileName = "ParseTest.dat";                                    
  parse_test = fopen(parse_test_fileName.c_str(), "w+");

  // Parse blocks for keywords
  for(int n=0; n<N_Blocks; n++){
    int block_length;
    int limit_per_G1=0;
    std::vector<AtomNumber> ex_pairs_n;
    std::vector<AtomNumber> g1_ids;
    std::vector<AtomNumber> g2_ids;
    std::vector<int> g1_resids;
    std::vector<int> g2_resids;
    std::vector<string> g1_names;
    std::vector<string> g2_names;
    // fprintf(parse_test, "Block Params (%d):\n", n);
    // fprintf(parse_test, "%s\n", block_params[n].c_str());
    std::vector<string> block_data=Tools::getWords(block_params[n]);
    // fprintf(parse_test, "Block Data (%d):\n", n);
    // for(int x=0; x<block_data.size(); x++){
    //   fprintf(parse_test, "%s\n", block_data[x].c_str());
    // }
    //vector<string> G1_data=Tools::getWords(block_data[0]);
    std::vector<string> G1_data;
    Tools::parseVector(block_data, "G1", G1_data);
    // fprintf(parse_test, "Block Data [0] (G1):\n", n);
    // for(int x=0; x<G1_data.size(); x++){
    //   fprintf(parse_test, "%s\n", G1_data[x].c_str());
    // }
    // fprintf(parse_test, "Block Data [1] (G2):\n", n);
    std::vector<string> G2_data;
    Tools::parseVector(block_data, "G2", G2_data);
    // for(int x=0; x<G2_data.size(); x++){
    //   fprintf(parse_test, "%s\n", G2_data[x].c_str());
    // }
    Tools::parseVector(G1_data,"ATOMID",g1_ids);
    // fprintf(parse_test, "g1_ids size: %d\n", g1_ids.size());
    for(int i=0; i<g1_ids.size(); i++){
      ID_list[n][0].push_back(g1_ids[i]);
      // fprintf(parse_test, "%d\n", g1_ids[i]);
    }

    Tools::parseVector(G2_data,"ATOMID",g2_ids);
    // fprintf(parse_test, "g2_ids size: %d\n", g2_ids.size());
    for(int i=0; i<g2_ids.size(); i++){
      ID_list[n][1].push_back(g2_ids[i]);
      // fprintf(parse_test, "%d\n", g2_ids[i]);
    }
    Tools::parseVector(G1_data,"RESID",g1_resids);
    // fprintf(parse_test, "g1_resids size: %d\n", g1_resids.size());
    for(int i=0; i<g1_resids.size(); i++){
      ResID_list[n][0].push_back(g1_resids[i]);
      // fprintf(parse_test, "%d\n", g1_resids[i]);
    }
    Tools::parseVector(G2_data,"RESID",g2_resids);
    // fprintf(parse_test, "g2_resids size: %d\n", g2_resids.size());
    for(int i=0; i<g2_resids.size(); i++){
      ResID_list[n][1].push_back(g2_resids[i]);
      // fprintf(parse_test, "%d\n", g2_resids[i]);
    }
    Tools::parseVector(G1_data,"NAME",g1_names);
    // fprintf(parse_test, "g1_names size: %d\n", g1_names.size());
    for(int i=0; i<g1_names.size(); i++){
      Name_list[n][0].push_back(g1_names[i]);
      // fprintf(parse_test, "%s\n", g1_names[i].c_str());
    }
    Tools::parseVector(G2_data,"NAME",g2_names);
    // fprintf(parse_test, "g2_names size: %d\n", g2_names.size());
    for(int i=0; i<g2_names.size(); i++){
      Name_list[n][1].push_back(g2_names[i]);
      // fprintf(parse_test, "%s\n", g2_names[i].c_str());
    }
    
    // fprintf(parse_test, "%d\n", g1_ids);
    // fprintf(parse_test, "%d\n", g2_ids);
    // fprintf(parse_test, "%d\n", g1_resids);
    // fprintf(parse_test, "%d\n", g2_resids);


    Tools::parse(block_data, "SIZE", block_length);
    block_lengths[n]=block_length;
    int buffer_pairs=block_length*2;
    Tools::parseVector(block_data,"EXCLUDE_PAIRS",ex_pairs_n);
    for(int i=0; i<ex_pairs_n.size()-1; i++){
      std::pair<AtomNumber,AtomNumber> excluded_pair = {ex_pairs_n[i],ex_pairs_n{i+1}};
      Exclude_Pairs[n].push_back(excluded_pair);
    }
    Tools::parse(block_data,"LIMIT_PER_G1",limit_per_G1);
    Tools::parse(block_data,"BUFFER",buffer_pairs);
    fprintf(parse_test, "Limit per G1: %d\n", limit_per_G1);
    fprintf(parse_test, "Buffer Pairs: %d\n", buffer_pairs);
    if(limit_per_G1>0){
      G1_limits[n]=limit_per_G1;
    }
    Buffer_Pairs[n]=buffer_pairs;
  }

  fprintf(parse_test, "N_Blocks = %d\n\n\n", N_Blocks);
  //fprintf(parse_test, "Block Params:\n");
  for(int n=0; n<N_Blocks; n++){
    //fprintf(parse_test, "%s\n", block_params[n].c_str());
    fprintf(parse_test, "Block %d\n", n);
    for(int g=0; g<2; g++){
      fprintf(parse_test, "ID List G%d:\n", (g+1));
      for(int f=0; f<ID_list[n][g].size(); f++){
        fprintf(parse_test, "%d\t", ID_list[n][g][f]);
      }
      fprintf(parse_test, "\n");
    }
    for(int g=0; g<2; g++){
      fprintf(parse_test, "RESID List G%d:\n", (g+1));
      for(int f=0; f<ResID_list[n][g].size(); f++){
        fprintf(parse_test, "%d\t", ResID_list[n][g][f]);
      }
      fprintf(parse_test, "\n");
    }
    for(int g=0; g<2; g++){
      fprintf(parse_test, "Name List G%d:\n", (g+1));
      for(int f=0; f<Name_list[n][g].size(); f++){
        fprintf(parse_test, "%s\t", Name_list[n][g][f].c_str());
      }
      fprintf(parse_test, "\n");
    }
    fprintf(parse_test, "Block Length: %d\n", block_lengths[n]);
    fprintf(parse_test, "Exclude Pairs: \n");
    for(int f=0; f<Exclude_Pairs[n].size(); f++){
      fprintf(parse_test, "%d\t", Exclude_Pairs[n][f]);
    }
    fprintf(parse_test, "\n");
  }
  fclose(parse_test);
  exit();

  // In previous section build filter lists to use in building atom lists
  // ID_list, Resid_list, and Name_list are all vectors of dimension N_blocks x 2 (G1/G2) x n components
  // The union of filters is then easy by checking for each block if more than one list has entries

  filters.resize(N_Blocks); //N blocks, two groups, three filters
  for(int n=0; n<N_Blocks; n++){
    filters[n].resize(2);
  }
  for(int n=0; n<N_Blocks; n++){
    for(int g=0; g<2; g++){
      filters[n][g].resize(3);
    }
  }
  for(int n=0; n<N_Blocks; n++){
    for(int g=0; g<2; g++){
      for(int f=0; f<3; f++){
        filters[n][g][f]=false;
      }
    }
  }
  for(int n=0; n<N_Blocks; n++){
    // pseudo-code ish
    for(int g=0; g<2; g++){
      if(ID_list[n][g].size()>0){
        filters[n][g][0]=true;
      }
      if(ResID_list[n][g].size()>0){
        filters[n][g][1]=true;
      }
      if(Name_list[n][g].size()>0){
        filters[n][g][2]=true;
      }
    }
  }

  // Not sure if the outer loop should be over n blocks or over i indices
  listall.clear();
  // Need to add code for global atom list
  for(int i; i<mypdb.getAtomNumbers().size(); i++){
    AtomNumber ind=mypdb.getAtomNumbers()[i];
    int resid=mypdb.getResidueNumber(ind);
    string atom_name=mypdb.getAtomName(ind);
    bool atom_added = false;
    for(int n=0; n<N_Blocks; n++){
      for(int g=0; g<2; g++){
        // Boolean to integer conversion
        bool id_check = false;
        bool res_check = false;
        bool name_check = false;
        int caseIndex = filters[n][g][2] * 4 + filters[n][g][1] * 2 + filters[n][g][0];
        switch(caseIndex) {
          case 0: // 000 in binary (false, false, false)
            // An error should be thrown here because it means
            // that a block was defined without any filter specifications
            error("Error: Block" + std::to_string(n) + ", Group" + std::to_string(g) + " has no identifier corresponding to ATOMID, RESID, or NAME");
            break;
          case 1: // 001 in binary (false, false, true)
            for(int f=0; f<ID_list[n][g].size(); f++){
              if(ind==ID_list[n][g][f]){
                block_groups_atom_list[n][g].push_back(ind);
                // Why am I erasing the ID here? Somethings the world may never know...
                ID_list[n][g].erase(ID_list[n][g].begin() + f);
                atom_added = true;
                break;
              }
            }
            break;
          case 2: // 010 in binary (false, true, false)
            for(int f=0; f<ResID_list[n][g].size(); f++){
              if(resid==ResID_list[n][g][f]){
                block_groups_atom_list[n][g].push_back(ind);
                atom_added = true;
                break;
              }
            }
            break;
          case 3: // 011 in binary (false, true, true)
            // This case shouldn't really happen because there would never be a
            // reason to add another identifier if the AtomID is present
            // but is included to make the code more flexible for
            // inexperienced users
            for(int f=0; f<ID_list[n][g].size(); f++){
              if(ind==ID_list[n][g][f]){
                id_check = true;
                break;
              }
            }
            if(id_check){
              block_groups_atom_list[n][g].push_back(ind);
              atom_added = true;
            }
            break;
            // for(int f=0; f<ResID_list[n][g].size(); f++){
            //   if(resid==ResID_list[n][g][f]){
            //     res_check = true;
            //     break;
            //   }
            // }
            // if(res_check && id_check){
            //   block_groups_atom_list[n][g].push_back(ind);
            //   atom_added = true;
            // }
            // break;
            break;
          case 4: // 100 in binary (true, false, false)
            for(int f=0; f<Name_list[n][g].size(); f++){
              if(atom_name==Name_list[n][g][f]){
                block_groups_atom_list[n][g].push_back(ind);
                atom_added = true;
                break;
              }
            }
            break;
          case 5: // 101 in binary (true, false, true)
            // This case shouldn't really happen because there would never be a
            // reason to add another identifier if the AtomID is present
            // but is included to make the code more flexible for
            // inexperienced users
            for(int f=0; f<ID_list[n][g].size(); f++){
              if(ind==ID_list[n][g][f]){
                id_check = true;
                break;
              }
            }
            if(id_check){
              block_groups_atom_list[n][g].push_back(ind);
              atom_added = true;
            }
            break;
            // for(int f=0; f<Name_list[n][g].size(); f++){
            //   if(atom_name==Name_list[n][g][f]){
            //     name_check = true;
            //     break;
            //   }
            // }
            // if(name_check && id_check){
            //   block_groups_atom_list[n][g].push_back(ind);
            //   atom_added = true;
            // }
            // break;
            break;
          case 6: // 110 in binary (true, true, false)
            for(int f=0; f<ResID_list[n][g].size(); f++){
              if(resid==ResID_list[n][g][f]){
                res_check = true;
                break;
              }
            }
            for(int f=0; f<Name_list[n][g].size(); f++){
              if(atom_name==Name_list[n][g][f]){
                name_check = true;
                break;
              }
            }
            if(res_check && name_check){
              block_groups_atom_list[n][g].push_back(ind);
              atom_added = true;
            }
            break;
          case 7: // 111 in binary (true, true, true)
            // This case shouldn't really happen because there would never be a
            // reason to add another identifier if the AtomID is present
            // but is included to make the code more flexible for
            // inexperienced users
            for(int f=0; f<ID_list[n][g].size(); f++){
              if(ind==ID_list[n][g][f]){
                id_check = true;
                break;
              }
            }
            if(id_check){
              block_groups_atom_list[n][g].push_back(ind);
              atom_added = true;
            }
            break;
            // for(int f=0; f<ResID_list[n][g].size(); f++){
            //   if(resid==ResID_list[n][g][f]){
            //     res_check = true;
            //     break;
            //   }
            // }
            // for(int f=0; f<Name_list[n][g].size(); f++){
            //   if(atom_name==Name_list[n][g][f]){
            //     name_check = true;
            //     break;
            //   }
            // }
            // if(res_check && name_check && ind_check){
            //   block_groups_atom_list[n][g].push_back(ind);
            //   atom_added = true;
            // }
            // break;
        }
      }
    }
    if(atom_added){
      listall.push_back(ind);
    }
  }

// Should I have a counter to determine stride updates? If a tolerance update is triggered a step before a stride
// update would have been triggered it is inefficient to immediately update again. The stride update should be the
// stride since the last update, whether from stride or tolerance.


// Should I use AtomNumber or index in block_groups_atom_list to specify atom pairs?
  if(getStep()==0){
  using AtomPair = std::pair<AtomNumber, AtomNumber>;
  using DistAtomPair = std::pair<double, AtomPair>;
  using MaxHeap = std::priority_queue<DistAtomPair, std::vector<DistAtomPair>, CompareDist>;
  using MaxHeapVector = std::vector<MaxHeap>;
  MaxHeapVector maxHeapVec;
    for(int n=0; n<N_Blocks; n++){
      int tot_num_pairs = block_lengths[n] + Buffer_Pairs[n];
      for(int i=0; i<block_groups_atom_list[n][0].size(); i++){
        Vector Pos0,Pos1;
        AtomNumber ind0=block_groups_atom_list[n][0][i];
        Pos0=mypdb.getPosition(ind0);
        for(int j=0; j<block_groups_atom_list[n][1].size(); j++){
          AtomNumber ind1=block_groups_atom_list[n][1][j];
          AtomPair test_pair={ind0,ind1};
          if(Exclude_Pairs[n].find(test_pair) != Exclude_Pairs[n].end()){
            continue;
          }
          Pos1=mypdb.getPosition(ind1);
          ddist=pbcDistance(Pos0,Pos1);
          mag=ddist.modulo();
          if(maxHeapVec[n].size() < tot_num_pairs){
            maxHeapVec[n].push_back({mag, {ind0,ind1}});
          } else if (mag < maxHeapVec[n].top().first) {
            maxHeapVec[n].pop();
            maxHeapVec[n].push_back({mag, {ind0,ind1}});
          }
        }
      }
    }
    int steps_since_update=1;
  }
  if (enableLog) {
    log << "Total Nlists: " << N_Blocks << " \n";
    for (unsigned n=0; n<N_Blocks; n++) {
      log << "  PIV Block " << n+1 << " + Buffer: " << maxHeapVec[n] << " Pairs\n";
    }
  }

  r00.resize(N_Blocks);
  sw.resize(N_Blocks);
  for (unsigned n=0; n<N_Blocks; n++) {
    if( !parseNumbered( "SWITCH", n+1, sw[n] ) ) break;
  }
  if(CompDer) {
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    if (enableLog) {
      log << "Switching Function Parameters \n";
    }
    sfs.resize(N_Blocks);
    std::string errors;
    for (unsigned n=0; n<N_Blocks; n++) {
      sfs[n].set(sw[n],errors);
      std::string num;
      Tools::convert(n+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[n]=sfs[n].get_r0();
      if (enableLog) {
        log << "  Swf: " << n << "  r0=" << (sfs[n].description()).c_str() << " \n";
      }
    }
  }


  // build the rPINES distances (transformation and sorting is done afterwards)
  if (enableLog) {
    if(CompDer) {
      log << "  PINES  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    }
  }
  checkRead();
  // Create components of PINES

  // cPINES hasn't been created yet so it can't be used in the loop.
  // The loop is set up generally as N(N-1)/2. 
  // The if/else statement accounts for there being >1 elements in the solute-solvent blocks, 
  // and expects that the interaction with solvent is the last block of the each solute atom's interactions
  // Total count keeps a running tally of the elements in the entire PINES so that there will be an equal number of components.
  int total_count=0
  for(int n=0; n<N_Blocks; n++){
    for(i=0; i<block_lengths[n]; i++){
      string comp = "ELEMENT-" + to_string(total_count);
      addComponentWithDerivatives(comp); 
      componentIsNotPeriodic(comp);
      total_count += 1;
    }
  }

  using AtomPair = std::pair<AtomNumber, AtomNumber>;
  using DistAtomPair = std::pair<double, AtomPair>;
  using MaxHeap = std::priority_queue<DistAtomPair, std::vector<DistAtomPair>, CompareDist>;
  using MaxHeapVector = std::vector<MaxHeap>;

  std::set<AtomNumber> uniqueIndices;

  for(n=0; n<N_Blocks; n++){
    MaxHeap tempHeap = maxHeapVec[n];
    while(!tempHeap.empty()){
      DistAtomPair topElement = tempHeap.top();
      tempHeap.pop();
      uniqueIndices.insert(topElement.second.first);
      uniqueIndices.insert(topElement.second.second);
    }
  }
  listreduced.clear();
  listreduced.insert(listreduced.end(), uniqueIndices.begin(), uniqueIndices.end());
  requestAtoms(listreduced);

  ann_deriv.resize(listreduced.size());

 // printf("TOTAL COUNT :%d\n\n\n", total_count);
  for (int i = 0; i < ann_deriv.size(); i++) {
    ann_deriv[i].resize(total_count);
  }

  ds_array.resize(total_count);
}

// Does this need to be implemented for the heaps?

// The following deallocates pointers
// PINES::~PINES()
// {
//   for (unsigned j=0; j<Nlist; j++) {
//     delete nl[j];
//   }
//   for(unsigned j=0; j<Nlist; j++) {
//     delete nl_small[j];
//   }
//   delete nlall;
//   delete nlreduced;
// }

// Would it look like this when implemented?

// PINES::~PINES()
// {
//   for (unsigned j=0; j<N_Blocks; j++) {
//     delete maxHeapVec[j];
//   }
// }


// SD request atoms in every frame.
void PINES::prepare() {
  // This should be changed from a (fixed) global perspective for frames
  // To a (dynamic) local perspective where the condition is triggered only nstride since
  // The last update
  // std::vector<std::pair<AtomNumber, AtomNumber> > pairlist(N_Blocks);
  if(getStep()!=0){
    int total_count=0;
    for(int n=0; n<N_Blocks; n++){
      total_count += block_lengths[n];
    }

    if(steps_since_update==0) {

      // delete maxHeapVec; // ???
      using AtomPair = std::pair<AtomNumber, AtomNumber>;
      using DistAtomPair = std::pair<double, AtomPair>;
      using MaxHeap = std::priority_queue<DistAtomPair, std::vector<DistAtomPair>, CompareDist>;
      using MaxHeapVector = std::vector<MaxHeap>;
      MaxHeapVector maxHeapVec;
      for(int n=0; n<N_Blocks; n++){
        int tot_num_pairs = block_lengths[n] + Buffer_Pairs[n];
        for(int i=0; i<block_groups_atom_list[n][0].size(); i++){
          Vector Pos0,Pos1;
          AtomNumber ind0=block_groups_atom_list[n][0][i];
          Pos0=getPosition(ind0);
          for(int j=0; j<block_groups_atom_list[n][1].size(); j++){
            AtomNumber ind1=block_groups_atom_list[n][1][j];
            AtomPair test_pair={ind0,ind1};
            if(Exclude_Pairs[n].find(test_pair) != Exclude_Pairs[n].end()){
              continue;
            }
            Pos1=mypdb.getPosition(ind1);
            ddist=pbcDistance(Pos0,Pos1);
            mag=ddist.modulo();
            if(maxHeapVec[n].size() < tot_num_pairs){
              maxHeapVec[n].push_back({mag, {ind0,ind1}});
            } else if (mag < maxHeapVec[n].top().first) {
              maxHeapVec[n].pop();
              maxHeapVec[n].push_back({mag, {ind0,ind1}});
            }
          }
        }
      }

      std::set<AtomNumber> uniqueIndices;

      for(n=0; n<N_Blocks; n++){
        MaxHeap tempHeap = maxHeapVec[n];
        while(!tempHeap.empty()){
          DistAtomPair topElement = tempHeap.top();
          tempHeap.pop();
          uniqueIndices.insert(topElement.second.first);
          uniqueIndices.insert(topElement.second.second);
        }
      }
      listreduced.clear();
      listreduced.insert(listreduced.end(), uniqueIndices.begin(), uniqueIndices.end());
      requestAtoms(listreduced);

      ann_deriv.resize(listreduced.size());

    // printf("TOTAL COUNT :%d\n\n\n", total_count);
      for (int i = 0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }

      ds_array.resize(total_count);
      steps_since_update += 1;
    } else if(steps_since_update>=nstride) {
      requestAtoms(listall);

      ann_deriv.resize(listall.size());

      int total_count=0
      for(int n=0; n<N_Blocks; n++){
        total_count += block_lengths[n];
      }
      for(unsigned i=0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }
      ds_array.resize(total_count);
      steps_since_update = 0;
    }
  }
}

void PINES::calculate()
{

  // Local variables
  // The following are probably needed as static arrays
  static int prev_stp=-1;
  static int init_stp=1;
  static std:: vector<std:: vector<Vector> > prev_pos(N_Blocks);
  static std:: vector<std:: vector<double> > cPINES(N_Blocks);
  static std:: vector<std:: vector<int> > Atom0(N_Blocks);
  static std:: vector<std:: vector<int> > Atom1(N_Blocks);
  std:: vector<std:: vector<int> > A0(Nprec);
  std:: vector<std:: vector<int> > A1(Nprec);
  //std:: vector<std::pair<double, std::pair<AtomNumber, AtomNumber> > > sortedHeapDict(N_Blocks);
  size_t stride=1;
  unsigned rank=0;

  // Serial vs parallelization?
  if(!serial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  } else {
    stride=1;
    rank=0;
  }

//Tolerance check and update
    // Do the sorting only once per timestep to avoid building the PINES N times for N rPINES PDB structures!
    //
    // build COMs from positions if requested
    // update neighbor lists when an atom moves out of the Neighbor list skin
  if (doneigh && ((getStep()+1)%nlall->getStride()==0)) {
    bool doupdate=false;
    // For the first step build previous positions = actual positions
    if (prev_stp==-1) {
      for (unsigned j=0; j<Nlist; j++) {
        for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
          Vector Pos;

          Pos=getPosition(nl[j]->getFullAtomList()[i].index());

          prev_pos[j].push_back(Pos);
        }
      }
      doupdate=true;
    }
    // Decide whether to update lists based on atom displacement, every stride
    std:: vector<std:: vector<Vector> > tmp_pos(Nlist);
    if (getStep() % nlall->getStride() ==0) {
      for (unsigned j=0; j<Nlist; j++) {
        for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
          Vector Pos;

          Pos=getPosition(nl[j]->getFullAtomList()[i].index());

          tmp_pos[j].push_back(Pos);
          if (pbcDistance(tmp_pos[j][i],prev_pos[j][i]).modulo()>=nl_skin[j]) {
            doupdate=true;
          }
        }
      }
    }
    // Update Nlists if needed
    if (doupdate==true) {
      for (unsigned j=0; j<Nlist; j++) {
        for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
          prev_pos[j][i]=tmp_pos[j][i];
        }
        nl[j]->update(prev_pos[j]);
        if (enableLog) {
          log << " Step " << getStep() << "  Neighbour lists updated " << nl[j]->size() << " \n";
        }
      }
    }
  }
  Vector ddist;

  // Build "Nlist" PINES blocks
  for(unsigned n=0; n<N_Blocks; n++) {
    if(maxHeapVec[n].size()!= 0) {
      int total_num_pairs = block_lengths[n] + Buffer_Pairs[n];
      for(int i=0; i<total_num_pairs; i++) {
        if(i>=Buffer_Pairs[n]){
          int block_ind = total_num_pairs - 1 - i;
          cPINES[n][block_ind] = maxHeapVec[n].top().first;
          pairlist[n][block_ind] = maxHeapVec[n].top().second;
        }
        //sortedHeapDict[n].push_back(maxHeapVec[n]);
        maxHeapVec[n].pop();
      }
    } else {

      
      // Rebuild MaxHeapVec by calculatin dist from pairlist inds
      using AtomPair = std::pair<AtomNumber, AtomNumber>;
      using DistAtomPair = std::pair<double, AtomPair>;
      using MaxHeap = std::priority_queue<DistAtomPair, std::vector<DistAtomPair>, CompareDist>;
      using MaxHeapVector = std::vector<MaxHeap>;
      MaxHeapVector maxHeapVec;

      int tot_num_pairs = block_lengths[n] + Buffer_Pairs[n];
      for(int i=0; i<tot_num_pairs; i++){
        Vector Pos0,Pos1;
        ind0=pairlist[n][i].first;
        ind1=pairlist[n][i].second;
        Pos0=getPosition(ind0);
        Pos1=getPosition(ind1);
        ddist=pbcDistance(Pos0,Pos1);
        mag=ddist.modulo();
        maxHeapVec[n].push_back({mag, {ind0,ind1}});
        //Need to calculate cPINES now?
      }
      for(int i=0; i<tot_num_pairs; i++){
        int block_ind = tot_num_pairs - 1 - i;
        if(i>=Buffer_Pairs[n]){
          cPINES[n][block_ind] = maxHeapVec[n].top().first;
          pairlist[n][block_ind] = maxHeapVec[n].top().second;
          maxHeapVec[n].pop();
        }
      }
      
    }
  }

  Vector distance;
  double dfunc=0.;

 // Create/open PV value files                                              
  FILE *PINES_rep_file = NULL;
  if (writestride) {
    if (getStep() % writePINESstride == 0) {                                                                                      
      string PINES_rep_fileName = "PINES_representation_" + to_string(getStep()) + ".dat";                                    
      PINES_rep_file = fopen(PINES_rep_fileName.c_str(), "w+"); 
    }
  }
  FILE *PINES_rep_file_traj = NULL;
  if (writePINEStraj) {
    string PINES_rep_fileName_traj = "PINES_representation_traj.dat";
    PINES_rep_file_traj = fopen(PINES_rep_fileName_traj.c_str(), "a");
  }
  
// Build ann_deriv
  for(unsigned j=0; j<ann_deriv.size(); j++) {
    for(unsigned i=0; i<ann_deriv[j].size(); i++) {
      for(unsigned k=0; k<3; k++) {ann_deriv[j][i][k]=0.;}
    }
  }
  for(unsigned j=0; j<3; j++) {
    for(unsigned k=0; k<3; k++) {
      m_virial[j][k]=0.;
    }
  }
  for(unsigned n=0; n<N_Blocks; n++) {
    for(unsigned i=0; i<block_lengths[n]; i++) {
      PIV[n][i]=0.;
    }
  }
  // resize vectors to the appropriate sizes and set starting values to zero --NH
  //PINES_Pair0.resize(ds_array.size());
  //PINES_Pair1.resize(ds_array.size());
  unsigned PINES_element=0;
  for(unsigned n=0; n<N_Blocks; n++) {
    for(unsigned i=0; i<block_lengths[n]; i++) {
      unsigned i0=0;
      unsigned i1=0;
      // // Atom0 and Atom1 are lists that index atoms for PINES elements
      i0=pairlist[n][i].first;
      // // Record the atom IDs for the PINES elements of interest --NH
      // PINES_Pair0[PINES_element] = i0;
      i1=pairlist[n][i].second;
      // PINES_Pair1[PINES_element] = i1;
      // Pos0 and Pos1 are 1x3 vectors that hold the xyz coordinates of the indexed atoms
      Vector Pos0,Pos1;
      Pos0=getPosition(i0);
      Pos1=getPosition(i1);
      // // distance is also a 1x3 vector of the xyz distances between the two atoms after consideration of the pbc
      distance=pbcDistance(Pos0,Pos1);
      dfunc=0.;
      // dm is a scalar value that is the magnitude of the distance between the atoms
      double dm=distance.modulo();
      // sfs[j] is the parameters for the switching function, which can be chosen to be different for different blocks
      // In this case, all blocks use the same switching function so all sfs[j] are the same function.
      // Used with .calculate(dm*Fvol, dfunc), the PINES element value is returned and the derivative stored in dfunc
      PIV[n][i] = sfs[n].calculate(cPINES[n][i]*Fvol, dfunc);

      double ds_element=0.;
      // Create the ds_array one element at a time --NH
      ds_element = scaling[n]*Fvol*Fvol*dfunc*cPINES[n][i];
      ds_array[PINES_element] = ds_element;
      // Create 1x3 vector of (dr/dx,dr/dy,dr/dz) --NH
      Vector dr_dcoord = distance/dm;
      
      // the xyz components of the distance between atoms is scaled by tmp and added or subtracted to reflect
      // that distance is calculated as Pos1 - Pos0
      // Calculate ann_deriv values for the current PINES element in the loop --NH
      ann_deriv[i0][PINES_element] = -ds_element*dr_dcoord;
      ann_deriv[i1][PINES_element] =  ds_element*dr_dcoord;

      // This m_virial is likely not correct but has been included in case it is necessary to test the code --NH
      m_virial    -= ds_element*Tensor(distance,distance); // Question
      PINES_element += 1;

    }
  }
  if (!serial && comm.initialized() ) {
    int count = 0;
    for(unsigned j=0; j<N_Blocks; j++) {
        for(unsigned i=0; i<cPINES[j].size(); i++) {
            count += 1;
        }
    }
    
    comm.Barrier();
    // SD -- This probably works because cPINES[j] size is variable for each j.
    for (unsigned j=0; j< N_Blocks; j++) {
      for (unsigned k=0; k<cPINES[j].size(); k++) {
        comm.Sum(cPINES[j][k]);
        cPINES[j][k] /= comm.Get_size();
      }
    }

    // SD -- This probably works because comm.Sum cannot handle 3D vectors.
    if(!ann_deriv.empty()) {
      for (unsigned i=0;  i<ann_deriv.size(); i++) {
        for (unsigned j=0; j<ann_deriv[j].size(); j++) {
          for (unsigned k=0; k<3; k++) {
            comm.Sum(ann_deriv[i][j][k]);
            ann_deriv[i][j][k] /= comm.Get_size();
          }
        }
      }
    }
    // SD -- this is probably not needed.
    comm.Sum(&m_virial[0][0],9);
  }
  prev_stp=getStep();

// PIV calculated, write values to files
  for(unsigned n=0; n<N_Blocks; n++) {
    int total_num_pairs = block_lengths[n] + Buffer_Pairs[n];
    for(int pair=total_num_pairs-1; pair<0; pair--) {
      maxHeapVec[n].pop();
      if(pair<block_lengths[n]){
        // use this pair/dist in PIV block n
        // hash the positional value of the pair in the block with its PIV value
      }
    }
  }
  if (writePINEStraj) {
    fprintf(PINES_rep_file_traj, "\n#END OF FRAME\n");
    fclose(PINES_rep_file_traj);
  }
  if (writestride) {
    if ( getStep() % writePINESstride == 0) {
      fprintf(PINES_rep_file, "\n#END OF FRAME: %lld \n", getStep());
      fclose(PINES_rep_file);
    }
  }


// Pass values and derivates to next stage
  unsigned total_count=0;
  for (unsigned j = 0; j < N_Blocks; j++) {
    // This should really be block_lengths + buffer and probably needs to be
    // corrected elsewhere
    for (unsigned i = 0; i < block_lengths[j]; i++) {
      string comp = "ELEMENT-" + to_string(total_count);
      Value* valueNew=getPntrToComponent(comp);
      valueNew -> set(cPINES[j][i]);
      // Pass the 3D array to the plumed core --NH
      // A 2D array is passed for each PINES element (component) --NH
      for(unsigned k=0; k<ann_deriv.size(); k++) {
        setAtomsDerivatives(valueNew, k, ann_deriv[k][total_count]);
      }
      total_count += 1;
    }
  }

}
//Close Namespaces at the very beginning
}
}
