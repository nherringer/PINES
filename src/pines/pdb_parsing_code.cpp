  
  
  
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

// Parse blocks for keywords
for(int n=0; n<N_Blocks; n++){
  int block_length;
  int limit_per_G1;
  std::vector<int> ex_pairs_n;
  std::vector<int> g1_ids;
  std::vector<int> g2_ids;
  std::vector<int> g1_resids;
  std::vector<int> g2_resids;
  vector<string> block_data=Tools::getWords(block_params[n]);
  vector<string> G1_data=Tools::getWords(block_data[0]);
  vector<string> G2_data=Tools::getWords(block_data[1]);
  parseVector(G1_data,"ATOMID",g1_ids);
  for(int i=0; i<g1_ids.size(); i++){
    ID_list[n][0].push_back(g1_ids[i]);
  }
  parseVector(G2_data,"ATOMID",g2_ids);
  for(int i=0; i<g2_ids.size(); i++){
    ID_list[n][1].push_back(g2_ids[i]));
  }
  parseVector(G1_data,"RESID",g1_resids);
  for(int i=0; i<g1_resids.size(); i++){
    ResID_list[n][0].push_back(g1_resids[i]));
  }
  parseVector(G2_data,"RESID",g2_resids);
  for(int i=0; i<g2_resids.size(); i++){
    ResID_list[n][1].push_back(g2_resids[i]));
  }
  parseVector(G1_data,"NAME",g1_names);
  for(int i=0; i<g1_names.size(); i++){
    Name_list[n][0].push_back(g1_names[i]);
  }
  parseVector(G2_data,"NAME",g2_names);
  for(int i=0; i<g2_names.size(); i++){
    Name_list[n][1].push_back(g2_names[i]);
  }
  Tools::parse(block_data, "SIZE", block_length);
  block_lengths[n]=block_length;
  parseVector(block_data,"EXCLUDE_PAIRS",ex_pairs_n);
  for(int i=0; i<ex_pairs_n.size(); i++){
    Exclude_Pairs[n].push_back(ex_pairs_n[i]);
  }
  Tools::parse(block_data,"LIMIT_PER_G1",limit_per_G1);
  if(limit_per_G1.size()>0){
    G1_limits[n]=limit_per_G1;
  }
}

// In previous section build filter lists to use in building atom lists
// ID_list, Resid_list, and Name_list are all vectors of dimension N_blocks x 2 (G1/G2) x n components
// The union of filters is then easy by checking for each block if more than one list has entries

// The parser will almost certainly require more work but for now I'm moving
// on with the assumption that I've been able to parse the plumed file
// and now I need to generate nl_all and nl_blocks efficiently
std::vector<bool> filters[N_Blocks][2][3]; //N blocks, two groups, three filters
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
      int caseIndex = filters[n][g][2] * 4 + filters[n][g][1] * 2 + filters[n][g][0];
      switch(caseIndex) {
        case 0: // 000 in binary (false, false, false)
          // An error should be thrown here because it means
          // that a block was defined without any filter specifications
          error("Error: Block" + n + ", Group" + g + " has no identifier corresponding to ATOMID, RESID, or NAME");
          break;
        case 1: // 001 in binary (false, false, true)
          for(int f=0; f<ID_list[n][g].size(); f++){
            if(ind==ID_list[n][g][f]){
              block_groups_atom_list[n][g].push_back(ind);
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
          bool id_check = false;
          bool res_check = false;
          bool name_check = false;
          for(int f=0; f<ID_list[n][g].size(); f++){
            if(id==ID_list[n][g][f]){
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
          bool id_check = false;
          bool res_check = false;
          bool name_check = false;
          for(int f=0; f<ID_list[n][g].size(); f++){
            if(id==ID_list[n][g][f]){
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
          bool res_check = false;
          bool name_check = false;
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
          bool id_check = false;
          bool res_check = false;
          bool name_check = false;
          for(int f=0; f<ID_list[n][g].size(); f++){
            if(id==ID_list[n][g][f]){
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

// // n blocks x group x components
// // This currently doesn't take into account multiple filters but this can
// // be added with a switch statement I beleive
// for(int i; i<mypdb.getAtomNumbers().size(); i++){
//   AtomNumber ind=mypdb.getAtomNumbers()[i];
//   int resid=mypdb.getResidueNumber(ind);
//   string atom_name=mypdb.getAtomName(ind);
//   for(int n=0; n<N_Blocks; n++){
//     bool g1_present=false;
//     bool g2_present=false;
//     for(int j=0; j<ID_list[n][0].size(); j++){
//       if(ind==ID_list[n][0][j]){
//         if(!g1_present){
//           G1_atom_list[n].push_back(ind);
//           g1_present=true;
//           break;
//         }
//       }
//       if(ind==ID_list[n][1][j]){
//         if(!g2_present){
//           G2_atom_list[n].push_back(ind);
//           g2_present=true;
//           break;
//         }
//       }
//     }
//     for(int j=0; j<ResID_list[n][0].size(); j++){
//       if(resid==ResID_list[n][0][j]){
//         if(!g1_present){
//           G1_atom_list[n].push_back(ind);
//           g1_present=true;
//           break;
//         }
//       }
//       if(resid==ResID_list[n][1][j]){
//         if(!g2_present){
//           G2_atom_list[n].push_back(ind);
//           g2_present=true;
//           break;
//         }
//       }
//     }
//     for(int j=0; j<Name_list[n][0].size(); j++){
//       if(atom_name==Name_list[n][0][j]){
//         if(!g1_present){
//           G1_atom_list[n].push_back(ind);
//           g1_present=true;
//           break;
//         }
//       }
//       if(atom_name==Name_list[n][1][j]){
//         if(!g2_present){
//           G2_atom_list[n].push_back(ind);
//           g2_present=true;
//           break;
//         }
//       }
//     }
//     if(g1_present || g2_present){
//       if(std::find(listall.begin(), listall.end(), ind) == listall.end()){
//         listall.push_back(ind);
//       }
//     }
//   }
// }

parseVector("NL_CUTOFF",nl_cut);
parseVector("NL_STRIDE",nl_st);
parseVector("NL_SKIN",nl_skin);
for (unsigned n=0; n<N_Blocks; n++) {
  if(nl_cut[n]<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
  if(nl_st[n]<=0) error("NL_STRIDE should be explicitly specified and positive");
  if(nl_skin[n]<=0.) error("NL_SKIN should be explicitly specified and positive");
  nl_cut[n]=nl_cut[n]+nl_skin[n];
}

nlall= new NeighborList(listall,true,pbc,getPbc(),comm,9999999,0);
// build heap neighbor lists
// heap_n-restricted block_groups_atom_list[n][i]

// N closest interactions, so take N

for(int n=0; n<N_Blocks; n++){
  for(int g1=0; g1<block_groups_atom_list[n][0].size(); g1++){
      int ind0 = block_groups_atom_list[n][0][g1]
      Pos0=getPosition(ind0);
    for(int g2=0; g2<block_groups_atom_list[n][1].size(); g2++){
      int ind1 = block_groups_atom_list[n][1][g2]
      Pos1=getPosition(ind1);
      dist=pbcDistance(Pos0,Pos1);
      heap_nl[n].push({ind0, ind1, dist})
    }
  }
  nl[n]= new NeighborList(block_groups_atom_list[n][0],block_groups_atom_list[n][1],true,false,pbc,getPbc(),comm,nl_cut[n],nl_st[n]);
}
// nlall is the full atomlist of all possible relevant atoms
// nl is the smaller atomlist (built from nlall info) of all interactions within a given cutoff distance + buffer distance
// nl_small is the smallest atomlist (built from nl info) of all interactions within a given proximal number + buffer number



void PINES::prepare() {
  if(nlall->getStride()>0) {
    if((getStep()+1)%nlall->getStride()==0) {
      requestAtoms(nlall->getFullAtomList());

      ann_deriv.resize(listall.size());

      int total_count=0;
      for(int i=0; i<block_lengths.size(); i++){
        total_count+=block_lengths[i];
      }
      for(unsigned i=0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }
      ds_array.resize(total_count);
    } else if((getStep()%nlall->getStride()==0) && (getStep()!=0)) {
      
      //initiate lists
      for(int n=0; n<N_Blocks; n++){
        for(int m=0; m<heap_nl[n][0].size(); m++){
          for(int l=0; l<heap_nl[n][1].size(); l++){
              Pos0=getPosition(id0);
              Pos1=getPosition(id1);

              ddist=pbcDistance(Pos0,Pos1);
          }
        }  
      }




      std::vector<vector<AtomNumber>> snn_list;
      std::vector<vector<AtomNumber>> Hydrogen_nn_list;
      std::vector<double> snn_mags;
      std::vector<AtomNumber> all_ids;
      std::vector<double> all_mags;
      std::vector<AtomNumber> total_waterOx_list;
      std::vector<AtomNumber> total_waterHydrogen_list;

      int buffer=10;
      int Nlist_count;
      snn_list.clear();
      Hydrogen_nn_list.clear();
      total_waterOx_list.clear();
      total_waterHydrogen_list.clear();

      if (solv_blocks == 3) {
        Hydrogen_nn_list.resize(Natm-2);
        snn_list.resize(Natm-2);
      } else {
        Hydrogen_nn_list.resize(Natm-1);
        snn_list.resize(Natm-1);
      }

      Vector ddist;
      double smallest_val;
      Nlist_count=0;

      Nlist_count=0;
      for(unsigned j=0; j<Natm; j++) {
        for(unsigned i=j+1; i<Natm; i++) {
          if( (atype[i] == "OW") || ((solv_blocks == 2) && (atype[i] == "HW1")) ) {
            all_mags.clear();
            all_ids.clear();

            for(unsigned k=0; k<nl[Nlist_count]->size(); k++) {
              unsigned id0=(nl[Nlist_count]->getClosePairAtomNumber(k).first).index();
              unsigned id1=(nl[Nlist_count]->getClosePairAtomNumber(k).second).index();
              Vector Pos0,Pos1;
              double mag;
              
              Pos0=getPosition(id0);
              Pos1=getPosition(id1);

              ddist=pbcDistance(Pos0,Pos1);
              mag=ddist.modulo();

              all_mags.push_back(mag);
              all_ids.push_back(listall[id0]);
            }
            for(unsigned x=0; x<NL_const_size+10; x++) {
              smallest_val = all_mags[0];
              int nl_pos = 0;
              for(unsigned y=0; y<all_mags.size(); y++) {
                if(smallest_val > all_mags[y]) {
                  smallest_val = all_mags[y];
                  nl_pos = y;
                }
              }

              snn_mags.push_back(smallest_val);
              snn_list[j].push_back(all_ids[nl_pos]);

              unsigned atom_indx = listall.size();
              if (solv_blocks > 1) {
                for(unsigned k=0; k<listall.size(); k++) {
                  if (listall[k] == all_ids[nl_pos]) {
                    atom_indx = k;
                  }
                }
                for(unsigned k=0; k<listall.size(); k++) {
                  if (AtomToResID_Dict[atom_indx] == AtomToResID_Dict[k]) {
                    if( (solv_blocks == 3) && (atom_indx != k) ){
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    } else if ( (solv_blocks == 2) ) {
                      Hydrogen_nn_list[j].push_back(listall[k]);
                    }
                  }
                }
              }


              all_mags[nl_pos] = 99999.;
              if(total_waterOx_list.empty()) {
                total_waterOx_list.push_back(all_ids[nl_pos]);
              } else {
                bool id_present=false;
                for(unsigned y=0; y<total_waterOx_list.size(); y++) {
                  if(all_ids[nl_pos] == total_waterOx_list[y]) {
                    id_present=true;
                  }
                }
                if(!id_present) {
                  total_waterOx_list.push_back(all_ids[nl_pos]);
                }
              }
              if (solv_blocks > 1) {
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-2]);
                }
                if (std::find(total_waterHydrogen_list.begin(), total_waterHydrogen_list.end(), Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]) == total_waterHydrogen_list.end()) {
                  total_waterHydrogen_list.push_back(Hydrogen_nn_list[j][Hydrogen_nn_list[j].size()-1]);
                }
              }
            }
          }
          Nlist_count += 1;
        }
      }
      if (solv_blocks > 0) {
        std::sort(total_waterOx_list.begin(), total_waterOx_list.end());
      }
      listreduced.clear();
      listreduced = listnonwater;

      if (solv_blocks == 1) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
        }
      } else if (solv_blocks == 2) {
        for(unsigned i=0; i<total_waterHydrogen_list.size(); i++) {
          listreduced.push_back(total_waterHydrogen_list[i]);
        }
      } else if (solv_blocks == 3) {
        for(unsigned i=0; i<total_waterOx_list.size(); i++) {
          listreduced.push_back(total_waterOx_list[i]);
          listreduced.push_back(total_waterHydrogen_list[i*2]);
          listreduced.push_back(total_waterHydrogen_list[i*2+1]);
        }
      }
  
      if (solv_blocks == 3) { 
        for(unsigned j=0; j<Natm-2; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if(snn_list[j][k]==listreduced[i]){
                AtomNumber atom_id;
                snn_list[j][k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
              }
              if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                AtomNumber atom_id;
                Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
              }
            }
          }
        }
      } else {
        for(unsigned j=0; j<Natm-1; j++) {
          for(unsigned k=0; k<NL_const_size+10; k++) {
            for(unsigned i=0; i<listreduced.size(); i++) {
              if (solv_blocks == 1) {
                if(snn_list[j][k]==listreduced[i]){
                  AtomNumber atom_id;
                  snn_list[j][k]=atom_id.setIndex(i);
                }
              }
              if (solv_blocks == 2) {
                if(Hydrogen_nn_list[j][2*k]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k]=atom_id.setIndex(i);
                }
                if(Hydrogen_nn_list[j][2*k+1]==listreduced[i]){
                  AtomNumber atom_id;
                  Hydrogen_nn_list[j][2*k+1]=atom_id.setIndex(i);
                }
              }
            }
          }
        }
      }
      unsigned count_nl_loop=0;
      for(unsigned j=0; j<Natm-1; j++) {
        for(unsigned i=j+1; i<Natm; i++) {
          if (count_nl_loop < Nlist) {
            // Only update solvent blocks (included IDs are dynamic)
            if(atype[i] == "OW") {
              delete nl_small[count_nl_loop];
              nl_small[count_nl_loop] = new NeighborList(snn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
            } else if (atype[i] == "HW1") {
              delete nl_small[count_nl_loop];
              nl_small[count_nl_loop] = new NeighborList(Hydrogen_nn_list[j],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[count_nl_loop],nl_st[count_nl_loop]);
            }
          }
          count_nl_loop += 1;
        }
      }

      
      delete nlreduced;
      nlreduced= new NeighborList(listreduced,true,pbc,getPbc(),comm,99999,nl_st[0]);

      requestAtoms(nlreduced->getFullAtomList());

      ann_deriv.resize(listreduced.size());

     // printf("ANN size %d",ann_deriv.size());
      int total_count=0;
      for(int i=0; i<block_lengths.size(); i++){
        total_count+=block_lengths[i];
      }
      for(unsigned i=0; i < ann_deriv.size(); i++) {
        ann_deriv[i].resize(total_count);
      }
      ds_array.resize(total_count);
    }
  }
}