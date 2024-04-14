//
// Created by genshen on 2019-01-10.
//

#include "period_lattice_list.h"
#include <iostream>
#include <logs/logs.h>
PeriodLatticeList::PeriodLatticeList(const LatListMeta meta) : LatticesList(meta) {}

_type_neighbour_status PeriodLatticeList::get1nnStatus(_type_lattice_coord x, _type_lattice_coord y,
                                                       _type_lattice_coord z) {
  return 0xFF;
}

_type_neighbour_status PeriodLatticeList::get2nnStatus(_type_lattice_coord x, _type_lattice_coord y,
                                                       _type_lattice_coord z) {
  return 0x2F;
}

int PeriodLatticeList::get1nn(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice **_1nn_list) {
  if (x % 2 == 0) {
    // unsigned long latti_id = getId(x - 1, y - 1, z - 1);
    _1nn_list[0] = &_lattices[z - 1][y - 1][x - 1];
    _1nn_list[1] = &_lattices[z + 0][y - 1][x - 1];
    _1nn_list[2] = &_lattices[z - 1][y + 0][x - 1];
    _1nn_list[3] = &_lattices[z + 0][y + 0][x - 1];
    _1nn_list[4] = &_lattices[z - 1][y - 1][x + 1];
    _1nn_list[5] = &_lattices[z + 0][y - 1][x + 1];
    _1nn_list[6] = &_lattices[z - 1][y + 0][x + 1];
    _1nn_list[7] = &_lattices[z + 0][y + 0][x + 1];
  } else {
    _1nn_list[0] = &_lattices[z + 0][y + 0][x - 1];
    _1nn_list[1] = &_lattices[z + 1][y + 0][x - 1];
    _1nn_list[2] = &_lattices[z + 0][y + 1][x - 1];
    _1nn_list[3] = &_lattices[z + 1][y + 1][x - 1];
    _1nn_list[4] = &_lattices[z + 0][y + 0][x + 1];
    _1nn_list[5] = &_lattices[z + 1][y + 0][x + 1];
    _1nn_list[6] = &_lattices[z + 0][y + 1][x + 1];
    _1nn_list[7] = &_lattices[z + 1][y + 1][x + 1];
  }
  return 8;
}

int PeriodLatticeList::get1nn2(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice *_1nn_list) {
  _type_lattice_id latti_id;
  if (x % 2 == 0) {
    latti_id = getId(x - 1, y - 1, z - 1);
    getnnlattice(latti_id, &_1nn_list[0]);
    latti_id = getId(x - 1, y - 1, z);
    getnnlattice(latti_id, &_1nn_list[1]);
    latti_id = getId(x - 1, y, z - 1);
    getnnlattice(latti_id, &_1nn_list[2]);
    latti_id = getId(x - 1, y, z);
    getnnlattice(latti_id, &_1nn_list[3]);
    latti_id = getId(x + 1, y - 1, z - 1);
    getnnlattice(latti_id, &_1nn_list[4]);
    latti_id = getId(x + 1, y - 1, z);
    getnnlattice(latti_id, &_1nn_list[5]);
    latti_id = getId(x + 1, y, z - 1);
    getnnlattice(latti_id, &_1nn_list[6]);
    latti_id = getId(x + 1, y, z);
    getnnlattice(latti_id, &_1nn_list[7]);
    // _1nn_list[0] = &_lattices[z - 1][y - 1][x - 1];
    // _1nn_list[1] = &_lattices[z + 0][y - 1][x - 1];
    // _1nn_list[2] = &_lattices[z - 1][y + 0][x - 1];
    // _1nn_list[3] = &_lattices[z + 0][y + 0][x - 1];
    // _1nn_list[4] = &_lattices[z - 1][y - 1][x + 1];
    // _1nn_list[5] = &_lattices[z + 0][y - 1][x + 1];
    // _1nn_list[6] = &_lattices[z - 1][y + 0][x + 1];
    // _1nn_list[7] = &_lattices[z + 0][y + 0][x + 1];
  } else {
    latti_id = getId(x - 1, y, z);
    getnnlattice(latti_id, &_1nn_list[0]);
    latti_id = getId(x - 1, y, z + 1);
    getnnlattice(latti_id, &_1nn_list[1]);
    latti_id = getId(x - 1, y + 1, z);
    getnnlattice(latti_id, &_1nn_list[2]);
    latti_id = getId(x - 1, y + 1, z + 1);
    getnnlattice(latti_id, &_1nn_list[3]);
    latti_id = getId(x + 1, y, z);
    getnnlattice(latti_id, &_1nn_list[4]);
    latti_id = getId(x + 1, y, z + 1);
    getnnlattice(latti_id, &_1nn_list[5]);
    latti_id = getId(x + 1, y + 1, z);
    getnnlattice(latti_id, &_1nn_list[6]);
    latti_id = getId(x + 1, y + 1, z + 1);
    getnnlattice(latti_id, &_1nn_list[7]);
    // _1nn_list[0] = &_lattices[z + 0][y + 0][x - 1];
    // _1nn_list[1] = &_lattices[z + 1][y + 0][x - 1];
    // _1nn_list[2] = &_lattices[z + 0][y + 1][x - 1];
    // _1nn_list[3] = &_lattices[z + 1][y + 1][x - 1];
    // _1nn_list[4] = &_lattices[z + 0][y + 0][x + 1];
    // _1nn_list[5] = &_lattices[z + 1][y + 0][x + 1];
    // _1nn_list[6] = &_lattices[z + 0][y + 1][x + 1];
    // _1nn_list[7] = &_lattices[z + 1][y + 1][x + 1];
  }
  return 8;
}

int PeriodLatticeList::get1nn3(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice *_1nn_list, _type_lattice_id source_id, LatticeTypes atom_type) {
  _type_lattice_id latti_id;
  if (x % 2 == 0) {
    latti_id = getId(x - 1, y - 1, z - 1);
    if(latti_id == source_id){
      _1nn_list[0].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[0]);
    }
    latti_id = getId(x - 1, y - 1, z);
    if(latti_id == source_id){
      _1nn_list[1].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[1]);
    }
    latti_id = getId(x - 1, y, z - 1);
    if(latti_id == source_id){
      _1nn_list[2].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[2]);
    }
    latti_id = getId(x - 1, y, z);
    if(latti_id == source_id){
      _1nn_list[3].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[3]);
    }
    latti_id = getId(x + 1, y - 1, z - 1);
    if(latti_id == source_id){
      _1nn_list[4].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[4]);
    }
    latti_id = getId(x + 1, y - 1, z);
    if(latti_id == source_id){
      _1nn_list[5].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[5]);
    }
    latti_id = getId(x + 1, y, z - 1);
    if(latti_id == source_id){
      _1nn_list[6].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[6]);
    }
    latti_id = getId(x + 1, y, z);
    if(latti_id == source_id){
      _1nn_list[7].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[7]);
    }
    // _1nn_list[0] = &_lattices[z - 1][y - 1][x - 1];
    // _1nn_list[1] = &_lattices[z + 0][y - 1][x - 1];
    // _1nn_list[2] = &_lattices[z - 1][y + 0][x - 1];
    // _1nn_list[3] = &_lattices[z + 0][y + 0][x - 1];
    // _1nn_list[4] = &_lattices[z - 1][y - 1][x + 1];
    // _1nn_list[5] = &_lattices[z + 0][y - 1][x + 1];
    // _1nn_list[6] = &_lattices[z - 1][y + 0][x + 1];
    // _1nn_list[7] = &_lattices[z + 0][y + 0][x + 1];
  } else {
    latti_id = getId(x - 1, y, z);
    if(latti_id == source_id){
      _1nn_list[0].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[0]);
    }
    latti_id = getId(x - 1, y, z + 1);
    if(latti_id == source_id){
      _1nn_list[1].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[1]);
    }
    latti_id = getId(x - 1, y + 1, z);
    if(latti_id == source_id){
      _1nn_list[2].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[2]);
    }
    latti_id = getId(x - 1, y + 1, z + 1);
    if(latti_id == source_id){
      _1nn_list[3].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[3]);
    }
    latti_id = getId(x + 1, y, z);
    if(latti_id == source_id){
      _1nn_list[4].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[4]);
    }
    latti_id = getId(x + 1, y, z + 1);
    if(latti_id == source_id){
      _1nn_list[5].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[5]);
    }
    latti_id = getId(x + 1, y + 1, z);
    if(latti_id == source_id){
      _1nn_list[6].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[6]);
    }
    latti_id = getId(x + 1, y + 1, z + 1);
    if(latti_id == source_id){
      _1nn_list[7].type = atom_type;
    }else{
      getnnlattice(latti_id, &_1nn_list[7]);
    }
    // _1nn_list[0] = &_lattices[z + 0][y + 0][x - 1];
    // _1nn_list[1] = &_lattices[z + 1][y + 0][x - 1];
    // _1nn_list[2] = &_lattices[z + 0][y + 1][x - 1];
    // _1nn_list[3] = &_lattices[z + 1][y + 1][x - 1];
    // _1nn_list[4] = &_lattices[z + 0][y + 0][x + 1];
    // _1nn_list[5] = &_lattices[z + 1][y + 0][x + 1];
    // _1nn_list[6] = &_lattices[z + 0][y + 1][x + 1];
    // _1nn_list[7] = &_lattices[z + 1][y + 1][x + 1];
  }
  return 8;
}

int PeriodLatticeList::get2nn(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice **_2nn_list) {
  _2nn_list[0] = &_lattices[z + 0][y + 0][x - 2];
  _2nn_list[1] = &_lattices[z + 0][y - 1][x + 0];
  _2nn_list[2] = &_lattices[z - 1][y + 0][x + 0];
  _2nn_list[3] = &_lattices[z + 1][y + 0][x + 0];
  _2nn_list[4] = &_lattices[z + 0][y + 1][x + 0];
  _2nn_list[5] = &_lattices[z + 0][y + 0][x + 2];
  return 6;
}

int PeriodLatticeList::get2nn2(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice *_2nn_list) {
  _type_lattice_id latti_id;
  latti_id = getId(x - 2, y, z);
  getnnlattice(latti_id, &_2nn_list[0]);
  // _2nn_list[0] = &_lattices[z + 0][y + 0][x - 2];
  latti_id = getId(x, y - 1, z);
  getnnlattice(latti_id, &_2nn_list[1]);
  // _2nn_list[1] = &_lattices[z + 0][y - 1][x + 0];
  latti_id = getId(x, y, z - 1);
  getnnlattice(latti_id, &_2nn_list[2]);
  // _2nn_list[2] = &_lattices[z - 1][y + 0][x + 0];
  latti_id = getId(x, y, z + 1);
  getnnlattice(latti_id, &_2nn_list[3]);
  // _2nn_list[3] = &_lattices[z + 1][y + 0][x + 0];
  latti_id = getId(x, y + 1, z);
  getnnlattice(latti_id, &_2nn_list[4]);
  // _2nn_list[4] = &_lattices[z + 0][y + 1][x + 0];
  latti_id = getId(x + 2, y, z);
  getnnlattice(latti_id, &_2nn_list[5]);
  // _2nn_list[5] = &_lattices[z + 0][y + 0][x + 2];
  return 6;
}

int PeriodLatticeList::get2nn3(_type_lattice_coord x, _type_lattice_coord y, _type_lattice_coord z,
                              Lattice *_2nn_list, _type_lattice_id source_id, LatticeTypes atom_type) {
  _type_lattice_id latti_id;
  
  latti_id = getId(x - 2, y, z);
  if(latti_id == source_id){
    _2nn_list[0].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[0]);
  }
  // _2nn_list[0] = &_lattices[z + 0][y + 0][x - 2];
  latti_id = getId(x, y - 1, z);
  if(latti_id == source_id){
    _2nn_list[1].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[1]);
  }
  // _2nn_list[1] = &_lattices[z + 0][y - 1][x + 0];
  latti_id = getId(x, y, z - 1);
  if(latti_id == source_id){
    _2nn_list[2].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[2]);
  }
  // _2nn_list[2] = &_lattices[z - 1][y + 0][x + 0];
  latti_id = getId(x, y, z + 1);
  if(latti_id == source_id){
    _2nn_list[3].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[3]);
  }
  // _2nn_list[3] = &_lattices[z + 1][y + 0][x + 0];
  latti_id = getId(x, y + 1, z);
  if(latti_id == source_id){
    _2nn_list[4].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[4]);
  }
  // _2nn_list[4] = &_lattices[z + 0][y + 1][x + 0];
  latti_id = getId(x + 2, y, z);
  if(latti_id == source_id){
    _2nn_list[5].type = atom_type;
  }else{
    getnnlattice(latti_id, &_2nn_list[5]);
  }
  // _2nn_list[5] = &_lattices[z + 0][y + 0][x + 2];
  return 6;
}

void PeriodLatticeList::getnnlattice(_type_lattice_id latti_id, Lattice *_1nn_list_lattice){
  _1nn_list_lattice->id = latti_id;
  auto it_vac = vac_hash.find(latti_id);
  auto it_cu = cu_hash.find(latti_id);
  auto it_mn = mn_hash.find(latti_id);
  auto it_ni = ni_hash.find(latti_id);
  if (it_vac != vac_hash.end()) {
    _1nn_list_lattice->type = LatticeTypes{LatticeTypes::V};
  }else if (it_cu != cu_hash.end()){
    _1nn_list_lattice->type = LatticeTypes{LatticeTypes::Cu};
  }else if (it_mn != mn_hash.end()){
    _1nn_list_lattice->type = LatticeTypes{LatticeTypes::Mn};
  }else if (it_ni != ni_hash.end()){
    _1nn_list_lattice->type = LatticeTypes{LatticeTypes::Ni};
  }else{
    _1nn_list_lattice->type = LatticeTypes{LatticeTypes::Fe};
  }
}
