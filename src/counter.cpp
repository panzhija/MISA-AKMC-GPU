//
// Created by genshen on 2019-03-08.
//

#include "counter.h"
#include "type_define.h"
#include <cassert>
#include <iostream>

void counter::setLatTypeToStrFunc(const fn_lat_type_to_str fn_lat_type_to_str) {
  this->func_lat_type_to_str = fn_lat_type_to_str;
}

int counter::getAtomCount(const LatticeTypes::lat_type tp) {
#ifdef KMC_DEBUG_MODE
  assert(LatticeTypes{tp}.isAtom() || LatticeTypes{tp}.isVacancy());
#endif
  if (LatticeTypes{tp}.isVacancy()) {
    return data[tp];
  }
  int count = 0;
  LatticeTypes::lat_type high, low;
  for (const auto pair : data) {
    high = LatticeTypes{pair.first}.getHighEnd();
    low = LatticeTypes{pair.first}.getLowEnd();
    if (high == tp) {
      count += pair.second;
    }
    if (low == tp) {
      count += pair.second;
    }
  }
  return count;
}

counter counter::newCounter(LatticesList *p_list) {
  counter c;
  for (_type_lattice_coord z = 0; z < p_list->meta.box_z; z++) {
    for (_type_lattice_coord y = 0; y < p_list->meta.box_y; y++) {
      // note: x is already doubled.
      for (_type_lattice_coord x = 0; x < p_list->meta.box_x; x++) {
        _type_lattice_coord latti_id = p_list->getId(p_list->meta.ghost_x + x, p_list->meta.ghost_y + y, p_list->meta.ghost_z + z);
        auto it_vac = p_list->vac_hash.find(latti_id);
        auto it_cu = p_list->cu_hash.find(latti_id);
        auto it_mn = p_list->mn_hash.find(latti_id);
        auto it_ni = p_list->ni_hash.find(latti_id);
        if (it_vac != p_list->vac_hash.end()) {
          c.add(LatticeTypes::V);
        }else if (it_cu != p_list->cu_hash.end()){
          c.add(LatticeTypes::Cu);
        }else if (it_mn != p_list->mn_hash.end()){
          c.add(LatticeTypes::Mn);
        }else if (it_ni != p_list->ni_hash.end()){
          c.add(LatticeTypes::Ni);
        }else{
          c.add(LatticeTypes::Fe);
        }
      }
    }
  }
  return c;
}

std::ostream &operator<<(std::ostream &os, const counter &counter) {
  if (counter.func_lat_type_to_str) {
    for (const auto &c : counter.data) {
      os << "[" << counter.func_lat_type_to_str(c.first) << "]: " << c.second << "\n";
    }
  } else {
    for (const auto &c : counter.data) {
      os << "[" << c.first << "]: " << c.second << "\n";
    }
  }
  return os;
}