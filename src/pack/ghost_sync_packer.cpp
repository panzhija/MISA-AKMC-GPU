//
// Created by genshen on 2019/10/14.
//

#include "ghost_sync_packer.h"
#include "utils/macros.h"
#include <iostream>
#include <omp.h>

GhostSyncPacker::GhostSyncPacker(LatticesList *lats_list) : lats(lats_list) {}

const unsigned long GhostSyncPacker::sendLength(const std::vector<comm::Region<pack_region_type>> send_regions,
                                                const int dimension, const int direction) {
  unsigned long size_ghost = 0;
  for (auto r : send_regions) {
    size_ghost += r.volume();
  }
  return BCC_DBX * size_ghost;
}

const int GhostSyncPacker::sendLength2(const int dimension,
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                                       std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {
  int size_send = 0;
  // for (auto r : send_regions) {
  //   // 返回体积 既发送总量
  //   size_send += r.volume();
  // }
  // return BCC_DBX * size_send;
  switch (dimension)
  {
  case 0:
    size_send = now_exchange_surface_x.size();
    break;
  case 1:
    size_send = now_exchange_surface_y.size();
    break;
  case 2:
    size_send = now_exchange_surface_z.size();
    break;
  default:
    assert(false);
    break;
  }
  return size_send;
}

void GhostSyncPacker::onSend(pack_date_type *buffer, const std::vector<comm::Region<pack_region_type>> send_regions,
                             const unsigned long send_len, const int dimension, const int direction) {
  unsigned long len = 0;
  for (auto &r : send_regions) {
    for (int z = r.z_low; z < r.z_high; z++) {
      for (int y = r.y_low; y < r.y_high; y++) {
        for (int x = r.x_low; x < r.x_high; x++) {
          unsigned long latti_id = lats->getId(BCC_DBX * x, y, z);
          buffer[len].id = latti_id;
          auto it_vac = lats->vac_hash.find(latti_id);
          auto it_cu = lats->cu_hash.find(latti_id);
          auto it_mn = lats->mn_hash.find(latti_id);
          auto it_ni = lats->ni_hash.find(latti_id);
          if (it_vac != lats->vac_hash.end()) {
            buffer[len].type = LatticeTypes{LatticeTypes::V};
          }else if (it_cu != lats->cu_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Cu};
          }else if (it_mn != lats->mn_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Mn};
          }else if (it_ni != lats->ni_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Ni};
          }else{
            buffer[len].type = LatticeTypes{LatticeTypes::Fe};
          }
          len++;

          unsigned long latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          buffer[len].id = latti_id2;
          auto it_vac2 = lats->vac_hash.find(latti_id2);
          auto it_cu2 = lats->cu_hash.find(latti_id2);
          auto it_mn2 = lats->mn_hash.find(latti_id2);
          auto it_ni2 = lats->ni_hash.find(latti_id2);
          if (it_vac2 != lats->vac_hash.end()) {
            buffer[len].type = LatticeTypes{LatticeTypes::V};
          }else if (it_cu2 != lats->cu_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Cu};
          }else if (it_mn2 != lats->mn_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Mn};
          }else if (it_ni2 != lats->ni_hash.end()){
            buffer[len].type = LatticeTypes{LatticeTypes::Ni};
          }else{
            buffer[len].type = LatticeTypes{LatticeTypes::Fe};
          }
          // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x];
          // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x + 1];
          len++;
        }
      }
    }
  }
}

void GhostSyncPacker::onSend2(ChangeLattice *buffer, const int send_len, const int dimension,
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                              std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {

  int len = 0;

  if(dimension == 0) {
    for (const auto& pair : now_exchange_surface_x) {
      _type_lattice_size latti_x = pair % lats->meta.size_x;
      _type_lattice_size latti_y = (pair / lats->meta.size_x) % lats->meta.size_y;
      _type_lattice_size latti_z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].x = latti_x;
      buffer[len].y = latti_y;
      buffer[len].z = latti_z;
      auto it_vac = lats->vac_hash.find(pair);
      auto it_cu = lats->cu_hash.find(pair);
      auto it_mn = lats->mn_hash.find(pair);
      auto it_ni = lats->ni_hash.find(pair);
      if (it_vac != lats->vac_hash.end()) {
        buffer[len].type = LatticeTypes{LatticeTypes::V};
      }else if (it_cu != lats->cu_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Cu};
      }else if (it_mn != lats->mn_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Mn};
      }else if (it_ni != lats->ni_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Ni};
      }else{
        buffer[len].type = LatticeTypes{LatticeTypes::Fe};
      }
      len++;
    }
  }else if(dimension == 1) {
    for (const auto& pair : now_exchange_surface_y) {
      _type_lattice_size latti_x = pair % lats->meta.size_x;
      _type_lattice_size latti_y = (pair / lats->meta.size_x) % lats->meta.size_y;
      _type_lattice_size latti_z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].x = latti_x;
      buffer[len].y = latti_y;
      buffer[len].z = latti_z;
      auto it_vac = lats->vac_hash.find(pair);
      auto it_cu = lats->cu_hash.find(pair);
      auto it_mn = lats->mn_hash.find(pair);
      auto it_ni = lats->ni_hash.find(pair);
      if (it_vac != lats->vac_hash.end()) {
        buffer[len].type = LatticeTypes{LatticeTypes::V};
      }else if (it_cu != lats->cu_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Cu};
      }else if (it_mn != lats->mn_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Mn};
      }else if (it_ni != lats->ni_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Ni};
      }else{
        buffer[len].type = LatticeTypes{LatticeTypes::Fe};
      }
      len++;
    }
  }else {
    for (const auto& pair : now_exchange_surface_z) {
      _type_lattice_size latti_x = pair % lats->meta.size_x;
      _type_lattice_size latti_y = (pair / lats->meta.size_x) % lats->meta.size_y;
      _type_lattice_size latti_z = pair / (lats->meta.size_x * lats->meta.size_y);
      buffer[len].x = latti_x;
      buffer[len].y = latti_y;
      buffer[len].z = latti_z;
      auto it_vac = lats->vac_hash.find(pair);
      auto it_cu = lats->cu_hash.find(pair);
      auto it_mn = lats->mn_hash.find(pair);
      auto it_ni = lats->ni_hash.find(pair);
      if (it_vac != lats->vac_hash.end()) {
        buffer[len].type = LatticeTypes{LatticeTypes::V};
      }else if (it_cu != lats->cu_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Cu};
      }else if (it_mn != lats->mn_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Mn};
      }else if (it_ni != lats->ni_hash.end()){
        buffer[len].type = LatticeTypes{LatticeTypes::Ni};
      }else{
        buffer[len].type = LatticeTypes{LatticeTypes::Fe};
      }
      len++;
    }
  }
  assert(len == send_len);
}

void GhostSyncPacker::onReceive(pack_date_type *buffer, const std::vector<comm::Region<pack_region_type>> recv_regions,
                                const unsigned long receive_len, const int dimension, const int direction) {
  unsigned long len = 0;
  for (auto &r : recv_regions) {
    for (int z = r.z_low; z < r.z_high; z++) {
      for (int y = r.y_low; y < r.y_high; y++) {
        for (int x = r.x_low; x < r.x_high; x++) {
          unsigned long latti_id = lats->getId(BCC_DBX * x, y, z);
          auto it_vac = lats->vac_hash.find(latti_id);
          auto it_cu = lats->cu_hash.find(latti_id);
          auto it_mn = lats->mn_hash.find(latti_id);
          auto it_ni = lats->ni_hash.find(latti_id);
          if (it_vac != lats->vac_hash.end()) {
            lats->vac_hash.erase(it_vac);
          }else if (it_cu != lats->cu_hash.end()){
            lats->cu_hash.erase(it_cu);
          }else if (it_mn != lats->mn_hash.end()){
            lats->mn_hash.erase(it_mn);
          }else if (it_ni != lats->ni_hash.end()){
            lats->ni_hash.erase(it_ni);
          }

          switch (buffer[len].type._type)
          {
          case LatticeTypes::V:
            lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
            // lats->vac_hash.emplace(latti_id);
            break;
          case LatticeTypes::Cu:
            lats->cu_hash.emplace(latti_id);
            break;
          case LatticeTypes::Ni:
            lats->ni_hash.emplace(latti_id);
            break;
          case LatticeTypes::Mn:
            lats->mn_hash.emplace(latti_id);
          default:
            break;
          }
          len++;

          unsigned long latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          auto it_vac2 = lats->vac_hash.find(latti_id2);
          auto it_cu2 = lats->cu_hash.find(latti_id2);
          auto it_mn2 = lats->mn_hash.find(latti_id2);
          auto it_ni2 = lats->ni_hash.find(latti_id2);
          if (it_vac2 != lats->vac_hash.end()) {
            lats->vac_hash.erase(it_vac2);
          }else if (it_cu2 != lats->cu_hash.end()){
            lats->cu_hash.erase(it_cu2);
          }else if (it_mn2 != lats->mn_hash.end()){
            lats->mn_hash.erase(it_mn2);
          }else if (it_ni2 != lats->ni_hash.end()){
            lats->ni_hash.erase(it_ni2);
          }

          switch (buffer[len].type._type)
          {
          case LatticeTypes::V:
            lats->vac_hash.emplace(std::make_pair(latti_id2, VacancyHash{}));
            // lats->vac_hash.emplace(latti_id2);
            break;
          case LatticeTypes::Cu:
            lats->cu_hash.emplace(latti_id2);
            break;
          case LatticeTypes::Ni:
            lats->ni_hash.emplace(latti_id2);
            break;
          case LatticeTypes::Mn:
            lats->mn_hash.emplace(latti_id2);
          default:
            break;
          }
          // lats->_lattices[z][y][BCC_DBX * x].type = buffer[len++].type;
          // lats->_lattices[z][y][BCC_DBX * x + 1].type = buffer[len++].type;
          len++;
        }
      }
    }
  }
}


void GhostSyncPacker::onReceive2(ChangeLattice *buffer, const int receive_len, const int dimension,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                                 const _type_lattice_count *sub_box_lattice_size, 
                                 const _type_lattice_count *neighbour_local_sub_box, unsigned int next_id,
                                 const comm::ColoredDomain *p_domain) {
  
  int len = 0;

  _type_lattice_id x, y, z;
  _type_lattice_id latti_id;
  if(dimension == 2){
    for(len = 0; len < receive_len; len++){
      x = buffer[len].x;
      y = buffer[len].y;
      if(next_id > 3){
        buffer[len].z += sub_box_lattice_size[dimension];
      }else{
        buffer[len].z -= sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]);
      }
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      auto it_vac = lats->vac_hash.find(latti_id);
      auto it_cu = lats->cu_hash.find(latti_id);
      auto it_mn = lats->mn_hash.find(latti_id);
      auto it_ni = lats->ni_hash.find(latti_id);
      if (it_vac != lats->vac_hash.end()) {
        lats->vac_hash.erase(it_vac);
      }else if (it_cu != lats->cu_hash.end()){
        lats->cu_hash.erase(it_cu);
      }else if (it_mn != lats->mn_hash.end()){
        lats->mn_hash.erase(it_mn);
      }else if (it_ni != lats->ni_hash.end()){
        lats->ni_hash.erase(it_ni);
      }

      switch (buffer[len].type._type)
      {
      case LatticeTypes::V:
        lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
        break;
      case LatticeTypes::Cu:
        lats->cu_hash.emplace(latti_id);
        break;
      case LatticeTypes::Ni:
        lats->ni_hash.emplace(latti_id);
        break;
      case LatticeTypes::Mn:
        lats->mn_hash.emplace(latti_id);
      default:
        break;
      }
    }
  }else if(dimension == 1){
    for(len = 0; len < receive_len; len++){
      x = buffer[len].x;
      if(next_id == 2 || next_id == 3 || next_id == 6 || next_id == 7){
        buffer[len].y += sub_box_lattice_size[dimension];
      }else{
        buffer[len].y -= sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]);
      }
      y = buffer[len].y;
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      auto it_vac = lats->vac_hash.find(latti_id);
      auto it_cu = lats->cu_hash.find(latti_id);
      auto it_mn = lats->mn_hash.find(latti_id);
      auto it_ni = lats->ni_hash.find(latti_id);
      if (it_vac != lats->vac_hash.end()) {
        lats->vac_hash.erase(it_vac);
      }else if (it_cu != lats->cu_hash.end()){
        lats->cu_hash.erase(it_cu);
      }else if (it_mn != lats->mn_hash.end()){
        lats->mn_hash.erase(it_mn);
      }else if (it_ni != lats->ni_hash.end()){
        lats->ni_hash.erase(it_ni);
      }

      switch (buffer[len].type._type)
      {
      case LatticeTypes::V:
        lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
        break;
      case LatticeTypes::Cu:
        lats->cu_hash.emplace(latti_id);
        break;
      case LatticeTypes::Ni:
        lats->ni_hash.emplace(latti_id);
        break;
      case LatticeTypes::Mn:
        lats->mn_hash.emplace(latti_id);
      default:
        break;
      }
    }
  }else if(dimension == 0){
    for(len = 0; len < receive_len; len++){
      if(next_id % 2 == 1){
        buffer[len].x += BCC_DBX * (sub_box_lattice_size[dimension]);
      }else{
        buffer[len].x -= BCC_DBX * (sub_box_lattice_size[dimension] + (neighbour_local_sub_box[dimension] - sub_box_lattice_size[dimension]));
      }
      x = buffer[len].x;
      y = buffer[len].y;
      z = buffer[len].z;
      latti_id = lats->getId(x, y, z);
      auto it_vac = lats->vac_hash.find(latti_id);
      auto it_cu = lats->cu_hash.find(latti_id);
      auto it_mn = lats->mn_hash.find(latti_id);
      auto it_ni = lats->ni_hash.find(latti_id);
      if (it_vac != lats->vac_hash.end()) {
        lats->vac_hash.erase(it_vac);
      }else if (it_cu != lats->cu_hash.end()){
        lats->cu_hash.erase(it_cu);
      }else if (it_mn != lats->mn_hash.end()){
        lats->mn_hash.erase(it_mn);
      }else if (it_ni != lats->ni_hash.end()){
        lats->ni_hash.erase(it_ni);
      }

      switch (buffer[len].type._type)
      {
      case LatticeTypes::V:
        lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
        break;
      case LatticeTypes::Cu:
        lats->cu_hash.emplace(latti_id);
        break;
      case LatticeTypes::Ni:
        lats->ni_hash.emplace(latti_id);
        break;
      case LatticeTypes::Mn:
        lats->mn_hash.emplace(latti_id);
      default:
        break;
      }
    }
  }else{
    assert(false);
  }

  // 存储转发
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  for (int sect = 0; sect < 8; sect++) {
    for (int d = 1; d < comm::DIMENSION_SIZE; d++) {
      send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
      for (auto &r : send_regions[d]) {
        _type_lattice_size x_low = BCC_DBX * r.x_low;
        _type_lattice_size y_low = r.y_low;
        _type_lattice_size z_low = r.z_low;
        _type_lattice_size x_high = BCC_DBX * r.x_high;
        _type_lattice_size y_high = r.y_high;
        _type_lattice_size z_high = r.z_high;

        for(len = 0; len < receive_len; len++) {
          x = buffer[len].x;
          y = buffer[len].y;
          z = buffer[len].z;
          if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
            _type_lattice_id latti_id = lats->getId(x, y, z);
            if(d == 0){
              auto it_from = exchange_surface_x[sect].find(latti_id);
              if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(latti_id);
            }else if(d == 1) {
              auto it_from = exchange_surface_y[sect].find(latti_id);
              if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(latti_id);
            }else {
              auto it_from = exchange_surface_z[sect].find(latti_id);
              if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(latti_id);
            }
          }
        }
      }
    }
  }
 
}