//
// Created by runchu on 2019-10-16.
//

#include "ghost_init_packer.h"
#include "utils/macros.h"
#include <comm/preset/comm_forwarding_region.h>
#include <iostream>
#include <omp.h>

GhostInitPacker::GhostInitPacker(const comm::ColoredDomain *p_domain, LatticesList *lats_list)
    : p_domain(p_domain), lats(lats_list) {}

const unsigned long GhostInitPacker::sendLength(const int dimension, const int direction) {
  comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  return BCC_DBX * send_region.volume();
}

int GhostInitPacker::sendLength2(const int dimension, const int direction) {
  comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;

  for (_type_lattice_id z = send_region.z_low; z < send_region.z_high; z++) {
    for (_type_lattice_id y = send_region.y_low; y < send_region.y_high; y++) {
      for (_type_lattice_id x = send_region.x_low; x < send_region.x_high; x++) {
        unsigned long latti_id = lats->getId(BCC_DBX * x, y, z);
        // lats->isUniqueLattiId(latti_id);
        auto it_vac = lats->vac_hash.find(latti_id);
        auto it_cu = lats->cu_hash.find(latti_id);
        auto it_mn = lats->mn_hash.find(latti_id);
        auto it_ni = lats->ni_hash.find(latti_id);
        if (!(it_vac == lats->vac_hash.end() && it_cu == lats->cu_hash.end() && it_mn == lats->mn_hash.end() && it_ni == lats->ni_hash.end())) {
          len++;
        }
        
        unsigned long latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
        // lats->isUniqueLattiId(latti_id2);
        auto it_vac2 = lats->vac_hash.find(latti_id2);
        auto it_cu2 = lats->cu_hash.find(latti_id2);
        auto it_mn2 = lats->mn_hash.find(latti_id2);
        auto it_ni2 = lats->ni_hash.find(latti_id2);
        if (!(it_vac2 == lats->vac_hash.end() && it_cu2 == lats->cu_hash.end() && it_mn2 == lats->mn_hash.end() && it_ni2 == lats->ni_hash.end())) {
          len++;
        }
      }
    }
  }

  return len;

}
void GhostInitPacker::onSend(Lattice *buffer, const unsigned long send_len, const int dimension,
              const int direction){
  const comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  unsigned long len = 0;
  for (int z = send_region.z_low; z < send_region.z_high; z++) {
    for (int y = send_region.y_low; y < send_region.y_high; y++) {
      for (int x = send_region.x_low; x < send_region.x_high; x++) {
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

void GhostInitPacker::onSend2(Lattice *buffer, int send_len, const int dimension,
                             const int direction) {
  const comm::Region<comm::_type_lattice_size> send_region = comm::fwCommLocalSendRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;
  _type_lattice_count index = 0;

  for (_type_lattice_id z = send_region.z_low; z < send_region.z_high; z++) {
    for (_type_lattice_id y = send_region.y_low; y < send_region.y_high; y++) {
      for (_type_lattice_id x = send_region.x_low; x < send_region.x_high; x++) {
        _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
        
        auto it_vac = lats->vac_hash.find(latti_id);
        auto it_cu = lats->cu_hash.find(latti_id);
        auto it_mn = lats->mn_hash.find(latti_id);
        auto it_ni = lats->ni_hash.find(latti_id);
        // lats->isUniqueLattiId(latti_id);
        // 这里的 id 存的是 本次传输的数组下标，而不是原子的局部 id
        if (it_vac != lats->vac_hash.end()) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::V};
          index++;
        }else if (it_cu != lats->cu_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Cu};
          index++;
        }else if (it_mn != lats->mn_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Mn};
          index++;
        }else if (it_ni != lats->ni_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Ni};
          index++;
        }else{
          // Fe
        }
        len++;
        
        _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
        auto it_vac2 = lats->vac_hash.find(latti_id2);
        auto it_cu2 = lats->cu_hash.find(latti_id2);
        auto it_mn2 = lats->mn_hash.find(latti_id2);
        auto it_ni2 = lats->ni_hash.find(latti_id2);
        // lats->isUniqueLattiId(latti_id2);
        if (it_vac2 != lats->vac_hash.end()) {
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::V};
          index++;
        }else if (it_cu2 != lats->cu_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Cu};
          index++;
        }else if (it_mn2 != lats->mn_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Mn};
          index++;
        }else if (it_ni2 != lats->ni_hash.end()){
          buffer[index].id = len;
          buffer[index].type = LatticeTypes{LatticeTypes::Ni};
          index++;
        }else{
          // Fe
        }
        len++;
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x];
        // buffer[len++] = lats->_lattices[z][y][BCC_DBX * x + 1];
      }
    }
  }
  assert(index == send_len);
}

void GhostInitPacker::onReceive(Lattice buffer[], const unsigned long receive_len, const int dimension,
                 const int direction){
  const comm::Region<comm::_type_lattice_size> recv_region = comm::fwCommLocalRecvRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);
  unsigned long len = 0;
  for (int z = recv_region.z_low; z < recv_region.z_high; z++) {
    for (int y = recv_region.y_low; y < recv_region.y_high; y++) {
      for (int x = recv_region.x_low; x < recv_region.x_high; x++) {
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
        case LatticeTypes::Fe:
          // 不做操作
          break;
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
        case LatticeTypes::Fe:
          // 不做操作
          break;
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

void GhostInitPacker::onReceive2(Lattice *buffer, const int receive_len, const int dimension, const int direction) {
  const comm::Region<comm::_type_lattice_size> recv_region = comm::fwCommLocalRecvRegion(
      p_domain->lattice_size_ghost, p_domain->local_sub_box_lattice_region, dimension, direction);

  int len = 0;
  _type_lattice_count index = 0;

  if(receive_len == 0){
    // 说明全是Fe原子
    for (_type_lattice_id z = recv_region.z_low; z < recv_region.z_high; z++) {
      for (_type_lattice_id y = recv_region.y_low; y < recv_region.y_high; y++) {
        for (_type_lattice_id x = recv_region.x_low; x < recv_region.x_high; x++) {
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          // lats->isUniqueLattiId(latti_id);
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
          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          // lats->isUniqueLattiId(latti_id2);
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
        }
      }
    }
  }else{
    
    // kiwi::logs::v(" ", " DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD .\n");
    for (_type_lattice_id z = recv_region.z_low; z < recv_region.z_high; z++) {
      for (_type_lattice_id y = recv_region.y_low; y < recv_region.y_high; y++) {
        for (_type_lattice_id x = recv_region.x_low; x < recv_region.x_high; x++) {
          
          _type_lattice_id latti_id = lats->getId(BCC_DBX * x, y, z);
          // lats->isUniqueLattiId(latti_id);
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

          if(buffer[index].id == len){
            // 在 buffer 中找到了
            switch (buffer[index].type._type)
            {
            case LatticeTypes::V:
              lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
              break;
            case LatticeTypes::Cu:
              lats->cu_hash.emplace(latti_id);
              break;
            case LatticeTypes::Mn:
              lats->mn_hash.emplace(latti_id);
              break;
            case LatticeTypes::Ni:
              lats->ni_hash.emplace(latti_id);
              break;
            default:
              break;
            }
            index++;
          }
          len++;

          _type_lattice_id latti_id2 = lats->getId(BCC_DBX * x + 1, y, z);
          // lats->isUniqueLattiId(latti_id);
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

          if(buffer[index].id == len){
            switch (buffer[index].type._type)
            {
            case LatticeTypes::V:
              lats->vac_hash.emplace(std::make_pair(latti_id2, VacancyHash{}));
              break;
            case LatticeTypes::Cu:
              lats->cu_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Mn:
              lats->mn_hash.emplace(latti_id2);
              break;
            case LatticeTypes::Ni:
              lats->ni_hash.emplace(latti_id2);
              break;
            default:
              break;
            }
            index++;
          }
          len++;
        }
      }
    }
  }
  assert(index == receive_len);
  // std::cout << "3333333333333333333333333333333333" << std::endl;
  // todo set more information if the lattice is vacancy or dumbbell.
}
