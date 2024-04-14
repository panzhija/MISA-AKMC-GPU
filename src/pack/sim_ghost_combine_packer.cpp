#include "sim_ghost_combine_packer.h"
#include <utils/macros.h>

SimGhostCombinePacker::SimGhostCombinePacker(LatticesList *lattice_list) : lats(lattice_list){}
const unsigned long SimGhostCombinePacker::sendLength(const std::vector<comm::Region<pack_region_type>> send_regions,
                                              const int dimension, const int direction) {}

void SimGhostCombinePacker::onSend(Lattice *buffer, const std::vector<comm::Region<pack_region_type>> send_regions,
                           const unsigned long send_len, const int dimension, const int direction) {}

void SimGhostCombinePacker::onReceive(Lattice *buffer, const std::vector<comm::Region<pack_region_type>> recv_regions,
                              const unsigned long receive_len, const int dimension, const int direction) {}

int SimGhostCombinePacker::sendLengthGhost(const int dimension, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, ChangeLattice *send_count) {
  //int size_send = 0;

  return exchange_ghost[dimension % 7].size();

  // switch (dimension)
  // {
  // case 0:
  //   size_send += exchange_ghost[6].size();
  //   size_send += exchange_ghost[5].size();
  //   size_send += exchange_ghost[3].size();
  //   size_send += exchange_ghost[2].size();
  //   send_count->x = exchange_ghost[6].size();
  //   send_count->y = exchange_ghost[5].size();
  //   send_count->z = exchange_ghost[3].size();
  //   size_send++;
  //   break;
  // case 1:
  //   size_send += exchange_ghost[4].size();
  //   size_send += exchange_ghost[1].size();
  //   send_count->x = exchange_ghost[4].size();
  //   send_count->y = exchange_ghost[1].size();
  //   send_count->z = 0;
  //   size_send++;
  //   break;
  // case 2:
  //   size_send += exchange_ghost[0].size();
  //   break;
  // default:
  //   assert(false);
  //   break;
  // }
  // return size_send;
}

int SimGhostCombinePacker::sendLengthCombine(std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                                             std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {
  int size_send = 0;
  size_send += exchange_ghost[0].size();
  size_send += now_exchange_surface_z.size();
  return size_send;
}

int SimGhostCombinePacker::sendLengthSurface(const int dimension,
                                             std::unordered_set<_type_lattice_id>& now_exchange_surface_x,
                                             std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                                             std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {
  switch (dimension)
  {
    case 0:
      return now_exchange_surface_x.size();
    case 1:
      return now_exchange_surface_y.size();
    case 2:
      return now_exchange_surface_z.size();
    default:
      assert(false);
      break;
  }
}

void SimGhostCombinePacker::onSendGhost(const int dimension, ChangeLattice *buffer, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, ChangeLattice send_count) {
  int len = 0;

  for (const auto& lattice : exchange_ghost[dimension % 7]) {
    buffer[len++] = lattice;
  }
  // if(dimension == 0) {
  //   for (const auto& lattice : exchange_ghost[6]) {
  //     buffer[len++] = lattice;
  //   }
  //   for (const auto& lattice : exchange_ghost[5]) {
  //     buffer[len++] = lattice;
  //   }
  //   for (const auto& lattice : exchange_ghost[3]) {
  //     buffer[len++] = lattice;
  //   }
  //   for (const auto& lattice : exchange_ghost[2]) {
  //     buffer[len++] = lattice;
  //   }
  //   buffer[len++] = send_count;
  // }else if(dimension == 1) {
  //   for (const auto& lattice : exchange_ghost[4]) {
  //     buffer[len++] = lattice;
  //   }
  //   for (const auto& lattice : exchange_ghost[1]) {
  //     buffer[len++] = lattice;
  //   }
  //   buffer[len++] = send_count;
  // }else if(dimension == 2) {
  //   for (const auto& lattice : exchange_ghost[0]) {
  //     buffer[len++] = lattice;
  //   }
  // }
}

void SimGhostCombinePacker::onSendCombine(ChangeLattice *buffer, 
                                          std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                                          std::unordered_set<_type_lattice_id>& now_exchange_surface_z) {
  int len = 0;

  _type_lattice_size latti_x;
  _type_lattice_size latti_y;
  _type_lattice_size latti_z;
  for (const auto& lattice : exchange_ghost[0]) {
    buffer[len++] = lattice;
  }
  for (const auto& pair : now_exchange_surface_z) {
    latti_x = pair % lats->meta.size_x;
    latti_y = (pair / lats->meta.size_x) % lats->meta.size_y;
    latti_z = pair / (lats->meta.size_x * lats->meta.size_y);
    buffer[len].x = latti_x;
    buffer[len].y = latti_y;
    buffer[len].z = latti_z;
    auto it_vac = lats->vac_hash.find(pair);
    auto it_cu = lats->cu_hash.find(pair);
    auto it_mn = lats->mn_hash.find(pair);
    auto it_ni = lats->ni_hash.find(pair);
    if (it_vac != lats->vac_hash.end()) {
      buffer[len].type._type = LatticeTypes::V;
    }else if (it_cu != lats->cu_hash.end()){
      buffer[len].type._type = LatticeTypes::Cu;
    }else if (it_mn != lats->mn_hash.end()){
      buffer[len].type._type = LatticeTypes::Mn;
    }else if (it_ni != lats->ni_hash.end()){
      buffer[len].type._type = LatticeTypes::Ni;
    }else{
      buffer[len].type._type = LatticeTypes::Fe;
    }
    len++;
  }
}

void SimGhostCombinePacker::onSendSurface(const int dimension, ChangeLattice *buffer,
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
  }
  else if(dimension == 2) {
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
}

void SimGhostCombinePacker::onReceiveGhost(const int dimension, ChangeLattice *buffer, const int receive_len, 
                                           std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, unsigned int id, 
                                           std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                           std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                                           std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z, const comm::ColoredDomain *p_domain) {
  int len = 0;

  _type_lattice_id x, y, z;
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  _type_lattice_size x_low;
  _type_lattice_size y_low;
  _type_lattice_size z_low;
  _type_lattice_size x_high;
  _type_lattice_size y_high;
  _type_lattice_size z_high;
  _type_lattice_id latti_id;
  int count;
  for(len = 0; len < receive_len; len++){
      // if(id % 2 == 0){
      //   buffer[len].x += BCC_DBX * (p_domain->sub_box_lattice_size[dimension]);
      // }else{
      //   buffer[len].x -= BCC_DBX * (p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]));
      // }
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
    for (int sect = 0; sect < 8; sect++) {
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
          send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                         p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
        }
        for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
          for (auto &r : send_regions[d]) {
            x_low = BCC_DBX * r.x_low;
            y_low = r.y_low;
            z_low = r.z_low;
            x_high = BCC_DBX * r.x_high;
            y_high = r.y_high;
            z_high = r.z_high;
            if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
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

  //}
  // if(dimension == 0){
  //   count = buffer[receive_len - 1].x;
  //   for(len = 0; len < receive_len - 1; len++){
  //     if(id % 2 == 0){
  //       buffer[len].x += BCC_DBX * (p_domain->sub_box_lattice_size[dimension]);
  //     }else{
  //       buffer[len].x -= BCC_DBX * (p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]));
  //     }
  //     x = buffer[len].x;
  //     y = buffer[len].y;
  //     z = buffer[len].z;
  //     latti_id = lats->getId(x, y, z);
  //     auto it_vac = lats->vac_hash.find(latti_id);
  //     auto it_cu = lats->cu_hash.find(latti_id);
  //     auto it_mn = lats->mn_hash.find(latti_id);
  //     auto it_ni = lats->ni_hash.find(latti_id);
  //     if (it_vac != lats->vac_hash.end()) {
  //       lats->vac_hash.erase(it_vac);
  //     }else if (it_cu != lats->cu_hash.end()){
  //       lats->cu_hash.erase(it_cu);
  //     }else if (it_mn != lats->mn_hash.end()){
  //       lats->mn_hash.erase(it_mn);
  //     }else if (it_ni != lats->ni_hash.end()){
  //       lats->ni_hash.erase(it_ni);
  //     }

  //     switch (buffer[len].type._type)
  //     {
  //     case LatticeTypes::V:
  //       lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
  //       break;
  //     case LatticeTypes::Cu:
  //       lats->cu_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Ni:
  //       lats->ni_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Mn:
  //       lats->mn_hash.emplace(latti_id);
  //     default:
  //       break;
  //     }

  //     for (int sect = 0; sect < 8; sect++) {
  //       for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
  //                                                        p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
  //       }
  //       for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         for (auto &r : send_regions[d]) {
  //           x_low = BCC_DBX * r.x_low;
  //           y_low = r.y_low;
  //           z_low = r.z_low;
  //           x_high = BCC_DBX * r.x_high;
  //           y_high = r.y_high;
  //           z_high = r.z_high;
  //           if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
  //             if(d == 0){
  //               auto it_from = exchange_surface_x[sect].find(latti_id);
  //               if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(latti_id);
  //             }else if(d == 1) {
  //               auto it_from = exchange_surface_y[sect].find(latti_id);
  //               if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(latti_id);
  //             }else {
  //               auto it_from = exchange_surface_z[sect].find(latti_id);
  //               if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(latti_id);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }else if(dimension == 1){
  //   count = buffer[receive_len - 1].x;
  //   for(len = 0; len < receive_len - 1; len++){
  //     x = buffer[len].x;
  //     if(id == 0 || id == 1 || id == 4 || id == 5){
  //       buffer[len].y += p_domain->sub_box_lattice_size[dimension];
  //     }else{
  //       buffer[len].y -= p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]);
  //     }
  //     y = buffer[len].y;
  //     z = buffer[len].z;
  //     latti_id = lats->getId(x, y, z);
  //     auto it_vac = lats->vac_hash.find(latti_id);
  //     auto it_cu = lats->cu_hash.find(latti_id);
  //     auto it_mn = lats->mn_hash.find(latti_id);
  //     auto it_ni = lats->ni_hash.find(latti_id);
  //     if (it_vac != lats->vac_hash.end()) {
  //       lats->vac_hash.erase(it_vac);
  //     }else if (it_cu != lats->cu_hash.end()){
  //       lats->cu_hash.erase(it_cu);
  //     }else if (it_mn != lats->mn_hash.end()){
  //       lats->mn_hash.erase(it_mn);
  //     }else if (it_ni != lats->ni_hash.end()){
  //       lats->ni_hash.erase(it_ni);
  //     }

  //     switch (buffer[len].type._type)
  //     {
  //     case LatticeTypes::V:
  //       lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
  //       // lats->vac_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Cu:
  //       lats->cu_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Ni:
  //       lats->ni_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Mn:
  //       lats->mn_hash.emplace(latti_id);
  //     default:
  //       break;
  //     }

  //     for (int sect = 0; sect < 8; sect++) {
  //       for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
  //                                                        p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
  //       }
  //       for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         for (auto &r : send_regions[d]) {
  //           x_low = BCC_DBX * r.x_low;
  //           y_low = r.y_low;
  //           z_low = r.z_low;
  //           x_high = BCC_DBX * r.x_high;
  //           y_high = r.y_high;
  //           z_high = r.z_high;
  //           if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
  //             if(d == 0){
  //               auto it_from = exchange_surface_x[sect].find(latti_id);
  //               if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(latti_id);
  //             }else if(d == 1) {
  //               auto it_from = exchange_surface_y[sect].find(latti_id);
  //               if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(latti_id);
  //             }else {
  //               auto it_from = exchange_surface_z[sect].find(latti_id);
  //               if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(latti_id);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }else if(dimension == 2){
  //   for(len = 0; len < receive_len; len++){
  //     x = buffer[len].x;
  //     y = buffer[len].y;
  //     if(id < 4){
  //       buffer[len].z += p_domain->sub_box_lattice_size[dimension];
  //     }else{
  //       buffer[len].z -= p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]);
  //     }
  //     z = buffer[len].z;
  //     latti_id = lats->getId(x, y, z);
  //     auto it_vac = lats->vac_hash.find(latti_id);
  //     auto it_cu = lats->cu_hash.find(latti_id);
  //     auto it_mn = lats->mn_hash.find(latti_id);
  //     auto it_ni = lats->ni_hash.find(latti_id);
  //     if (it_vac != lats->vac_hash.end()) {
  //       lats->vac_hash.erase(it_vac);
  //     }else if (it_cu != lats->cu_hash.end()){
  //       lats->cu_hash.erase(it_cu);
  //     }else if (it_mn != lats->mn_hash.end()){
  //       lats->mn_hash.erase(it_mn);
  //     }else if (it_ni != lats->ni_hash.end()){
  //       lats->ni_hash.erase(it_ni);
  //     }

  //     switch (buffer[len].type._type)
  //     {
  //     case LatticeTypes::V:
  //       lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
  //       break;
  //     case LatticeTypes::Cu:
  //       lats->cu_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Ni:
  //       lats->ni_hash.emplace(latti_id);
  //       break;
  //     case LatticeTypes::Mn:
  //       lats->mn_hash.emplace(latti_id);
  //     default:
  //       break;
  //     }

  //     for (int sect = 0; sect < 8; sect++) {
  //       for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
  //                                                        p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
  //       }
  //       for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
  //         for (auto &r : send_regions[d]) {
  //           x_low = BCC_DBX * r.x_low;
  //           y_low = r.y_low;
  //           z_low = r.z_low;
  //           x_high = BCC_DBX * r.x_high;
  //           y_high = r.y_high;
  //           z_high = r.z_high;
  //           if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
  //             if(d == 0){
  //               auto it_from = exchange_surface_x[sect].find(latti_id);
  //               if(it_from == exchange_surface_x[sect].end()) exchange_surface_x[sect].emplace(latti_id);
  //             }else if(d == 1) {
  //               auto it_from = exchange_surface_y[sect].find(latti_id);
  //               if(it_from == exchange_surface_y[sect].end()) exchange_surface_y[sect].emplace(latti_id);
  //             }else {
  //               auto it_from = exchange_surface_z[sect].find(latti_id);
  //               if(it_from == exchange_surface_z[sect].end()) exchange_surface_z[sect].emplace(latti_id);
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // 存储转发
  // if(dimension == 0){
  //   int last_count = receive_len - 1 - buffer[receive_len - 1].x - buffer[receive_len - 1].y - buffer[receive_len - 1].z;
  //   if(buffer[receive_len - 1].y > 0) {
  //     std::copy(buffer + count, buffer + count + buffer[receive_len - 1].y, std::back_inserter(exchange_ghost[4]));
  //     count += buffer[receive_len - 1].y;
  //   }
  //   if(buffer[receive_len - 1].z > 0) {
  //     std::copy(buffer + count, buffer + count + buffer[receive_len - 1].z, std::back_inserter(exchange_ghost[1]));
  //     count += buffer[receive_len - 1].z;
  //   }
  //   if(last_count > 0) {
  //     std::copy(buffer + count, buffer + count + last_count, std::back_inserter(exchange_ghost[0]));
  //   }
  // }else if(dimension == 1){
  //   if(buffer[receive_len - 1].y > 0){
  //     std::copy(buffer + count, buffer + count + buffer[receive_len - 1].y, std::back_inserter(exchange_ghost[0]));
  //   }
  // }
}

void SimGhostCombinePacker::onReceiveCombine(const int dimension, ChangeLattice *buffer, const int receive_len, unsigned int id,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                                             const comm::ColoredDomain *p_domain){
  int len = 0;

  _type_lattice_id x, y, z;
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  _type_lattice_size x_low;
  _type_lattice_size y_low;
  _type_lattice_size z_low;
  _type_lattice_size x_high;
  _type_lattice_size y_high;
  _type_lattice_size z_high;
  _type_lattice_id latti_id;
  for(len = 0; len < receive_len; len++){
    x = buffer[len].x;
    y = buffer[len].y;
    if(id < 4){
      buffer[len].z += p_domain->sub_box_lattice_size[dimension];
    }else{
      buffer[len].z -= p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]);
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
    for (int sect = 0; sect < 8; sect++) {
      for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
        send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                       p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
      }
      for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
        for (auto &r : send_regions[d]) {
          x_low = BCC_DBX * r.x_low;
          y_low = r.y_low;
          z_low = r.z_low;
          x_high = BCC_DBX * r.x_high;
          y_high = r.y_high;
          z_high = r.z_high;
          if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
            latti_id = lats->getId(x, y, z);
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

void SimGhostCombinePacker::onReceiveSurface(const int dimension, ChangeLattice *buffer,
                                             const int receive_len, unsigned int next_id,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                                             std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                                             const comm::ColoredDomain *p_domain) {
  int len = 0;

  _type_lattice_id x, y, z;
  const int dims[comm::DIMENSION_SIZE] = {comm::DIM_X, comm::DIM_Y, comm::DIM_Z};
  std::array<std::vector<comm::Region<comm::_type_lattice_coord>>, comm::DIMENSION_SIZE> send_regions; // send regions in each dimension
  _type_lattice_size x_low;
  _type_lattice_size y_low;
  _type_lattice_size z_low;
  _type_lattice_size x_high;
  _type_lattice_size y_high;
  _type_lattice_size z_high;
  _type_lattice_id latti_id;
  if(dimension == 0){
    for(len = 0; len < receive_len; len++){
      if(next_id % 2 == 1){
        buffer[len].x += BCC_DBX * (p_domain->sub_box_lattice_size[dimension]);
      }else{
        buffer[len].x -= BCC_DBX * (p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]));
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

      for (int sect = 0; sect < 8; sect++) {
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
          send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                     p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
        }
        for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
          for (auto &r : send_regions[d]) {
            x_low = BCC_DBX * r.x_low;
            y_low = r.y_low;
            z_low = r.z_low;
            x_high = BCC_DBX * r.x_high;
            y_high = r.y_high;
            z_high = r.z_high;
            if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
              latti_id = lats->getId(x, y, z);
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
  }else if(dimension == 1){
    for(len = 0; len < receive_len; len++){
      x = buffer[len].x;
      if(next_id == 2 || next_id == 3 || next_id == 6 || next_id == 7){
        buffer[len].y += p_domain->sub_box_lattice_size[dimension];
      }else{
        buffer[len].y -= p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]);
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

      for (int sect = 0; sect < 8; sect++) {
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
          send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                         p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
        }
        for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
          for (auto &r : send_regions[d]) {
            x_low = BCC_DBX * r.x_low;
            y_low = r.y_low;
            z_low = r.z_low;
            x_high = BCC_DBX * r.x_high;
            y_high = r.y_high;
            z_high = r.z_high;
            if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
              latti_id = lats->getId(x, y, z);
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
  }else if(dimension == 2){
    for(len = 0; len < receive_len; len++){
      x = buffer[len].x;
      y = buffer[len].y;
      if(next_id > 3){
        buffer[len].z += p_domain->sub_box_lattice_size[dimension];
      }else{
        buffer[len].z -= p_domain->sub_box_lattice_size[dimension] + (p_domain->neighbour_local_sub_box[dimension] - p_domain->sub_box_lattice_size[dimension]);
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

      for (int sect = 0; sect < 8; sect++) {
        for (int d = 0; d < comm::DIMENSION_SIZE; d++) {
          send_regions[d] = comm::fwCommSectorSendRegion(sect, dims[d], p_domain->lattice_size_ghost,
                                                         p_domain->local_split_coord, p_domain->local_sub_box_lattice_region);
        }
        for(int d = 0; d < comm::DIMENSION_SIZE; d++) {
          for (auto &r : send_regions[d]) {
            x_low = BCC_DBX * r.x_low;
            y_low = r.y_low;
            z_low = r.z_low;
            x_high = BCC_DBX * r.x_high;
            y_high = r.y_high;
            z_high = r.z_high;
            if(x_low <= x && x < x_high && y_low <= y && y < y_high && z_low <= z && z < z_high){
              latti_id = lats->getId(x, y, z);
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
}