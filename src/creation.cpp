//
// Created by genshen on 2019/10/6.
//

#include "creation.h"
#include "type_define.h"
#include "utils/macros.h"
#include "utils/random/random.h"
#include "utils/random/rng_type.h"
#include "utils/simulation_domain.h"
#include <cassert>
#include <fstream>
#include <utils/bundle.h>
#include <utils/mpi_utils.h>
#include <iostream>
#include <logs/logs.h>

void creation::createAtomsRandom(int64_t seed_create_types, LatticesList *lats,
                                 const std::vector<LatticeTypes::lat_type> types,
                                 const std::vector<int64_t> types_ratio, const comm::ColoredDomain *p_domain) {
  int64_t ratio_total = 0;
  for (int64_t i = 0; i < types.size(); i++) {
    ratio_total += types_ratio[i];
  }
  std::uniform_int_distribution<> types_dis(1, ratio_total);
  r::type_rng types_rng(seed_create_types);
  // todo skip ghost lattices.
  int64_t n = 0;

  for (_type_lattice_size z = 0; z < lats->meta.size_z; z++) {
    for (_type_lattice_size y = 0; y < lats->meta.size_y; y++) {
      for (_type_lattice_size x = 0; x < lats->meta.size_x; x++) {
        assert(n == x + lats->meta.size_x * (y + z * lats->meta.size_y));
        const int64_t rand_hit = types_dis(types_rng);
        LatticeTypes::lat_type latti_type = LatticeTypes::randomAtomsType(types.data(), types_ratio.data(), types.size(), rand_hit);
        // lats->_lattices[z][y][x].type._type = latti_type;
        switch (latti_type)
        {
        case LatticeTypes::Cu:
          lats->cu_hash.emplace(n);
          break;
        case LatticeTypes::Ni:
          lats->ni_hash.emplace(n);
          break;
        case LatticeTypes::Mn:
          lats->mn_hash.emplace(n);
          break;
        case LatticeTypes::Fe:
          // 不存储Fe
          break;
        default:
          break;
        }
        n++;
      }
    }
  }
  // unsigned long cu_count = 0;
  // unsigned long ni_count = 0;
  // unsigned long mn_count = 0;
  // kiwi::logs::v(" ", " Number of elements in cu_hash is : {} .\n", lats->cu_hash.size());
  // kiwi::logs::v(" ", " Number of elements in ni_hash is : {} .\n", lats->ni_hash.size());
  // kiwi::logs::v(" ", " Number of elements in mn_hash is : {} .\n", lats->mn_hash.size());
  // for (_type_lattice_size z = 0; z < lats->meta.size_z; z++) {
  //   for (_type_lattice_size y = 0; y < lats->meta.size_y; y++) {
  //     for (_type_lattice_size x = 0; x < lats->meta.size_x; x++) {
  //       switch (lats->_lattices[z][y][x].type._type)
  //       {
  //       case LatticeTypes::Cu:
  //         cu_count++;
  //         break;
  //       case LatticeTypes::Ni:
  //         ni_count++;
  //         break;
  //       case LatticeTypes::Mn:
  //         mn_count++;
  //         break;
  //       default:
  //         break;
  //       }
  //     }
  //   }
  // }
  // kiwi::logs::v(" ", " Number of cu elements in latti is : {} .\n", cu_count);
  // kiwi::logs::v(" ", " Number of ni elements in latti is : {} .\n", ni_count);
  // kiwi::logs::v(" ", " Number of mn elements in latti is : {} .\n", mn_count);
  // lats->forAllLattices(
  //     [&](const _type_lattice_coord x, const _type_lattice_coord y, const _type_lattice_coord z, Lattice &lattice) {
  //       // set lattice types randomly.
  //       const unsigned int rand_hit = types_dis(types_rng);
  //       lattice.type._type = LatticeTypes::randomAtomsType(types.data(), types_ratio.data(), types.size(), rand_hit);
  //       return true;
  //     });
}

void creation::createRandom(int64_t seed_create_types, int64_t seed_create_vacancy, LatticesList *lats,
                            VacancyList *va_list, const std::vector<LatticeTypes::lat_type> types,
                            const std::vector<int64_t> types_ratio, const int64_t va_count,
                            const comm::ColoredDomain *p_domain) {
  // create lattice types
  createAtomsRandom(seed_create_types, lats, types, types_ratio, p_domain);

  // create vacancy
  // to use SimulationDomain, make sure this function is called after
  // simulation::createDomain()
  // const comm::_type_lattice_size max_lattice_size = (lats->meta.box_x - 2 * lats->meta.ghost_x) * (lats->meta.box_y - 2 * lats->meta.ghost_y) * (lats->meta.box_z - 2 * lats->meta.ghost_z);
  // const int64_t cu_local =
  //     types_ratio[1] / SimulationDomain::comm_sim_pro.all_ranks +
  //     (SimulationDomain::comm_sim_pro.own_rank < (types_ratio[1] % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  // std::uniform_int_distribution<long long> types_dis(0, max_lattice_size - 1);
  // r::type_rng types_rng(seed_create_types);
  // for(_type_lattice_count cu_i = 0; cu_i < cu_local;) {
  //   _type_lattice_id local_id = types_dis(types_rng);
  //   _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
  //   local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
  //   _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
  //   _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
  //   _type_lattice_id latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);

  //   auto it_cu = lats->cu_hash.find(latti_id);
  //   if(it_cu == lats->cu_hash.end()){
  //     lats->cu_hash.emplace(latti_id);
  //     cu_i++;
  //   }
  // }
  const comm::_type_lattice_size max_lattice_size = lats->meta.box_x * lats->meta.box_y * lats->meta.box_z;
  const int64_t va_local =
      va_count / SimulationDomain::comm_sim_pro.all_ranks +
      (SimulationDomain::comm_sim_pro.own_rank < (va_count % SimulationDomain::comm_sim_pro.all_ranks) ? 1 : 0);
  std::uniform_int_distribution<long long> va_dis(0, max_lattice_size - 1);
  r::type_rng va_rng(seed_create_vacancy);
#ifdef KMC_DEBUG_MODE
  assert(max_lattice_size >= 1);
#endif
  for (int64_t va_i = 0; va_i < va_local;) {
    int64_t local_id = va_dis(va_rng);
    // get the coordinate of the lattice in simulation box
    // _type_lattice_coord x = local_id % (lats->meta.box_x - 2 * lats->meta.ghost_x);
    // local_id = local_id / (lats->meta.box_x - 2 * lats->meta.ghost_x);
    // _type_lattice_coord y = local_id % (lats->meta.box_y - 2 * lats->meta.ghost_y);
    // _type_lattice_coord z = local_id / (lats->meta.box_y - 2 * lats->meta.ghost_y);
    // int64_t latti_id = lats->getId(x + 2 * lats->meta.ghost_x, y + 2 * lats->meta.ghost_y, z + 2 * lats->meta.ghost_z);
    _type_lattice_coord x = local_id % lats->meta.box_x;
    local_id = local_id / lats->meta.box_x;
    _type_lattice_coord y = local_id % lats->meta.box_y;
    _type_lattice_coord z = local_id / lats->meta.box_y;
    unsigned long latti_id = lats->getId(x + lats->meta.ghost_x, y + lats->meta.ghost_y, z + lats->meta.ghost_z);
    //Lattice &lat = lats->getLat(x + lats->meta.ghost_x, y + lats->meta.ghost_y, z + lats->meta.ghost_z);
    auto it = lats->vac_hash.find(latti_id);
    if (it == lats->vac_hash.end()) {
      auto it_cu = lats->cu_hash.find(latti_id);
      if (it_cu != lats->cu_hash.end()){
        lats->cu_hash.erase(it_cu);
      }else{
        // Fe不存
      }
      lats->vac_hash.emplace(std::make_pair(latti_id, VacancyHash{}));
      va_i++;
    }
  }
 // if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " Vac count is : {} Cu count is {} .\n", va_local, cu_local);
  // todo init interval list
  // va_list->reindex(lats, p_domain->local_sub_box_lattice_region);
}

void creation::setGlobalId(LatticesList *lats_list, const comm::Region<comm::_type_lattice_coord> lbr,
                           const comm::Region<comm::_type_lattice_coord> gbr,
                           std::array<int64_t, comm::DIMENSION_SIZE> phase_space) {
  for (comm::_type_lattice_coord z = lbr.z_low; z < lbr.z_high; z++) {
    for (comm::_type_lattice_coord y = lbr.y_low; y < lbr.y_high; y++) {
      for (comm::_type_lattice_coord x = BCC_DBX * lbr.x_low; x < BCC_DBX * lbr.x_high; x++) {
        Lattice &lattice = lats_list->getLat(x, y, z);
        const comm::_type_lattice_coord gx = (x - lats_list->meta.ghost_x) + BCC_DBX * gbr.x_low;
        const comm::_type_lattice_coord gy = (y - lats_list->meta.ghost_y) + gbr.y_low;
        const comm::_type_lattice_coord gz = (z - lats_list->meta.ghost_z) + gbr.z_low;
        // convert xyz to global id.
        const comm::_type_lattice_coord gid = gx + BCC_DBX * phase_space[0] * (gy + gz * phase_space[1]);
        lattice.id = gid;
      }
    }
  }
}

void creation::createFromPile(const std::string pipe_file, int64_t seed_create_types, LatticesList *lats,
                              VacancyList *va_list, const std::vector<LatticeTypes::lat_type> types,
                              const std::vector<int64_t> types_ratio, const comm::ColoredDomain *p_domain) {
  createAtomsRandom(seed_create_types, lats, types, types_ratio, p_domain);

  // open vacancies file
  std::vector<std::array<_type_lattice_coord, 3>> vac_positions{};
  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
    std::ifstream vac_file;
    vac_file.open(pipe_file);
    _type_lattice_coord x, y, z;
    while (vac_file >> x >> y >> z) {
      vac_positions.push_back({x, y, z});
    }
  }

  kiwi::Bundle bundle = kiwi::Bundle();
  bundle.newPackBuffer(1024 * 1024); // make sure it is large enough

  if (kiwi::mpiUtils::global_process.own_rank == MASTER_PROCESSOR) {
    bundle.put(vac_positions);
  }
  MPI_Bcast(bundle.getPackedData(), bundle.getPackedDataCap(), MPI_BYTE, MASTER_PROCESSOR,
            MPI_COMM_WORLD); // synchronize data

  if (kiwi::mpiUtils::global_process.own_rank != MASTER_PROCESSOR) { // unpack
    int cursor = 0;
    bundle.get(cursor, vac_positions);
  }
  bundle.freePackBuffer();

  // place vacancies
  for (auto vac : vac_positions) {
    _type_lattice_coord x = vac[0]; // vac[0] is doubled.
    _type_lattice_coord y = vac[1];
    _type_lattice_coord z = vac[2];
    if (p_domain->local_sub_box_lattice_region.isIn(x / 2, y, z)) {
      Lattice &lat = lats->getLat(x + lats->meta.ghost_x, y + lats->meta.ghost_y, z + lats->meta.ghost_z);
      if (!lat.type.isVacancy()) {
        lat.type._type = LatticeTypes::V;
      }
    }
  }

  // update vacancies list
  va_list->reindex(lats, p_domain->local_sub_box_lattice_region);
}