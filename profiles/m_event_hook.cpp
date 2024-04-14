//
// Created by genshen on 2019/12/20.
//

#include "m_event_hook.h"
#include "io/lattice_dump.h"
#include "utils/simulation_domain.h"
#include <logs/logs.h>
#include <string>
#include "../models/abvi/kmc.h"
#include <mpi.h>
#include "utils/simulation_domain.h"

MEventHook::MEventHook(const conf::Output &output_conf, LatticesList *p_lattice_list, counter *p_counter)
    : output_config(output_conf), p_lattice_list(p_lattice_list), p_counter(p_counter) {}

void MEventHook::onStepFinished(unsigned long step) {
  // 输出孤立 Cu 原子的个数
  if (step % 100 == 0){
    comm::Region<_type_lattice_coord> region(p_lattice_list->meta.ghost_x, p_lattice_list->meta.ghost_y, p_lattice_list->meta.ghost_z,
                                             p_lattice_list->meta.ghost_x + p_lattice_list->meta.box_x,
                                             p_lattice_list->meta.ghost_y + p_lattice_list->meta.box_y,
                                             p_lattice_list->meta.ghost_z + p_lattice_list->meta.box_z);

    Lattice _1nn_list[LatticesList::MAX_1NN];
    int single_cu = 0;
    int total = 0;
    for (_type_lattice_coord z = region.z_low; z < region.z_high; z++) {
      for (_type_lattice_coord y = region.y_low; y < region.y_high; y++) {
        for (_type_lattice_coord x = region.x_low; x < region.x_high; x++) {
          unsigned long latti_id = p_lattice_list->getId(x, y, z);
          auto it_cu = p_lattice_list->cu_hash.find(latti_id);
          if (it_cu != p_lattice_list->cu_hash.end()) {
            p_lattice_list->get1nn2(x, y, z, _1nn_list);
            const _type_neighbour_status nei_status = p_lattice_list->get1nnBoundaryStatus(x, y, z);
            bool has_neighbor_cu = false;
            // travel its 1nn  neighbor
            for (int b = 0; b < LatticesList::MAX_NEI_BITS; b++) {
              // the neighbour lattice is available.
              // and the neighbour lattice is also Cu
              if (((nei_status >> b) & 1) && (_1nn_list[b].type._type == LatticeTypes::Cu)) {
                has_neighbor_cu = true;
                break;
              }
            }
            if (!has_neighbor_cu) {
              single_cu++;
            }
          }
        }
      }
    }
    MPI_Reduce(&single_cu, &total, 1, MPI_INT, MPI_SUM, 0, SimulationDomain::comm_sim_pro.comm);
    if(SimulationDomain::comm_sim_pro.own_rank == 0) kiwi::logs::v(" ", " step {} finished. now single_cu is {}.\n", step, total);
  }

  // if (output_config.logs_interval != 0) {
  //   if (step % 200 == 0 && SimulationDomain::comm_sim_pro.own_rank == 0) {
  //     kiwi::logs::d("hook", "finish step: {}\n", step);
  //   }
  // }

  if (output_config.dump_interval != 0) {
    if (step % output_config.dump_interval == 0) {
    // if (step == 199) {
    // if (step == 0 || step == 49999 || step == 99999 || step == 149999 || step == 199999 || step == 249999 || step == 299999 || step == 349999 || step == 399999 || step == 449999 || step == 499999 ||step == 999999 || step == 2499999 || step == 4999999 || step == 9999999) {
      std::string dump_file_path = fmt::format(output_config.dump_file_path, step);
      const kiwi::mpi_process process = SimulationDomain::kiwi_sim_pro;
      if (process.all_ranks != 1) {
        dump_file_path = fmt::format(output_config.dump_file_path + ".{}", step, process.own_rank);
      }
      kiwi::logs::v("dump", "dumping to file {} at step {}.\n", dump_file_path, step);
      LatticeDump dump;
      dump.dump(dump_file_path, p_lattice_list, step);
    }
  }

}

void MEventHook::onAllDone() {}
