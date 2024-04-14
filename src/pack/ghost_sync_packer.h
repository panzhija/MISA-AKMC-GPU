//
// Created by genshen on 2019/10/14.
//

#ifndef MISA_KMC_GHOST_SYNC_PACKER_H
#define MISA_KMC_GHOST_SYNC_PACKER_H

#include "algorithms/lattice_region_packer.h"
#include "lattice/lattice.h"
#include "lattice/lattices_list.h"
#include "utils/mpi_types.h"
#include <cassert>
#include <comm/domain/colored_domain.h>
#include <logs/logs.h>
#include <vector>
#include <string.h>
#include <comm/preset/sector_forwarding_region.h>
/**
 * \brief ghost sync means: before performing computing on the simulation area,
 *       the ghost area must be received for its neighbor processes.
 * \note: the region type(pack_region_type) is always comm::_type_lattice_coord.
 */
class GhostSyncPacker : public LatticeRegionPacker<Lattice> {
public:
  explicit GhostSyncPacker(LatticesList *lats_list);

  const unsigned long sendLength(const std::vector<comm::Region<pack_region_type>> send_regions, const int dimension,
                                 const int direction) override;

  const int sendLength2(const int dimension,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  void onSend(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> send_regions,
              const unsigned long send_len, const int dimension, const int direction) override;

  void onSend2(ChangeLattice buffer[], const int send_len, const int dimension,
               std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
               std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
               std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  void onReceive(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> recv_regions,
                 const unsigned long receive_len, const int dimension, const int direction) override;

  void onReceive2(ChangeLattice buffer[], const int receive_len, const int dimension,
                  std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                  std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                  std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                  const _type_lattice_count *sub_box_lattice_size, 
                  const _type_lattice_count *neighbour_local_sub_box, unsigned int next_id,
                  const comm::ColoredDomain *p_domain);

  static inline MPI_Datatype getMPI_DataType() { return mpi_types::_mpi_type_lattice_data; }

  static inline MPI_Datatype getMPI_DataTypeChange() { return mpi_types::_mpi_type_lattice_change; }

private:
  LatticesList *lats = nullptr;
};

#endif // MISA_KMC_GHOST_SYNC_PACKER_H
