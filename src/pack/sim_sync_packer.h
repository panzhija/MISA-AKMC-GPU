//
// Created by genshen on 2019/11/15.
//

#ifndef MISA_KMC_SIM_SYNC_PACKER_H
#define MISA_KMC_SIM_SYNC_PACKER_H

#include "algorithms/lattice_region_packer.h"
#include "utils/mpi_types.h"
#include <lattice/lattices_list.h>
#include <array>
#include <vector>
#include "lattice/lattice.h"
#include <comm/domain/colored_domain.h>
#include <cassert>
#include <comm/preset/sector_forwarding_region.h>
class SimSyncPacker : public LatticeRegionPacker<Lattice> {

public:
  explicit SimSyncPacker(LatticesList *lattice_list);

  const unsigned long sendLength(const std::vector<comm::Region<pack_region_type>> send_regions, const int dimension,
                                 const int direction) override;

  void onSend(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> send_regions,
              const unsigned long send_len, const int dimension, const int direction) override;

  void onReceive(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> recv_regions,
                 const unsigned long receive_len, const int dimension, const int direction) override;

  void onReceive2(ChangeLattice *buffer, std::vector<comm::Region<pack_region_type>> recv_regions,
                 const int receive_len, const int dimension, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                 const _type_lattice_count *sub_box_lattice_size, const _type_lattice_count *neighbour_local_sub_box, unsigned int id, 
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                 std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                 const comm::ColoredDomain *p_domain);

  const int sendLength2(const int dimension, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_x, ChangeLattice *send_count);

  void onSend2(ChangeLattice *buffer, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, const int dimension, 
               std::unordered_set<_type_lattice_id>& now_exchange_surface_x, ChangeLattice *send_count);

  static inline MPI_Datatype getMPI_DataType() { return mpi_types::_mpi_type_lattice_data; }
  
  static inline MPI_Datatype getMPI_DataTypeChange() { return mpi_types::_mpi_type_lattice_change; }

public:
  LatticesList *lats = nullptr;
  // std::array<std::vector<Lattice>, 7> &exchange_ghost_ref;
  // const comm::ColoredDomain *p_domain = nullptr;
  // int twice_count;
};

#endif // MISA_KMC_SIM_SYNC_PACKER_H
