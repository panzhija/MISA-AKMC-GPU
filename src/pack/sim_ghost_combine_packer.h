#ifndef MISA_KMC_SIM_GHOSY_COMBINE_PACKER_H
#define MISA_KMC_SIM_GHOSY_COMBINE_PACKER_H

#include "algorithms/lattice_region_packer.h"
#include "utils/mpi_types.h"
#include <lattice/lattices_list.h>
#include <array>
#include <vector>
#include "lattice/lattice.h"
#include <comm/domain/colored_domain.h>
#include <cassert>
#include <comm/preset/sector_forwarding_region.h>
class SimGhostCombinePacker : public LatticeRegionPacker<Lattice> {

public:
  explicit SimGhostCombinePacker(LatticesList *lattice_list);

  const unsigned long sendLength(const std::vector<comm::Region<pack_region_type>> send_regions, const int dimension,
                                 const int direction) override;

  void onSend(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> send_regions,
              const unsigned long send_len, const int dimension, const int direction) override;

  void onReceive(pack_date_type buffer[], const std::vector<comm::Region<pack_region_type>> recv_regions,
                 const unsigned long receive_len, const int dimension, const int direction) override;
  
  int sendLengthGhost(const int dimension, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, ChangeLattice *send_count);

  int sendLengthCombine(std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  int sendLengthSurface(const int dimension,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_x,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                        std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  void onSendGhost(const int dimension, ChangeLattice *buffer, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, ChangeLattice send_count);

  void onSendCombine(ChangeLattice *buffer, std::array<std::vector<ChangeLattice>, 7>& exchange_ghost,
                     std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  void onSendSurface(const int dimension, ChangeLattice *buffer,
                     std::unordered_set<_type_lattice_id>& now_exchange_surface_x, 
                     std::unordered_set<_type_lattice_id>& now_exchange_surface_y,
                     std::unordered_set<_type_lattice_id>& now_exchange_surface_z);

  void onReceiveGhost(const int dimension, ChangeLattice *buffer, const int receive_len, 
                      std::array<std::vector<ChangeLattice>, 7>& exchange_ghost, unsigned int id,
                      std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                      std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                      std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z, const comm::ColoredDomain *p_domain);
  
  void onReceiveCombine(const int dimension, ChangeLattice *buffer, const int receive_len, unsigned int id,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                        const comm::ColoredDomain *p_domain);

  void onReceiveSurface(const int dimension, ChangeLattice *buffer,
                        const int receive_len, unsigned int next_id,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_x,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_y,
                        std::array<std::unordered_set<_type_lattice_id>, 8>& exchange_surface_z,
                        const comm::ColoredDomain *p_domain);
  static inline MPI_Datatype getMPI_DataTypeChange() { return mpi_types::_mpi_type_lattice_change; }

public:
  LatticesList *lats = nullptr;
};

#endif // MISA_KMC_SIM_GHOSY_COMBINE_PACKER_H