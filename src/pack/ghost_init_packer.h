//
// Created by runchu on 2019-10-16.
//

#ifndef MISA_KMC_GHOST_INIT_PACKER_H
#define MISA_KMC_GHOST_INIT_PACKER_H

#include "lattice/lattices_list.h"
#include <comm/comm.hpp>
#include <comm/domain/colored_domain.h>
#include "utils/mpi_types.h"
#include <string.h>
#include <cassert>
#include <array>
#include <vector>
#include "lattice/lattice.h"
#include "algorithms/lattice_region_packer.h"
#include <logs/logs.h>

class GhostInitPacker : public comm::Packer<Lattice> {
public:
  typedef Lattice buffer_data_type;

  /**
   *
   * \param p_domain pointer to domain.
   * \param lats_list lattice list.
   */
  explicit GhostInitPacker(const comm::ColoredDomain *p_domain, LatticesList *lats_list);

  const unsigned long sendLength(const int dimension, const int direction) override;

  int sendLength2(const int dimension, const int direction);

  void onSend(buffer_data_type buffer[], const unsigned long send_len, const int dimension,
              const int direction) override;

  void onSend2(buffer_data_type buffer[], int send_len, const int dimension,
              const int direction);

  void onReceive(buffer_data_type buffer[], const unsigned long receive_len, const int dimension,
                 const int direction) override;
 
  void onReceive2(buffer_data_type buffer[], const int receive_len, const int dimension, const int direction);

  static inline MPI_Datatype getMPI_DataType() { return mpi_types::_mpi_type_lattice_data; }
public:
  const comm::ColoredDomain *p_domain;
  LatticesList *lats;
};

#endif // MISA_KMC_GHOST_INIT_PACKER_H
